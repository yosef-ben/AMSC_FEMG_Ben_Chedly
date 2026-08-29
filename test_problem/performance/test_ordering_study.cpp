// Sequential ordering study on the deterministic Corti-83 reference problem.
//
// The production solver of fisher_kolmogorov_problem factorizes the
// semi-implicit time-step matrix
//
//     K = M + (dt/2) H - (dt/2) W(c_extrapolated)
//
// once per time step with Eigen::SparseLU, whose default fill-reducing
// ordering is COLAMD. This driver replays exactly that time loop through the
// public accessors of the class, without touching the library, and times the
// three phases of every step (matrix rebuild, factorization, triangular
// solves) for one solver-ordering variant at a time:
//
//     lu_colamd     Eigen::SparseLU + COLAMD, the production configuration
//     lu_natural    Eigen::SparseLU + the natural (interiors-first) ordering
//     lu_amd        Eigen::SparseLU + AMD
//     ldlt_amd      Eigen::SimplicialLDLT + AMD (the matrix is symmetric)
//     ldlt_natural  Eigen::SimplicialLDLT + the natural ordering
//     lu_rcm        reverse Cuthill-McKee permutation applied to the matrix,
//                   then Eigen::SparseLU with the natural ordering
//
// For each variant the driver reports the bandwidth and profile of the matrix
// handed to the solver, the number of nonzeros of the computed factors, the
// median per-step times of the three phases, the median of a numeric-only
// refactorization (the production loop repeats the symbolic analysis at every
// step because the reaction weight changes; the pattern does not), and the
// maximum difference of the final state against the lu_colamd trajectory run
// for the same number of steps. With --validate the driver instead checks
// that its replayed lu_colamd loop reproduces the solve() of the library
// class to round-off.
//
// The model is the deterministic regional-rate problem of benchmark 21, the
// reference problem of the sequential-performance record: consistent mass,
// semi-implicit scheme, dt = 0.2 years over 20 years, D_e = w_e / max(w),
// the seven regional conversion rates read from the stored
// reaction_coefficients.csv and the entorhinal seed 0.10 over the 0.01
// background.

#include "fisher_kolmogorov_problem.hpp"

#include <Eigen/OrderingMethods>
#include <Eigen/SparseCholesky>
#include <Eigen/SparseLU>

#include <algorithm>
#include <chrono>
#include <cmath>
#include <cstddef>
#include <fstream>
#include <iostream>
#include <limits>
#include <queue>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>

namespace {

using femg::SparseMatrix;
using femg::Vector;
using clock_type = std::chrono::steady_clock;

double seconds_since(const clock_type::time_point &start) {
	return std::chrono::duration<double>(clock_type::now() - start).count();
}

std::vector<std::string> split_csv_line(const std::string &line) {
	std::vector<std::string> fields;
	std::stringstream stream(line);
	std::string field;
	while (std::getline(stream, field, ',')) {
		if (!field.empty() && field.back() == '\r') {
			field.pop_back();
		}
		fields.push_back(field);
	}
	return fields;
}

bool contains(const std::string &text, const std::string &token) {
	return text.find(token) != std::string::npos;
}

std::string lower_case(std::string text) {
	std::transform(text.begin(), text.end(), text.begin(),
		[](unsigned char character) {
			return static_cast<char>(std::tolower(character));
		});
	return text;
}

// node_id,name,region,alpha of the stored benchmark-21 run: the exact
// regional rate field of the reference problem, no classification redone.
struct RegionalField {
	std::vector<double> alpha;
	std::vector<double> initial;
};

RegionalField read_regional_field(const std::string &path) {
	std::ifstream input(path);
	if (!input) {
		throw std::runtime_error("Unable to open " + path);
	}
	std::string line;
	std::getline(input, line);
	RegionalField field;
	while (std::getline(input, line)) {
		const auto fields = split_csv_line(line);
		if (fields.size() != 4) {
			throw std::runtime_error("Unexpected coefficient CSV row.");
		}
		const auto id = static_cast<std::size_t>(std::stoul(fields[0]));
		if (id != field.alpha.size()) {
			throw std::runtime_error("Coefficient IDs are not contiguous.");
		}
		field.alpha.push_back(std::stod(fields[3]));
		field.initial.push_back(
			contains(lower_case(fields[1]), "entorhinal") ? 0.10 : 0.01);
	}
	return field;
}

std::vector<double> read_normalized_edge_weights(const std::string &path) {
	std::ifstream input(path);
	if (!input) {
		throw std::runtime_error("Unable to open " + path);
	}
	std::string line;
	std::getline(input, line);
	std::vector<double> weights;
	while (std::getline(input, line)) {
		const auto fields = split_csv_line(line);
		if (fields.size() != 7) {
			throw std::runtime_error("Unexpected edge CSV row.");
		}
		weights.push_back(std::stod(fields[6]));
	}
	const double maximum = *std::max_element(weights.begin(), weights.end());
	for (double &weight : weights) {
		weight /= maximum;
	}
	return weights;
}

// Lower-triangle bandwidth and profile (envelope) of a compressed matrix.
void bandwidth_profile(const SparseMatrix &matrix,
	long &bandwidth, long &profile) {
	bandwidth = 0;
	profile = 0;
	std::vector<long> reach(static_cast<std::size_t>(matrix.rows()), 0);
	for (int column = 0; column < matrix.outerSize(); ++column) {
		for (SparseMatrix::InnerIterator entry(matrix, column); entry;
			++entry) {
			const long distance =
				static_cast<long>(entry.row()) - static_cast<long>(entry.col());
			if (distance > 0) {
				auto &row_reach =
					reach[static_cast<std::size_t>(entry.row())];
				row_reach = std::max(row_reach, distance);
			}
		}
	}
	for (const long row_reach : reach) {
		bandwidth = std::max(bandwidth, row_reach);
		profile += row_reach;
	}
}

// Reverse Cuthill-McKee on the pattern of a symmetric matrix: breadth-first
// search from a minimum-degree vertex, neighbours visited by increasing
// degree, and the resulting order reversed.
Eigen::PermutationMatrix<Eigen::Dynamic, Eigen::Dynamic, int>
reverse_cuthill_mckee(const SparseMatrix &matrix) {
	const auto n = static_cast<std::size_t>(matrix.rows());
	std::vector<std::vector<int>> adjacency(n);
	for (int column = 0; column < matrix.outerSize(); ++column) {
		for (SparseMatrix::InnerIterator entry(matrix, column); entry;
			++entry) {
			if (entry.row() != entry.col()) {
				adjacency[static_cast<std::size_t>(entry.col())]
					.push_back(static_cast<int>(entry.row()));
			}
		}
	}
	std::vector<int> order;
	order.reserve(n);
	std::vector<bool> visited(n, false);
	auto degree = [&adjacency](int vertex) {
		return adjacency[static_cast<std::size_t>(vertex)].size();
	};
	for (;;) {
		int start = -1;
		for (std::size_t vertex = 0; vertex < n; ++vertex) {
			if (!visited[vertex]
				&& (start < 0
					|| degree(static_cast<int>(vertex))
						< degree(start))) {
				start = static_cast<int>(vertex);
			}
		}
		if (start < 0) {
			break;
		}
		std::queue<int> frontier;
		frontier.push(start);
		visited[static_cast<std::size_t>(start)] = true;
		while (!frontier.empty()) {
			const int vertex = frontier.front();
			frontier.pop();
			order.push_back(vertex);
			auto neighbours = adjacency[static_cast<std::size_t>(vertex)];
			std::sort(neighbours.begin(), neighbours.end(),
				[&degree](int a, int b) { return degree(a) < degree(b); });
			for (const int neighbour : neighbours) {
				if (!visited[static_cast<std::size_t>(neighbour)]) {
					visited[static_cast<std::size_t>(neighbour)] = true;
					frontier.push(neighbour);
				}
			}
		}
	}
	std::reverse(order.begin(), order.end());
	Eigen::PermutationMatrix<Eigen::Dynamic, Eigen::Dynamic, int>
		permutation(static_cast<int>(n));
	for (std::size_t position = 0; position < n; ++position) {
		// New position of old index order[position].
		permutation.indices()[order[position]] =
			static_cast<int>(position);
	}
	return permutation;
}

struct StepTimes {
	std::vector<double> rebuild;
	std::vector<double> factor;
	std::vector<double> refactor;
	std::vector<double> solve;
	long factor_nnz = 0;
	long bandwidth = 0;
	long profile = 0;
	std::size_t steps = 0;
	bool failed = false;
	Vector final_state;
};

// Factor nonzeros: SparseLU exposes nnzL/nnzU, SimplicialLDLT exposes the
// lower factor; resolved by overload at compile time (C++17, no concepts).
template <class Solver>
auto factor_nonzeros(const Solver &solver, long, int)
	-> decltype(solver.nnzL() + solver.nnzU()) {
	return solver.nnzL() + solver.nnzU();
}

template <class Solver>
long factor_nonzeros(const Solver &solver, long rows, long) {
	return solver.matrixL().nestedExpression().nonZeros() + rows;
}

double median(std::vector<double> values) {
	if (values.empty()) {
		return std::numeric_limits<double>::quiet_NaN();
	}
	std::sort(values.begin(), values.end());
	return values[values.size() / 2];
}

double total(const std::vector<double> &values) {
	double sum = 0.0;
	for (const double value : values) {
		sum += value;
	}
	return sum;
}

// One semi-implicit time loop, identical to the solve() of the library
// class, with one solver object per step exactly as the production code
// calls compute() on the rebuilt matrix. Solver is any Eigen sparse solver;
// Permute optionally twists the system by a fixed permutation (lu_rcm).
template <class Solver>
StepTimes run_variant(const femg::fisher_kolmogorov_problem &problem,
	double time_step, std::size_t max_steps, double step_time_limit,
	const Eigen::PermutationMatrix<Eigen::Dynamic, Eigen::Dynamic, int>
		*permutation) {
	const SparseMatrix &M = problem.mass_matrix();
	const SparseMatrix &H = problem.diffusion_matrix();
	StepTimes times;
	Vector current = problem.solution();
	Vector previous = current;
	bool metrics_recorded = false;

	for (std::size_t step = 1; step <= max_steps; ++step) {
		const Vector extrapolated = (step == 1)
			? current
			: (1.5 * current - 0.5 * previous).eval();

		auto start = clock_type::now();
		const SparseMatrix reaction =
			problem.assemble_reaction_weight_matrix(extrapolated, 1.0);
		SparseMatrix lhs = M + 0.5 * time_step * H
			- 0.5 * time_step * reaction;
		SparseMatrix rhs = M - 0.5 * time_step * H
			+ 0.5 * time_step * reaction;
		lhs.makeCompressed();
		rhs.makeCompressed();
		Vector rhs_vector = rhs * current;
		if (permutation != nullptr) {
			lhs = lhs.twistedBy(*permutation);
			lhs.makeCompressed();
			rhs_vector = (*permutation) * rhs_vector;
		}
		times.rebuild.push_back(seconds_since(start));

		if (!metrics_recorded) {
			bandwidth_profile(lhs, times.bandwidth, times.profile);
		}

		start = clock_type::now();
		Solver solver;
		solver.compute(lhs);
		times.factor.push_back(seconds_since(start));
		if (solver.info() != Eigen::Success) {
			times.failed = true;
			return times;
		}
		if (!metrics_recorded) {
			times.factor_nnz = factor_nonzeros(solver, lhs.rows(), 0);
			// The pattern of K is fixed across the steps, so a solver that
			// keeps its symbolic analysis pays only the numeric phase.
			Solver refactor_solver;
			refactor_solver.analyzePattern(lhs);
			start = clock_type::now();
			refactor_solver.factorize(lhs);
			times.refactor.push_back(seconds_since(start));
			metrics_recorded = true;
		}

		start = clock_type::now();
		Vector next = solver.solve(rhs_vector);
		times.solve.push_back(seconds_since(start));
		if (solver.info() != Eigen::Success) {
			times.failed = true;
			return times;
		}
		if (permutation != nullptr) {
			next = permutation->inverse() * next;
		}

		previous = current;
		current = next;
		times.steps = step;
		if (times.factor.back() + times.solve.back() > step_time_limit
			&& step >= 5) {
			break;
		}
	}
	times.final_state = current;
	return times;
}

using LuColamd = Eigen::SparseLU<SparseMatrix, Eigen::COLAMDOrdering<int>>;
using LuNatural = Eigen::SparseLU<SparseMatrix, Eigen::NaturalOrdering<int>>;
using LuAmd = Eigen::SparseLU<SparseMatrix, Eigen::AMDOrdering<int>>;
using LdltAmd = Eigen::SimplicialLDLT<SparseMatrix, Eigen::Lower,
	Eigen::AMDOrdering<int>>;
using LdltNatural = Eigen::SimplicialLDLT<SparseMatrix, Eigen::Lower,
	Eigen::NaturalOrdering<int>>;

} // namespace

int main(int argc, char *argv[]) {
	try {
		if (argc < 3) {
			std::cerr << "usage: " << argv[0]
				<< " <graph_fem.txt> <variant|--validate>"
				<< " [max_steps] [step_time_limit_s]\n";
			return 1;
		}
		std::string graph_file = argv[1];
		const std::string variant = argv[2];
		const std::size_t max_steps =
			argc >= 4 ? std::stoul(argv[3]) : 100;
		const double step_time_limit =
			argc >= 5 ? std::stod(argv[4]) : 5.0;
		const double time_step = 0.2;
		const double final_time = time_step * static_cast<double>(max_steps);

		const auto field = read_regional_field(
			"benchmarks/21_fisher_kolmogorov_corti83/results"
			"/reaction_coefficients.csv");
		const auto diffusion = read_normalized_edge_weights(
			"data/connectome/fornari83/edges.csv");

		femg::fisher_kolmogorov_problem problem(final_time, time_step);
		char *graph_argv[] = {argv[0], graph_file.data()};
		problem.init(2, graph_argv);
		problem.set_edge_diffusion_coefficients(diffusion);
		problem.set_vertex_reaction_coefficients(field.alpha);
		problem.set_vertex_initial_condition(field.initial);
		problem.set_time_scheme(
			femg::fisher_kolmogorov_problem::TimeScheme::corti_semi_implicit);
		problem.set_output_enabled(false);
		problem.set_verbose(false);
		problem.set_coefficients();
		problem.assemble_matrices();

		if (variant == "--validate") {
			// The replayed loop must reproduce the library class to
			// round-off: same matrices, same expressions, same solver.
			const StepTimes replayed = run_variant<LuColamd>(
				problem, time_step, max_steps,
				std::numeric_limits<double>::infinity(), nullptr);
			femg::fisher_kolmogorov_problem reference(final_time, time_step);
			reference.init(2, graph_argv);
			reference.set_edge_diffusion_coefficients(diffusion);
			reference.set_vertex_reaction_coefficients(field.alpha);
			reference.set_vertex_initial_condition(field.initial);
			reference.set_time_scheme(femg::fisher_kolmogorov_problem
				::TimeScheme::corti_semi_implicit);
			reference.set_output_enabled(false);
			reference.set_verbose(false);
			reference.set_coefficients();
			reference.assemble_matrices();
			reference.solve();
			const double difference =
				(replayed.final_state - reference.solution())
					.cwiseAbs().maxCoeff();
			std::cout << "VALIDATION," << problem.number_of_dofs() << ","
				<< max_steps << "," << difference << "\n";
			return 0;
		}

		// The reference trajectory for the numerical comparison: the
		// production configuration under the same adaptive step limit.
		const StepTimes reference = run_variant<LuColamd>(problem,
			time_step, max_steps, step_time_limit, nullptr);

		StepTimes times;
		if (variant == "lu_colamd") {
			times = reference;
		} else if (variant == "lu_natural") {
			times = run_variant<LuNatural>(problem, time_step, max_steps,
				step_time_limit, nullptr);
		} else if (variant == "lu_amd") {
			times = run_variant<LuAmd>(problem, time_step, max_steps,
				step_time_limit, nullptr);
		} else if (variant == "ldlt_amd") {
			times = run_variant<LdltAmd>(problem, time_step, max_steps,
				step_time_limit, nullptr);
		} else if (variant == "ldlt_natural") {
			times = run_variant<LdltNatural>(problem, time_step, max_steps,
				step_time_limit, nullptr);
		} else if (variant == "lu_rcm") {
			const SparseMatrix &M = problem.mass_matrix();
			const SparseMatrix &H = problem.diffusion_matrix();
			SparseMatrix pattern = M + H;
			pattern.makeCompressed();
			const auto permutation = reverse_cuthill_mckee(pattern);
			times = run_variant<LuNatural>(problem, time_step, max_steps,
				step_time_limit, &permutation);
		} else {
			throw std::invalid_argument("Unknown variant: " + variant);
		}

		double difference = std::numeric_limits<double>::quiet_NaN();
		if (!times.failed && times.steps > 0) {
			// Compare the states at the last step both runs reached.
			const std::size_t shared_steps =
				std::min(times.steps, reference.steps);
			if (times.steps == reference.steps) {
				difference = (times.final_state - reference.final_state)
					.cwiseAbs().maxCoeff();
			} else {
				const StepTimes shortened = run_variant<LuColamd>(problem,
					time_step, shared_steps,
					std::numeric_limits<double>::infinity(), nullptr);
				const Vector &variant_state =
					times.steps == shared_steps
						? times.final_state
						: run_variant<LuColamd>(problem, time_step,
							shared_steps,
							std::numeric_limits<double>::infinity(),
							nullptr).final_state;
				difference = (variant_state - shortened.final_state)
					.cwiseAbs().maxCoeff();
			}
		}

		SparseMatrix k_pattern =
			problem.mass_matrix() + problem.diffusion_matrix();
		k_pattern.makeCompressed();
		std::cout << "ORDERING," << variant << ","
			<< problem.number_of_dofs() << ","
			<< k_pattern.nonZeros() << ","
			<< times.bandwidth << "," << times.profile << ","
			<< times.factor_nnz << ","
			<< times.steps << "," << (times.failed ? 1 : 0) << ","
			<< median(times.rebuild) << ","
			<< median(times.factor) << ","
			<< median(times.refactor) << ","
			<< median(times.solve) << ","
			<< median(times.rebuild) + median(times.factor)
				+ median(times.solve) << ","
			<< total(times.rebuild) + total(times.factor)
				+ total(times.solve) << ","
			<< difference << "\n";
		return 0;
	} catch (const std::exception &error) {
		std::cerr << "Error: " << error.what() << "\n";
		return 1;
	}
}
