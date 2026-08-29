// Sequential static condensation against the full-system reference.
//
// The driver advances the deterministic Corti-83 reference problem with two
// interchangeable engines and never touches the library:
//
//     condensed     the per-edge static condensation of condensed_stepper.hpp
//     ldlt_natural  the full-system loop of the ordering study: rebuild of
//                   the global matrices through the public accessors of the
//                   class and Eigen::SimplicialLDLT on the natural ordering,
//                   the sequential reference selected by benchmark 28
//
// With --validate the condensed trajectory is compared step by step against
// the solve() of the production class over the whole horizon. In benchmark
// mode both engines run the same steps and one CSV record reports the phase
// medians of each and the final-state difference between the two.

#include "fisher_kolmogorov_problem.hpp"

#include "condensed_stepper.hpp"

#include <Eigen/OrderingMethods>
#include <Eigen/SparseCholesky>

#include <algorithm>
#include <chrono>
#include <cstddef>
#include <fstream>
#include <iostream>
#include <limits>
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

// The graph file of the library: header, optional coordinates, one line per
// edge with endpoints, length and cells. Offsets follow the file order,
// which is the edge-index order of the DoF numbering.
std::vector<condensed::EdgeData> read_edges(const std::string &path,
	const std::vector<double> &diffusion, std::size_t &n_vertices) {
	std::ifstream input(path);
	if (!input) {
		throw std::runtime_error("Unable to open " + path);
	}
	std::size_t n_edges = 0;
	if (!(input >> n_vertices >> n_edges)) {
		throw std::runtime_error("Invalid header in " + path);
	}
	std::vector<double> data;
	double value = 0.0;
	while (input >> value) {
		data.push_back(value);
	}
	std::size_t offset_in_file = 0;
	if (data.size() == 2 * n_vertices + 4 * n_edges) {
		offset_in_file = 2 * n_vertices;
	} else if (data.size() != 4 * n_edges) {
		throw std::runtime_error("Invalid data size in " + path);
	}
	if (diffusion.size() != n_edges) {
		throw std::runtime_error("One diffusivity per edge is required.");
	}
	std::vector<condensed::EdgeData> edges(n_edges);
	std::size_t dof_offset = 0;
	for (std::size_t e = 0; e < n_edges; ++e) {
		condensed::EdgeData &edge = edges[e];
		edge.source = static_cast<int>(data[offset_in_file + 4 * e]);
		edge.target = static_cast<int>(data[offset_in_file + 4 * e + 1]);
		const double length = data[offset_in_file + 4 * e + 2];
		edge.cells = static_cast<int>(data[offset_in_file + 4 * e + 3]);
		edge.h = length / static_cast<double>(edge.cells);
		edge.diffusion = diffusion[e];
		edge.offset = dof_offset;
		dof_offset += static_cast<std::size_t>(edge.cells - 1);
	}
	return edges;
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

std::vector<double> to_std(const Vector &vector) {
	return std::vector<double>(vector.data(), vector.data() + vector.size());
}

// The condensed loop: extrapolation exactly as the production scheme.
struct CondensedRun {
	std::vector<double> local;
	std::vector<double> interface_solve;
	std::vector<double> back;
	std::vector<double> final_state;
	std::vector<std::vector<double>> trajectory;   // filled when recorded
};

CondensedRun run_condensed(condensed::Stepper &stepper,
	const std::vector<double> &initial, std::size_t n_steps,
	bool record_trajectory) {
	CondensedRun run;
	std::vector<double> state = initial;
	std::vector<double> previous = initial;
	std::vector<double> extrapolated(state.size());
	for (std::size_t step = 1; step <= n_steps; ++step) {
		if (step == 1) {
			extrapolated = state;
		} else {
			for (std::size_t k = 0; k < state.size(); ++k) {
				extrapolated[k] = 1.5 * state[k] - 0.5 * previous[k];
			}
		}
		previous = state;
		const condensed::StepTimes times =
			stepper.step(extrapolated, state);
		run.local.push_back(times.local);
		run.interface_solve.push_back(times.interface);
		run.back.push_back(times.back);
		if (record_trajectory) {
			run.trajectory.push_back(state);
		}
	}
	run.final_state = state;
	return run;
}

using LdltNatural = Eigen::SimplicialLDLT<SparseMatrix, Eigen::Lower,
	Eigen::NaturalOrdering<int>>;

struct FullRun {
	std::vector<double> rebuild;
	std::vector<double> factor;
	std::vector<double> solve;
	Vector final_state;
};

FullRun run_full(const femg::fisher_kolmogorov_problem &problem,
	double time_step, std::size_t n_steps) {
	const SparseMatrix &M = problem.mass_matrix();
	const SparseMatrix &H = problem.diffusion_matrix();
	FullRun run;
	Vector current = problem.solution();
	Vector previous = current;
	for (std::size_t step = 1; step <= n_steps; ++step) {
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
		const Vector rhs_vector = rhs * current;
		run.rebuild.push_back(seconds_since(start));

		start = clock_type::now();
		LdltNatural solver;
		solver.compute(lhs);
		run.factor.push_back(seconds_since(start));
		if (solver.info() != Eigen::Success) {
			throw std::runtime_error("LDLT factorization failed.");
		}

		start = clock_type::now();
		const Vector next = solver.solve(rhs_vector);
		run.solve.push_back(seconds_since(start));

		previous = current;
		current = next;
	}
	run.final_state = current;
	return run;
}

} // namespace

int main(int argc, char *argv[]) {
	try {
		if (argc < 3) {
			std::cerr << "usage: " << argv[0]
				<< " <graph_fem.txt> <benchmark|--validate> [max_steps]\n";
			return 1;
		}
		std::string graph_file = argv[1];
		const std::string mode = argv[2];
		const std::size_t max_steps =
			argc >= 4 ? std::stoul(argv[3]) : 100;
		const double time_step = 0.2;
		const double final_time =
			time_step * static_cast<double>(max_steps);

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
		problem.set_time_scheme(femg::fisher_kolmogorov_problem
			::TimeScheme::corti_semi_implicit);
		problem.set_output_enabled(false);
		problem.set_verbose(false);
		problem.set_coefficients();
		problem.assemble_matrices();

		std::size_t n_vertices = 0;
		auto edges = read_edges(graph_file, diffusion, n_vertices);
		condensed::Stepper stepper(std::move(edges), n_vertices,
			to_std(problem.reaction_coefficients()), time_step);
		if (stepper.n_dofs() != problem.number_of_dofs()) {
			throw std::runtime_error("DoF layouts disagree.");
		}
		const std::vector<double> initial = to_std(problem.solution());

		if (mode == "--validate") {
			// The class records its own trajectory through the callback;
			// the condensed loop must match it at every step.
			std::vector<Vector> reference;
			reference.reserve(max_steps);
			problem.set_time_step_callback(
				[&reference](std::size_t step, double, const Vector &u) {
					if (step > 0) {
						reference.push_back(u);
					}
				});
			problem.solve();
			const CondensedRun run = run_condensed(stepper, initial,
				max_steps, true);
			double difference = 0.0;
			for (std::size_t step = 0; step < max_steps; ++step) {
				const Vector &u = reference.at(step);
				const std::vector<double> &v = run.trajectory.at(step);
				for (Eigen::Index k = 0; k < u.size(); ++k) {
					difference = std::max(difference,
						std::abs(u(k) - v[static_cast<std::size_t>(k)]));
				}
			}
			std::cout << "VALIDATION," << problem.number_of_dofs() << ","
				<< max_steps << "," << difference << "\n";
			return 0;
		}

		const CondensedRun run = run_condensed(stepper, initial,
			max_steps, false);
		const FullRun full = run_full(problem, time_step, max_steps);
		double difference = 0.0;
		for (Eigen::Index k = 0; k < full.final_state.size(); ++k) {
			difference = std::max(difference,
				std::abs(full.final_state(k)
					- run.final_state[static_cast<std::size_t>(k)]));
		}

		const double condensed_step = median(run.local)
			+ median(run.interface_solve) + median(run.back);
		const double full_step = median(full.rebuild)
			+ median(full.factor) + median(full.solve);
		std::cout << "CONDENSE," << problem.number_of_dofs() << ","
			<< max_steps << ","
			<< median(run.local) << "," << median(run.interface_solve)
			<< "," << median(run.back) << "," << condensed_step << ","
			<< total(run.local) + total(run.interface_solve)
				+ total(run.back) << ","
			<< median(full.rebuild) << "," << median(full.factor) << ","
			<< median(full.solve) << "," << full_step << ","
			<< total(full.rebuild) + total(full.factor)
				+ total(full.solve) << ","
			<< difference << "\n";
		return 0;
	} catch (const std::exception &error) {
		std::cerr << "Error: " << error.what() << "\n";
		return 1;
	}
}
