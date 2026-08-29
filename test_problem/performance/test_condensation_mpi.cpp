// MPI condensation with edge partitioning and a replicated interface.
//
// The connectome is not vertex-partitioned: every rank owns a contiguous
// block of edges, balanced by interior unknowns, and with it the interior
// DoFs of those edges. The 83 original vertices are the replicated
// interface. Per time step every rank condenses its own edges into a
// partial vertex system, two MPI_Allreduce calls sum the 83x83 Schur
// complement and the condensed right-hand side (a constant payload of
// about 55 kB),
// every rank solves the small interface system redundantly and
// back-substitutes only its own interiors. The mathematics is the one of
// the sequential stepper; only the summation order of the vertex system
// depends on the rank count.
//
// With --validate rank 0 also advances the full sequential stepper and the
// distributed state is gathered and compared at every step. In benchmark
// mode the per-step phase medians (local condensation, Allreduce,
// redundant interface solve, back-substitution) are reduced with MPI_MAX
// to rank 0, which prints one CSV record.

#include "fisher_kolmogorov_problem.hpp"

#include "condensed_stepper.hpp"

#include <mpi.h>
#ifdef _OPENMP
#include <omp.h>
#endif

#include <algorithm>
#include <cstddef>
#include <fstream>
#include <iostream>
#include <limits>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>

namespace {

using femg::Vector;

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

void extrapolate(std::size_t step, const std::vector<double> &state,
	const std::vector<double> &previous, std::vector<double> &result) {
	if (step == 1) {
		result = state;
		return;
	}
	for (std::size_t k = 0; k < state.size(); ++k) {
		result[k] = 1.5 * state[k] - 0.5 * previous[k];
	}
}

// Contiguous edge blocks balanced by interior unknowns: the edge whose
// cumulative interior count starts at p belongs to rank
// floor(p ranks / total), which is non-decreasing along the edge order,
// so every rank receives one contiguous block and the blocks cover all
// edges exactly once.
void partition_edges(const std::vector<condensed::EdgeData> &edges,
	int ranks, int rank, std::size_t &first, std::size_t &last) {
	std::size_t total = 0;
	for (const condensed::EdgeData &edge : edges) {
		total += static_cast<std::size_t>(edge.cells - 1);
	}
	first = edges.size();
	last = 0;
	std::size_t prefix = 0;
	for (std::size_t e = 0; e < edges.size(); ++e) {
		const int owner = total == 0
			? static_cast<int>((e * static_cast<std::size_t>(ranks))
				/ edges.size())
			: std::min(ranks - 1, static_cast<int>(
				(prefix * static_cast<std::size_t>(ranks)) / total));
		if (owner == rank) {
			first = std::min(first, e);
			last = e + 1;
		}
		prefix += static_cast<std::size_t>(edges[e].cells - 1);
	}
	if (first > last) {
		first = 0;
		last = 0;
	}
}

} // namespace

int main(int argc, char *argv[]) {
	MPI_Init(&argc, &argv);
	int rank = 0;
	int ranks = 1;
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	MPI_Comm_size(MPI_COMM_WORLD, &ranks);
	int status = 0;
	try {
		if (argc < 3) {
			if (rank == 0) {
				std::cerr << "usage: " << argv[0]
					<< " <graph_fem.txt> <benchmark|--validate>"
					<< " [max_steps]\n";
			}
			MPI_Finalize();
			return 1;
		}
		std::string graph_file = argv[1];
		const std::string mode = argv[2];
		const std::size_t max_steps =
			argc >= 4 ? std::stoul(argv[3]) : 100;
		int threads = 1;
#ifdef _OPENMP
		threads = argc >= 5 ? std::stoi(argv[4]) : 1;
		omp_set_num_threads(threads);
#endif
		(void)threads;
		const double time_step = 0.2;

		const auto field = read_regional_field(
			"benchmarks/21_fisher_kolmogorov_corti83/results"
			"/reaction_coefficients.csv");
		const auto diffusion = read_normalized_edge_weights(
			"data/connectome/fornari83/edges.csv");

		femg::fisher_kolmogorov_problem problem(
			time_step * static_cast<double>(max_steps), time_step);
		char *graph_argv[] = {argv[0], graph_file.data()};
		problem.init(2, graph_argv);
		problem.set_edge_diffusion_coefficients(diffusion);
		problem.set_vertex_reaction_coefficients(field.alpha);
		problem.set_vertex_initial_condition(field.initial);
		problem.set_output_enabled(false);
		problem.set_verbose(false);
		problem.set_coefficients();

		std::size_t n_vertices = 0;
		const auto all_edges = read_edges(graph_file, diffusion, n_vertices);
		const Vector &alpha = problem.reaction_coefficients();
		const Vector &u0 = problem.solution();
		const std::size_t global_interiors =
			static_cast<std::size_t>(alpha.size()) - n_vertices;

		std::size_t first_edge = 0;
		std::size_t last_edge = 0;
		partition_edges(all_edges, ranks, rank, first_edge, last_edge);
		const std::size_t g0 = first_edge < last_edge
			? all_edges[first_edge].offset : 0;
		std::vector<condensed::EdgeData> local_edges(
			all_edges.begin() + static_cast<std::ptrdiff_t>(first_edge),
			all_edges.begin() + static_cast<std::ptrdiff_t>(last_edge));
		std::size_t local_interiors = 0;
		for (condensed::EdgeData &edge : local_edges) {
			edge.offset -= g0;
			local_interiors += static_cast<std::size_t>(edge.cells - 1);
		}

		// Local layout: [own interiors][all vertices].
		auto localize = [&](const double *global) {
			std::vector<double> local(local_interiors + n_vertices);
			for (std::size_t k = 0; k < local_interiors; ++k) {
				local[k] = global[g0 + k];
			}
			for (std::size_t v = 0; v < n_vertices; ++v) {
				local[local_interiors + v] = global[global_interiors + v];
			}
			return local;
		};
		condensed::Stepper stepper(local_edges, n_vertices,
			localize(alpha.data()), time_step);
		const std::vector<double> initial = localize(u0.data());

		// Gather sizes for the validation assembly on rank 0.
		std::vector<int> counts(static_cast<std::size_t>(ranks));
		std::vector<int> offsets(static_cast<std::size_t>(ranks));
		{
			const int mine = static_cast<int>(local_interiors);
			MPI_Gather(&mine, 1, MPI_INT, counts.data(), 1, MPI_INT, 0,
				MPI_COMM_WORLD);
			int running = 0;
			for (std::size_t r = 0; r < counts.size(); ++r) {
				offsets[r] = running;
				running += counts[r];
			}
		}

		auto do_step = [&](std::size_t step, std::vector<double> &state,
			std::vector<double> &previous, std::vector<double> &chat,
			double *phase_times) {
			extrapolate(step, state, previous, chat);
			previous = state;
			double t0 = MPI_Wtime();
#ifdef _OPENMP
			if (threads > 1) {
				stepper.condense_omp(chat, state);
			} else {
				stepper.condense(chat, state);
			}
#else
			stepper.condense(chat, state);
#endif
			double t1 = MPI_Wtime();
			MPI_Allreduce(MPI_IN_PLACE, stepper.schur().data(),
				static_cast<int>(n_vertices * n_vertices), MPI_DOUBLE,
				MPI_SUM, MPI_COMM_WORLD);
			MPI_Allreduce(MPI_IN_PLACE, stepper.condensed_rhs().data(),
				static_cast<int>(n_vertices), MPI_DOUBLE, MPI_SUM,
				MPI_COMM_WORLD);
			double t2 = MPI_Wtime();
			const Eigen::VectorXd vertices = stepper.solve_interface();
			stepper.scatter_vertices(vertices, state);
			double t3 = MPI_Wtime();
#ifdef _OPENMP
			if (threads > 1) {
				stepper.back_substitute_omp(vertices, state);
			} else {
				stepper.back_substitute(vertices, state);
			}
#else
			stepper.back_substitute(vertices, state);
#endif
			double t4 = MPI_Wtime();
			phase_times[0] = t1 - t0;
			phase_times[1] = t2 - t1;
			phase_times[2] = t3 - t2;
			phase_times[3] = t4 - t3;
		};

		if (mode == "--validate") {
			condensed::Stepper reference_stepper(all_edges, n_vertices,
				std::vector<double>(alpha.data(),
					alpha.data() + alpha.size()), time_step);
			std::vector<double> reference(u0.data(),
				u0.data() + u0.size());
			std::vector<double> reference_prev = reference;
			std::vector<double> reference_chat(reference.size());

			std::vector<double> state = initial;
			std::vector<double> previous = initial;
			std::vector<double> chat(initial.size());
			std::vector<double> gathered(global_interiors + n_vertices);
			double difference = 0.0;
			double phases[4];
			for (std::size_t step = 1; step <= max_steps; ++step) {
				do_step(step, state, previous, chat, phases);
				MPI_Gatherv(state.data(),
					static_cast<int>(local_interiors), MPI_DOUBLE,
					gathered.data(), counts.data(), offsets.data(),
					MPI_DOUBLE, 0, MPI_COMM_WORLD);
				if (rank == 0) {
					for (std::size_t v = 0; v < n_vertices; ++v) {
						gathered[global_interiors + v] =
							state[local_interiors + v];
					}
					extrapolate(step, reference, reference_prev,
						reference_chat);
					reference_prev = reference;
					reference_stepper.step(reference_chat, reference);
					for (std::size_t k = 0; k < gathered.size(); ++k) {
						difference = std::max(difference,
							std::abs(gathered[k] - reference[k]));
					}
				}
			}
			if (rank == 0) {
#ifdef _OPENMP
				std::cout << "VALIDATION,"
					<< global_interiors + n_vertices << "," << ranks
					<< "," << threads << "," << max_steps << ","
					<< difference << "\n";
#else
				std::cout << "VALIDATION,"
					<< global_interiors + n_vertices << "," << ranks
					<< "," << max_steps << "," << difference << "\n";
#endif
			}
			MPI_Finalize();
			return 0;
		}

		std::vector<double> state = initial;
		std::vector<double> previous = initial;
		std::vector<double> chat(initial.size());
		std::vector<double> local_times;
		std::vector<double> comm_times;
		std::vector<double> interface_times;
		std::vector<double> back_times;
		const double loop_start = MPI_Wtime();
		double phases[4];
		for (std::size_t step = 1; step <= max_steps; ++step) {
			do_step(step, state, previous, chat, phases);
			local_times.push_back(phases[0]);
			comm_times.push_back(phases[1]);
			interface_times.push_back(phases[2]);
			back_times.push_back(phases[3]);
		}
		const double loop_total = MPI_Wtime() - loop_start;

		double medians[4] = {median(local_times), median(comm_times),
			median(interface_times), median(back_times)};
		double reduced[4];
		MPI_Reduce(medians, reduced, 4, MPI_DOUBLE, MPI_MAX, 0,
			MPI_COMM_WORLD);
		double total_reduced = 0.0;
		MPI_Reduce(&loop_total, &total_reduced, 1, MPI_DOUBLE, MPI_MAX, 0,
			MPI_COMM_WORLD);
		if (rank == 0) {
#ifdef _OPENMP
			std::cout << "HYBRID," << global_interiors + n_vertices
				<< "," << ranks << "," << threads << "," << max_steps
				<< "," << reduced[0] << "," << reduced[1] << ","
				<< reduced[2] << "," << reduced[3] << ","
				<< reduced[0] + reduced[1] + reduced[2] + reduced[3]
				<< "," << total_reduced << "\n";
#else
			std::cout << "MPI," << global_interiors + n_vertices << ","
				<< ranks << "," << max_steps << ","
				<< reduced[0] << "," << reduced[1] << "," << reduced[2]
				<< "," << reduced[3] << ","
				<< reduced[0] + reduced[1] + reduced[2] + reduced[3]
				<< "," << total_reduced << "\n";
#endif
		}
	} catch (const std::exception &error) {
		std::cerr << "Error on rank " << rank << ": " << error.what()
			<< "\n";
		status = 1;
	}
	MPI_Finalize();
	return status;
}
