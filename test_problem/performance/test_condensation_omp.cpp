// OpenMP condensation against the sequential condensed reference.
//
// Same model, same stepper, two engines: the sequential step() and the
// step_omp() whose edge loops run under OpenMP with thread-local vertex
// accumulators and one critical merge per thread. The mathematics is
// unchanged; only the summation order of the vertex system differs, so the
// two trajectories agree to round-off.
//
// With --validate the two engines advance in lockstep and every step is
// compared in the maximum norm. In benchmark mode both run the full
// horizon in the same process and one CSV record carries the phase medians
// of each, the number of threads and the final-state difference.

#include "fisher_kolmogorov_problem.hpp"

#include "condensed_stepper.hpp"

#include <omp.h>

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

double total(const std::vector<double> &values) {
	double sum = 0.0;
	for (const double value : values) {
		sum += value;
	}
	return sum;
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

} // namespace

int main(int argc, char *argv[]) {
	try {
		if (argc < 4) {
			std::cerr << "usage: " << argv[0]
				<< " <graph_fem.txt> <benchmark|--validate> <threads>"
				<< " [max_steps]\n";
			return 1;
		}
		std::string graph_file = argv[1];
		const std::string mode = argv[2];
		const int threads = std::stoi(argv[3]);
		const std::size_t max_steps =
			argc >= 5 ? std::stoul(argv[4]) : 100;
		const double time_step = 0.2;
		omp_set_num_threads(threads);

		const auto field = read_regional_field(
			"benchmarks/21_fisher_kolmogorov_corti83/results"
			"/reaction_coefficients.csv");
		const auto diffusion = read_normalized_edge_weights(
			"data/connectome/fornari83/edges.csv");

		// The class interpolates the vertex data to the DoFs; its solver
		// is never called here.
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
		auto edges = read_edges(graph_file, diffusion, n_vertices);
		const Vector &alpha = problem.reaction_coefficients();
		condensed::Stepper stepper(std::move(edges), n_vertices,
			std::vector<double>(alpha.data(), alpha.data() + alpha.size()),
			time_step);
		const Vector &u0 = problem.solution();
		const std::vector<double> initial(u0.data(), u0.data() + u0.size());

		if (mode == "--validate") {
			// Lockstep: both engines advance from the same history and
			// every step is compared before moving on.
			std::vector<double> seq = initial;
			std::vector<double> seq_prev = initial;
			std::vector<double> par = initial;
			std::vector<double> par_prev = initial;
			std::vector<double> chat(initial.size());
			double difference = 0.0;
			for (std::size_t step = 1; step <= max_steps; ++step) {
				extrapolate(step, seq, seq_prev, chat);
				seq_prev = seq;
				stepper.step(chat, seq);
				extrapolate(step, par, par_prev, chat);
				par_prev = par;
				stepper.step_omp(chat, par);
				for (std::size_t k = 0; k < seq.size(); ++k) {
					difference = std::max(difference,
						std::abs(seq[k] - par[k]));
				}
			}
			std::cout << "VALIDATION," << stepper.n_dofs() << ","
				<< threads << "," << max_steps << ","
				<< difference << "\n";
			return 0;
		}

		// Sequential reference in the same process, then the OpenMP run.
		std::vector<double> state = initial;
		std::vector<double> previous = initial;
		std::vector<double> chat(initial.size());
		std::vector<double> seq_local;
		std::vector<double> seq_interface;
		std::vector<double> seq_back;
		for (std::size_t step = 1; step <= max_steps; ++step) {
			extrapolate(step, state, previous, chat);
			previous = state;
			const condensed::StepTimes times = stepper.step(chat, state);
			seq_local.push_back(times.local);
			seq_interface.push_back(times.interface);
			seq_back.push_back(times.back);
		}
		const std::vector<double> seq_final = state;

		state = initial;
		previous = initial;
		std::vector<double> par_local;
		std::vector<double> par_reduce;
		std::vector<double> par_interface;
		std::vector<double> par_back;
		for (std::size_t step = 1; step <= max_steps; ++step) {
			extrapolate(step, state, previous, chat);
			previous = state;
			const condensed::Stepper::OmpTimes times =
				stepper.step_omp(chat, state);
			par_local.push_back(times.local);
			par_reduce.push_back(times.reduce);
			par_interface.push_back(times.interface);
			par_back.push_back(times.back);
		}
		double difference = 0.0;
		for (std::size_t k = 0; k < state.size(); ++k) {
			difference = std::max(difference,
				std::abs(state[k] - seq_final[k]));
		}

		const double seq_step = median(seq_local)
			+ median(seq_interface) + median(seq_back);
		const double par_step = median(par_local) + median(par_reduce)
			+ median(par_interface) + median(par_back);
		std::cout << "OMP," << stepper.n_dofs() << "," << threads << ","
			<< max_steps << ","
			<< median(seq_local) << "," << median(seq_interface) << ","
			<< median(seq_back) << "," << seq_step << ","
			<< total(seq_local) + total(seq_interface)
				+ total(seq_back) << ","
			<< median(par_local) << "," << median(par_reduce) << ","
			<< median(par_interface) << "," << median(par_back) << ","
			<< par_step << ","
			<< total(par_local) + total(par_reduce)
				+ total(par_interface) + total(par_back) << ","
			<< difference << "\n";
		return 0;
	} catch (const std::exception &error) {
		std::cerr << "Error: " << error.what() << "\n";
		return 1;
	}
}
