#include "fisher_kolmogorov_problem.hpp"

#include <algorithm>
#include <cmath>
#include <iostream>
#include <stdexcept>
#include <string>
#include <vector>

namespace {

double solve_uniform_logistic(
	femg::fisher_kolmogorov_problem::TimeScheme scheme,
	double time_step,
	const std::string &graph_file) {
	const double final_time = 1.0;
	const double alpha = 0.2;
	const double initial_value = 0.1;
	char executable[] = "test_fisher_kolmogorov_logistic";
	char *arguments[] = {
		executable,
		const_cast<char *>(graph_file.c_str())
	};

	femg::fisher_kolmogorov_problem problem(final_time, time_step);
	problem.init(2, arguments);
	problem.set_edge_diffusion_coefficients(
		std::vector<double>(problem.number_of_edges(), 1.0));
	problem.set_vertex_reaction_coefficients(
		std::vector<double>(problem.number_of_vertices(), alpha));
	problem.set_vertex_initial_condition(
		std::vector<double>(problem.number_of_vertices(), initial_value));
	problem.set_time_scheme(scheme);
	problem.set_output_enabled(false);
	problem.set_verbose(false);
	problem.set_coefficients();
	problem.assemble_matrices();
	problem.solve();

	const double exponential = std::exp(alpha * final_time);
	const double exact = initial_value * exponential
		/ (1.0 + initial_value * (exponential - 1.0));
	double maximum_error = 0.0;
	for (std::size_t vertex = 0;
		vertex < problem.number_of_vertices(); ++vertex) {
		maximum_error = std::max(
			maximum_error,
			std::abs(problem.vertex_value(vertex) - exact));
	}
	return maximum_error;
}

void verify_scheme(
	const std::string &name,
	femg::fisher_kolmogorov_problem::TimeScheme scheme,
	double minimum_rate,
	double maximum_rate,
	const std::string &graph_file) {
	const std::vector<double> time_steps = {0.02, 0.01, 0.005};
	std::vector<double> errors;
	std::cout << name << "\n";
	for (const double time_step : time_steps) {
		const double error =
			solve_uniform_logistic(scheme, time_step, graph_file);
		errors.push_back(error);
		std::cout << "  dt = " << time_step << ", error = " << error;
		if (errors.size() > 1) {
			const double rate = std::log(
				errors[errors.size() - 2] / errors.back()) / std::log(2.0);
			std::cout << ", rate = " << rate;
			if (rate < minimum_rate || rate > maximum_rate) {
				throw std::runtime_error(
					name + " has an unexpected convergence rate.");
			}
		}
		std::cout << "\n";
	}
}

} // namespace

int main(int argc, char *argv[]) {
	try {
		const std::string graph_file =
			argc >= 2 ? argv[1] : "data/star_4.txt";
		verify_scheme(
			"Corti semi-implicit scheme",
			femg::fisher_kolmogorov_problem::TimeScheme::corti_semi_implicit,
			1.8,
			2.2,
			graph_file);
		verify_scheme(
			"Backward Euler with Newton",
			femg::fisher_kolmogorov_problem::TimeScheme::backward_euler,
			0.9,
			1.1,
			graph_file);
		return 0;
	} catch (const std::exception &error) {
		std::cerr << "Error: " << error.what() << "\n";
		return 1;
	}
}
