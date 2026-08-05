#include "parabolic_problem.hpp"
#include "star_sine_functions.hpp"

#include <cmath>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <memory>
#include <string>
#include <vector>

int main() {
	const double diffusion = 1.0;
	const double T = 1.0;
	const double theta = 1.0;

	const std::vector<int> n_cells_values = {10, 20, 40, 80};

	std::filesystem::create_directories("output/convergence");
	std::ofstream convergence_file(
		"output/convergence/star_sine_space_convergence.csv");
	convergence_file << "h,L2,rate\n";

	double previous_error = 0.0;
	double previous_h = 0.0;

	std::cout << std::setw(12) << "h"
						<< std::setw(18) << "L2 error"
						<< std::setw(12) << "rate" << "\n";

	for (const int n_cells : n_cells_values) {
		const double h = 1.0 / static_cast<double>(n_cells);
		const double deltat = 0.025 * h * h;

		std::string graph_file =
			"data/star_4_" + std::to_string(n_cells) + ".txt";
		std::string executable_name = "test_heat_star_sine_convergence";
		char *argv[] = {executable_name.data(), graph_file.data()};

		femg::parabolic_problem problem(diffusion, T, deltat, theta);
		problem.set_initial_condition(std::make_shared<StarSineInitialCondition>());
		problem.set_source_function(std::make_shared<StarSineForcingTerm>());
		problem.set_dirichlet_vertices(
			{0, 1, 2, 3, 4},
			[](std::size_t, double) {
				return 0.0;
			});
		problem.set_output_enabled(false);
		problem.set_verbose(false);

		problem.init(2, argv);
		problem.set_coefficients();
		problem.assemble_matrices();
		problem.solve();

		const StarSineExactSolution exact_solution;
		const double error_l2 = problem.compute_l2_error(exact_solution);

		double rate = 0.0;
		if (previous_error > 0.0) {
			rate = std::log(previous_error / error_l2)
				/ std::log(previous_h / h);
		}

		std::cout << std::setw(12) << h
							<< std::setw(18) << error_l2;
		if (previous_error > 0.0) {
			std::cout << std::setw(12) << rate;
		} else {
			std::cout << std::setw(12) << "-";
		}
		std::cout << "\n";

		convergence_file << h << "," << error_l2 << ",";
		if (previous_error > 0.0) {
			convergence_file << rate;
		}
		convergence_file << "\n";

		previous_error = error_l2;
		previous_h = h;
	}

	return 0;
}
