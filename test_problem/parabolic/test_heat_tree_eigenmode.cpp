#include "eigenmode_heat_functions.hpp"
#include "parabolic_problem.hpp"
#include "spectral_problem.hpp"

#include <cmath>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <memory>
#include <string>

int main(int argc, char *argv[]) {
	try {
		const double diffusion = 1.0;
		const double T = 1.0;
		const double deltat = 0.01;
		const double theta = 0.5;
		const std::size_t mode = 5;

		std::string graph_file = "data/tree_15.txt";
		char *local_argv[] = {argv[0], graph_file.data()};
		const int graph_argc = (argc >= 2) ? argc : 2;
		char **graph_argv = (argc >= 2) ? argv : local_argv;

		femg::spectral_problem eigen_problem(10);
		eigen_problem.init(graph_argc, graph_argv);
		eigen_problem.set_coefficients();
		eigen_problem.assembly();
		eigen_problem.solve();

		const double lambda = eigen_problem.eigenvalue(mode);
		EigenmodeHeatData eigenmode_data(
			eigen_problem.edge_nodal_values(mode),
			eigen_problem.edge_lengths(),
			lambda,
			diffusion);

		femg::parabolic_problem problem(diffusion, T, deltat, theta);
		problem.set_initial_condition(
			std::make_shared<EigenmodeHeatInitialCondition>(eigenmode_data));
		problem.set_output_directory("output/visualization/tree_eigenmode");
		problem.set_verbose(false);

		problem.init(graph_argc, graph_argv);
		problem.set_coefficients();
		problem.assemble_matrices();
		problem.print_matrix_summary(std::cout);

		std::cout << "Using eigenmode " << mode
							<< " with lambda = " << lambda << "\n";

		std::filesystem::create_directories("output/convergence");
		std::ofstream decay_out("output/convergence/tree_eigenmode_decay.csv");
		decay_out << "time,numerical_l2,exact_l2\n";

		const double initial_l2 = std::sqrt(
			problem.solution().dot(problem.mass_matrix() * problem.solution()));

		problem.set_time_step_callback(
			[&](std::size_t, double time, const femg::Vector &solution) {
				const double numerical_l2 = std::sqrt(
					solution.dot(problem.mass_matrix() * solution));
				const double exact_l2 = initial_l2 *
					std::exp(-diffusion * lambda * time);
				decay_out << std::setprecision(16)
									<< time << ","
									<< numerical_l2 << ","
									<< exact_l2 << "\n";
			});

		problem.solve();

		const EigenmodeHeatExactSolution exact_solution(eigenmode_data);
		const double error_l2 = problem.compute_l2_error(exact_solution);
		std::cout << "Final L2 error against exp(-mu*lambda*t)*phi: "
							<< error_l2 << "\n";

		return 0;
	} catch (const std::exception &error) {
		std::cerr << "Error: " << error.what() << "\n";
		return 1;
	}
}
