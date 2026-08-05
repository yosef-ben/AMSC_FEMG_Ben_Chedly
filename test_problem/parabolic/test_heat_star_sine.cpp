#include "parabolic_problem.hpp"
#include "star_sine_functions.hpp"

#include <iostream>
#include <memory>
#include <string>

int main(int argc, char *argv[]) {
	const double diffusion = 1.0;
	const double T = 1.0;
	const double deltat = 0.005;
	const double theta = 1.0;

	std::string graph_file = "data/star_4.txt";
	char *local_argv[] = {argv[0], graph_file.data()};

	femg::parabolic_problem problem(diffusion, T, deltat, theta);
	problem.set_initial_condition(std::make_shared<StarSineInitialCondition>());
	problem.set_source_function(std::make_shared<StarSineForcingTerm>());
	problem.set_dirichlet_vertices(
		{0, 1, 2, 3, 4},
		[](std::size_t, double) {
			return 0.0;
		});
	problem.set_output_directory("output/visualization/star_sine");

	if (argc >= 2) {
		problem.init(argc, argv);
	} else {
		problem.init(2, local_argv);
	}

	problem.set_coefficients();
	problem.assemble_matrices();
	problem.print_matrix_summary(std::cout);
	problem.solve();

	const StarSineExactSolution exact_solution;
	const double error_l2 = problem.compute_l2_error(exact_solution);
	std::cout << "Final L2 error against u = sin(2*pi*x): "
						<< error_l2 << "\n";

	return 0;
}
