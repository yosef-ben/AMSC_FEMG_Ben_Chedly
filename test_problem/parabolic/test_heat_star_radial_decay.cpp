#include "parabolic_problem.hpp"
#include "star_radial_decay_functions.hpp"

#include <iostream>
#include <memory>
#include <string>

int main(int argc, char *argv[]) {
	const double diffusion = 1.0;
	const double T = 1.0;
	const double deltat = 0.05;
	const double theta = 1.0;

	std::string graph_file = "data/star_4.txt";
	char *local_argv[] = {argv[0], graph_file.data()};

	femg::parabolic_problem problem(diffusion, T, deltat, theta);
	problem.set_reaction_coefficient(star_radial_decay_reaction);
	problem.set_initial_condition(
		std::make_shared<StarRadialDecayInitialCondition>());
	problem.set_dirichlet_vertices(
		{1, 2, 3, 4},
		[](std::size_t, double) {
			return 0.0;
		});
	problem.set_output_directory("output/visualization/star_radial_decay");

	if (argc >= 2) {
		problem.init(argc, argv);
	} else {
		problem.init(2, local_argv);
	}

	problem.set_coefficients();
	problem.assemble_matrices();
	problem.print_matrix_summary(std::cout);
	problem.solve();

	const StarRadialDecayExactSolution exact_solution;
	const double error_l2 = problem.compute_l2_error(exact_solution);
	std::cout << "Final L2 error against radial decay solution: "
						<< error_l2 << "\n";

	return 0;
}
