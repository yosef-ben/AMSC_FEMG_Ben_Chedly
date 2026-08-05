#include "parabolic_problem.hpp"
#include "star_localized_heat_functions.hpp"

#include <iostream>
#include <memory>
#include <string>

int main(int argc, char *argv[]) {
	const double diffusion = 5.0;
	const double T = 1.0;
	const double deltat = 0.05;
	const double theta = 0.5;

	std::string graph_file = "data/star_4.txt";
	char *local_argv[] = {argv[0], graph_file.data()};

	femg::parabolic_problem problem(diffusion, T, deltat, theta);
	problem.set_reaction_coefficient(star_localized_heat_reaction);
	problem.set_initial_condition(
		std::make_shared<StarLocalizedHeatInitialCondition>(diffusion));
	problem.set_output_directory("output/visualization/star_localized");

	if (argc >= 2) {
		problem.init(argc, argv);
	} else {
		problem.init(2, local_argv);
	}

	problem.set_coefficients();
	problem.assemble_matrices();
	problem.print_matrix_summary(std::cout);
	problem.solve();

	const StarLocalizedHeatExactSolution exact_solution(diffusion);
	const double error_l2 = problem.compute_l2_error(exact_solution);
	std::cout << "Final L2 error against eigenmode relaxation: "
						<< error_l2 << "\n";

	return 0;
}
