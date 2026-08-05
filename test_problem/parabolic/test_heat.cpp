#include "parabolic_problem.hpp"

#include <iostream>
#include <string>

int main(int argc, char *argv[]) {
	const double diffusion = 1.0;
	const double T = 1.0;
	const double deltat = 0.05;
	const double theta = 1.0;

	std::string graph_file = "data/interval_1d.txt";
	char *local_argv[] = {argv[0], graph_file.data()};

	femg::parabolic_problem problem(diffusion, T, deltat, theta);

	if (argc >= 2) {
		problem.init(argc, argv);
	} else {
		problem.init(2, local_argv);
	}

	problem.set_coefficients();
	problem.assemble_matrices();
	problem.print_matrix_summary(std::cout);
	problem.solve();

	return 0;
}
