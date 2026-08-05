#include "parabolic_problem.hpp"

#include <cmath>
#include <iostream>
#include <memory>
#include <string>

class ConstantInitialCondition
	: public femg::parabolic_problem::FunctionU0 {
public:
	double value(std::size_t, double, double) const override {
		return 5.0;
	}
};

class ConstantExactSolution
	: public femg::parabolic_problem::ExactSolution {
public:
	double value(std::size_t, double, double) const override {
		return 5.0;
	}
};

int main(int argc, char *argv[]) {
	const double diffusion = 1.0;
	const double T = 1.0;
	const double deltat = 0.05;
	const double theta = 1.0;

	std::string graph_file = "data/star_4.txt";
	char *local_argv[] = {argv[0], graph_file.data()};

	femg::parabolic_problem problem(diffusion, T, deltat, theta);
	problem.set_initial_condition(std::make_shared<ConstantInitialCondition>());
	problem.set_output_directory("output/visualization/star_constant");

	if (argc >= 2) {
		problem.init(argc, argv);
	} else {
		problem.init(2, local_argv);
	}

	problem.set_coefficients();
	problem.assemble_matrices();
	problem.print_matrix_summary(std::cout);
	problem.solve();

	const ConstantExactSolution exact_solution;
	const double error_l2 = problem.compute_l2_error(exact_solution);
	std::cout << "Final L2 error against u = 5: " << error_l2 << "\n";

	return 0;
}
