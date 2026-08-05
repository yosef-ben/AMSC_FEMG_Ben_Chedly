#include "parabolic_problem.hpp"

#include <cmath>
#include <iostream>
#include <memory>
#include <string>
#include <vector>

double linear_star_value(std::size_t edge_index, double s) {
	if (edge_index == 0 || edge_index == 1) {
		return s;
	}
	return -s;
}

double linear_star_vertex_value(std::size_t vertex_index) {
	if (vertex_index == 1 || vertex_index == 2) {
		return 1.0;
	}
	return -1.0;
}

class LinearInitialCondition
	: public femg::parabolic_problem::FunctionU0 {
public:
	double value(std::size_t edge_index, double s, double) const override {
		return linear_star_value(edge_index, s);
	}
};

class LinearExactSolution
	: public femg::parabolic_problem::ExactSolution {
public:
	double value(std::size_t edge_index, double s, double) const override {
		return linear_star_value(edge_index, s);
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
	problem.set_initial_condition(std::make_shared<LinearInitialCondition>());
	problem.set_dirichlet_vertices(
		{1, 2, 3, 4},
		[](std::size_t vertex_index, double) {
			return linear_star_vertex_value(vertex_index);
		});
	problem.set_output_directory("output/visualization/star_linear");

	if (argc >= 2) {
		problem.init(argc, argv);
	} else {
		problem.init(2, local_argv);
	}

	problem.set_coefficients();
	problem.assemble_matrices();
	problem.print_matrix_summary(std::cout);
	problem.solve();

	const LinearExactSolution exact_solution;
	const double error_l2 = problem.compute_l2_error(exact_solution);
	std::cout << "Final L2 error against u = x + y: " << error_l2 << "\n";

	return 0;
}
