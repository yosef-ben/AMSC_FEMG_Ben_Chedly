#include "spectral_problem.hpp"

#include <iostream>

int main(int argc, char *argv[]) {
	try {
		femg::spectral_problem problem(10);
		problem.init(argc, argv);
		problem.set_output_directory("output/spectral/star");
		problem.set_eigenfunction_scale(0.25);

		problem.set_coefficients();
		problem.assembly();
		problem.solve();
		problem.sol_export();

		problem.print_eigenvalues(std::cout);
		std::cout << "Results written to output/spectral/star\n";
		return 0;
	} catch (const std::exception &error) {
		std::cerr << "Error: " << error.what() << "\n";
		return 1;
	}
}
