#include "parabolic_problem.hpp"

#include <cmath>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <string>
#include <vector>

namespace {

double run_problem(double time_step, double theta) {
	const double diffusion = 1.0;
	const double final_time = 0.1;

	std::string executable_name = "test_heat_time_convergence";
	std::string graph_file = "data/interval_320.txt";
	char *argv[] = {executable_name.data(), graph_file.data()};

	femg::parabolic_problem problem(
		diffusion,
		final_time,
		time_step,
		theta);

	problem.set_output_enabled(false);
	problem.set_verbose(false);
	problem.init(2, argv);
	problem.set_coefficients();
	problem.assemble_matrices();
	problem.solve();

	const femg::parabolic_problem::ExactSolution exact_solution;
	return problem.compute_l2_error(exact_solution);
}

} // namespace

int main() {
	const std::vector<double> time_steps = {0.05, 0.025, 0.0125, 0.00625};

	std::filesystem::create_directories("output/convergence");
	std::ofstream convergence_file("output/convergence/time_convergence.csv");
	convergence_file << "dt,L2_BE,rate_BE,L2_CN,rate_CN\n";

	double previous_error_be = 0.0;
	double previous_error_cn = 0.0;
	double previous_dt = 0.0;

	std::cout << std::setw(12) << "dt"
						<< std::setw(18) << "BE L2"
						<< std::setw(12) << "BE rate"
						<< std::setw(18) << "CN L2"
						<< std::setw(12) << "CN rate" << "\n";

	for (const double dt : time_steps) {
		const double error_be = run_problem(dt, 1.0);
		const double error_cn = run_problem(dt, 0.5);

		double rate_be = 0.0;
		double rate_cn = 0.0;
		if (previous_error_be > 0.0) {
			rate_be = std::log(previous_error_be / error_be)
				/ std::log(previous_dt / dt);
			rate_cn = std::log(previous_error_cn / error_cn)
				/ std::log(previous_dt / dt);
		}

		std::cout << std::setw(12) << dt
							<< std::setw(18) << error_be;
		if (previous_error_be > 0.0) {
			std::cout << std::setw(12) << rate_be;
		} else {
			std::cout << std::setw(12) << "-";
		}

		std::cout << std::setw(18) << error_cn;
		if (previous_error_cn > 0.0) {
			std::cout << std::setw(12) << rate_cn;
		} else {
			std::cout << std::setw(12) << "-";
		}
		std::cout << "\n";

		convergence_file << dt << "," << error_be << ",";
		if (previous_error_be > 0.0) {
			convergence_file << rate_be;
		}
		convergence_file << "," << error_cn << ",";
		if (previous_error_cn > 0.0) {
			convergence_file << rate_cn;
		}
		convergence_file << "\n";

		previous_error_be = error_be;
		previous_error_cn = error_cn;
		previous_dt = dt;
	}

	return 0;
}
