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
#include <vector>

namespace {

double run_case(
	const EigenmodeHeatData &eigenmode_data,
	const std::string &graph_file,
	double diffusion,
	double final_time,
	double deltat,
	double theta) {
	char executable_name[] = "test_heat_graphene_eigenmode_time_convergence";
	std::string local_graph_file = graph_file;
	char *argv[] = {executable_name, local_graph_file.data()};

	femg::parabolic_problem problem(diffusion, final_time, deltat, theta);
	problem.set_initial_condition(
		std::make_shared<EigenmodeHeatInitialCondition>(eigenmode_data));
	problem.set_output_enabled(false);
	problem.set_verbose(false);

	problem.init(2, argv);
	problem.set_coefficients();
	problem.assemble_matrices();
	problem.solve();

	const EigenmodeHeatExactSolution exact_solution(eigenmode_data);
	return problem.compute_l2_error(exact_solution);
}

void write_method_rows(
	std::ofstream &out,
	const std::string &method,
	double theta,
	const EigenmodeHeatData &eigenmode_data,
	const std::string &graph_file,
	double diffusion,
	double final_time,
	const std::vector<double> &time_steps) {
	double previous_dt = 0.0;
	double previous_error = 0.0;

	for (const double deltat : time_steps) {
		const double error = run_case(
			eigenmode_data, graph_file, diffusion, final_time, deltat, theta);

		out << method << "," << theta << "," << deltat << ","
				<< std::setprecision(16) << error << ",";

		if (previous_dt > 0.0 && previous_error > 0.0 && error > 0.0) {
			const double rate =
				std::log(previous_error / error) / std::log(previous_dt / deltat);
			out << rate;
		}

		out << "\n";

		previous_dt = deltat;
		previous_error = error;
	}
}

} // namespace

int main() {
	try {
		const std::string graph_file = "data/graphene_13.txt";
		const double diffusion = 1.0;
		const double final_time = 1.0;
		const std::size_t mode = 1;
		const std::vector<double> time_steps = {0.2, 0.1, 0.05, 0.025};

		char executable_name[] = "test_heat_graphene_eigenmode_time_convergence";
		std::string local_graph_file = graph_file;
		char *argv[] = {executable_name, local_graph_file.data()};

		femg::spectral_problem eigen_problem(10);
		eigen_problem.init(2, argv);
		eigen_problem.set_coefficients();
		eigen_problem.assembly();
		eigen_problem.solve();

		const double lambda = eigen_problem.eigenvalue(mode);
		const EigenmodeHeatData eigenmode_data(
			eigen_problem.edge_nodal_values(mode),
			eigen_problem.edge_lengths(),
			lambda,
			diffusion);

		std::filesystem::create_directories("output/convergence");
		std::ofstream out(
			"output/convergence/graphene_eigenmode_time_convergence.csv");
		out << "method,theta,dt,L2_error,rate\n";

		write_method_rows(
			out, "BE", 1.0, eigenmode_data, graph_file,
			diffusion, final_time, time_steps);
		write_method_rows(
			out, "CN", 0.5, eigenmode_data, graph_file,
			diffusion, final_time, time_steps);

		std::cout << "Using eigenmode " << mode
							<< " with lambda = " << lambda << "\n";
		std::cout << "Wrote output/convergence/"
							<< "graphene_eigenmode_time_convergence.csv\n";

		return 0;
	} catch (const std::exception &error) {
		std::cerr << "Error: " << error.what() << "\n";
		return 1;
	}
}
