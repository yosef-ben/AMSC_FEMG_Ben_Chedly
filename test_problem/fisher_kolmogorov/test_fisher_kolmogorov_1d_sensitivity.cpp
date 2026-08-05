#include "fisher_kolmogorov_problem.hpp"

#include <cmath>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <stdexcept>
#include <string>
#include <vector>

int main(int argc, char *argv[]) {
	try {
		const std::string graph_file = argc >= 2
			? argv[1]
			: "data/interval_weickenmeier_200.txt";
		const std::string output_file = argc >= 3
			? argv[2]
			: "output/fisher_kolmogorov/one_dimensional_sensitivity/profiles.csv";
		const double final_time = 20.0;
		const double time_step = 0.1;
		const double baseline_diffusion = 1.0e-4;
		const double initial_peak = 0.1;
		const int n_cells = 200;
		const std::vector<double> multipliers = {1.0, 2.0, 4.0};

		const std::filesystem::path path(output_file);
		if (path.has_parent_path()) {
			std::filesystem::create_directories(path.parent_path());
		}
		std::ofstream output(output_file);
		if (!output) {
			throw std::runtime_error("Unable to write " + output_file);
		}
		output << "diffusion,alpha,time,x,c\n";
		output << std::setprecision(16);

		for (const double diffusion_multiplier : multipliers) {
			for (const double alpha : multipliers) {
				char *local_argv[] = {
					argv[0],
					const_cast<char *>(graph_file.c_str())
				};
				femg::fisher_kolmogorov_problem problem(
					final_time, time_step);
				problem.init(2, local_argv);
				problem.set_edge_diffusion_coefficients(
					std::vector<double>(
						problem.number_of_edges(),
						baseline_diffusion * diffusion_multiplier));
				problem.set_vertex_reaction_coefficients(
					std::vector<double>(
						problem.number_of_vertices(), alpha));
				problem.set_vertex_initial_condition(
					std::vector<double>(
						problem.number_of_vertices(), 0.0));
				problem.set_time_scheme(
					femg::fisher_kolmogorov_problem::TimeScheme::backward_euler);
				problem.set_newton_parameters(1.0e-11, 30);
				problem.set_output_enabled(false);
				problem.set_verbose(false);
				problem.set_coefficients();

				std::vector<double> initial_values(
					static_cast<std::size_t>(n_cells + 1), 0.0);
				initial_values[static_cast<std::size_t>(n_cells / 2)] =
					initial_peak;
				problem.set_edge_initial_values(0, initial_values);
				problem.assemble_matrices();

				problem.set_time_step_callback(
					[&](std::size_t step, double time, const femg::Vector &) {
						const auto values = problem.edge_values(0);
						for (int node = 0; node <= n_cells; ++node) {
							const double x = -1.0
								+ 2.0 * static_cast<double>(node)
									/ static_cast<double>(n_cells);
							output
								<< baseline_diffusion * diffusion_multiplier << ","
								<< alpha << "," << time << "," << x << ","
								<< values[static_cast<std::size_t>(node)] << "\n";
						}
					});
				problem.solve();

				const double minimum = problem.solution().minCoeff();
				const double maximum = problem.solution().maxCoeff();
				if (!std::isfinite(minimum) || !std::isfinite(maximum)
					|| minimum < -1.0e-8 || maximum > 1.0 + 1.0e-8) {
					throw std::runtime_error(
						"The numerical concentration left the physical range.");
				}
				std::cout << "d = "
					<< baseline_diffusion * diffusion_multiplier
					<< ", alpha = " << alpha
					<< ", final range = [" << minimum
					<< ", " << maximum << "]\n";
			}
		}

		std::cout << "Profiles written to " << output_file << "\n";
		return 0;
	} catch (const std::exception &error) {
		std::cerr << "Error: " << error.what() << "\n";
		return 1;
	}
}
