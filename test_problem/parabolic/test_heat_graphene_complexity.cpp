#include "parabolic_problem.hpp"

#include <chrono>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <string>
#include <vector>

namespace {

std::string write_graphene_file(int n_cells) {
	std::filesystem::create_directories("output/timing/graphs");

	const std::string filename =
		"output/timing/graphs/graphene_13_" + std::to_string(n_cells) + ".txt";

	std::ofstream out(filename);
	if (!out) {
		throw std::runtime_error("Unable to write graph file: " + filename);
	}

	out << "12 13\n";
	out << "0.0 1.7320508075688772\n";
	out << "1.0 1.7320508075688772\n";
	out << "1.5 0.8660254037844386\n";
	out << "1.0 0.0\n";
	out << "0.0 0.0\n";
	out << "-0.5 0.8660254037844386\n";
	out << "2.5 0.8660254037844386\n";
	out << "3.0 1.7320508075688772\n";
	out << "4.0 1.7320508075688772\n";
	out << "4.5 0.8660254037844386\n";
	out << "4.0 0.0\n";
	out << "3.0 0.0\n";

	const int edges[13][2] = {
		{0, 1}, {1, 2}, {2, 3}, {3, 4}, {4, 5}, {5, 0}, {2, 6},
		{6, 7}, {7, 8}, {8, 9}, {9, 10}, {10, 11}, {11, 6}
	};

	for (const auto &edge : edges) {
		out << edge[0] << " " << edge[1] << " 1.0 " << n_cells << "\n";
	}

	return filename;
}

struct TimingRow {
	int n_cells = 0;
	int n_dofs = 0;
	int n_time_steps = 0;
	int repetitions = 0;
	double init_seconds = 0.0;
	double assembly_seconds = 0.0;
	double factorization_seconds = 0.0;
	double time_stepping_seconds = 0.0;
	double total_solve_seconds = 0.0;
	double total_seconds = 0.0;
};

TimingRow run_case(int n_cells) {
	const std::string graph_file = write_graphene_file(n_cells);
	std::string executable_name = "test_heat_graphene_complexity";
	std::string local_graph_file = graph_file;
	char *argv[] = {executable_name.data(), local_graph_file.data()};

	const double diffusion = 1.0;
	const double final_time = 0.1;
	const double time_step = 0.01;
	const double theta = 0.5;

	const auto total_start = std::chrono::steady_clock::now();

	femg::parabolic_problem problem(diffusion, final_time, time_step, theta);
	problem.set_output_enabled(false);
	problem.set_verbose(false);

	const auto init_start = std::chrono::steady_clock::now();
	problem.init(2, argv);
	problem.set_coefficients();
	const auto init_end = std::chrono::steady_clock::now();

	const auto assembly_start = std::chrono::steady_clock::now();
	problem.assemble_matrices();
	const auto assembly_end = std::chrono::steady_clock::now();

	const auto solve_timing = problem.solve_with_timing();
	const auto total_end = std::chrono::steady_clock::now();

	TimingRow row;
	row.n_cells = n_cells;
	row.n_dofs = static_cast<int>(problem.solution().size());
	row.n_time_steps = static_cast<int>(solve_timing.n_time_steps);
	row.init_seconds =
		std::chrono::duration<double>(init_end - init_start).count();
	row.assembly_seconds =
		std::chrono::duration<double>(assembly_end - assembly_start).count();
	row.factorization_seconds = solve_timing.factorization_seconds;
	row.time_stepping_seconds = solve_timing.time_stepping_seconds;
	row.total_solve_seconds = solve_timing.total_solve_seconds;
	row.total_seconds =
		std::chrono::duration<double>(total_end - total_start).count();

	return row;
}

TimingRow average_case(int n_cells, int repetitions) {
	TimingRow average;
	average.n_cells = n_cells;
	average.repetitions = repetitions;

	for (int rep = 0; rep < repetitions; ++rep) {
		const TimingRow current = run_case(n_cells);
		average.n_dofs = current.n_dofs;
		average.n_time_steps = current.n_time_steps;
		average.init_seconds += current.init_seconds;
		average.assembly_seconds += current.assembly_seconds;
		average.factorization_seconds += current.factorization_seconds;
		average.time_stepping_seconds += current.time_stepping_seconds;
		average.total_solve_seconds += current.total_solve_seconds;
		average.total_seconds += current.total_seconds;
	}

	const double scale = 1.0 / static_cast<double>(repetitions);
	average.init_seconds *= scale;
	average.assembly_seconds *= scale;
	average.factorization_seconds *= scale;
	average.time_stepping_seconds *= scale;
	average.total_solve_seconds *= scale;
	average.total_seconds *= scale;

	return average;
}

} // namespace

int main() {
	try {
		std::filesystem::create_directories("output/timing");

		const std::vector<int> cell_counts = {100, 200, 400, 800, 1600};
		const int repetitions = 3;

		std::ofstream out("output/timing/graphene_parabolic_complexity.csv");
		out << "n_cells_per_edge,n_dofs,n_time_steps,repetitions,"
				<< "init_seconds,assembly_seconds,factorization_seconds,"
				<< "time_stepping_seconds,total_solve_seconds,total_seconds\n";

		for (const int n_cells : cell_counts) {
			const TimingRow row = average_case(n_cells, repetitions);

			out << row.n_cells << ","
					<< row.n_dofs << ","
					<< row.n_time_steps << ","
					<< row.repetitions << ","
					<< std::setprecision(16)
					<< row.init_seconds << ","
					<< row.assembly_seconds << ","
					<< row.factorization_seconds << ","
					<< row.time_stepping_seconds << ","
					<< row.total_solve_seconds << ","
					<< row.total_seconds << "\n";

			std::cout << "n_cells = " << row.n_cells
								<< ", dofs = " << row.n_dofs
								<< ", average total = "
								<< row.total_seconds << " s\n";
		}

		std::cout << "Wrote output/timing/graphene_parabolic_complexity.csv\n";
		return 0;
	} catch (const std::exception &error) {
		std::cerr << "Error: " << error.what() << "\n";
		return 1;
	}
}
