#include "spectral_problem.hpp"

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
		"output/timing/graphs/graphene_13_spectral_"
		+ std::to_string(n_cells) + ".txt";

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
	double init_seconds = 0.0;
	double assembly_seconds = 0.0;
	double solve_seconds = 0.0;
	double total_seconds = 0.0;
};

TimingRow run_case(int n_cells) {
	const std::string graph_file = write_graphene_file(n_cells);
	std::string executable_name = "test_spectral_graphene_complexity";
	std::string local_graph_file = graph_file;
	char *argv[] = {executable_name.data(), local_graph_file.data()};

	const auto total_start = std::chrono::steady_clock::now();

	femg::spectral_problem problem(12);

	const auto init_start = std::chrono::steady_clock::now();
	problem.init(2, argv);
	problem.set_coefficients();
	const auto init_end = std::chrono::steady_clock::now();

	const auto assembly_start = std::chrono::steady_clock::now();
	problem.assembly();
	const auto assembly_end = std::chrono::steady_clock::now();

	const auto solve_start = std::chrono::steady_clock::now();
	problem.solve();
	const auto solve_end = std::chrono::steady_clock::now();

	const auto total_end = std::chrono::steady_clock::now();

	TimingRow row;
	row.n_cells = n_cells;
	row.n_dofs = 13 * n_cells - 1;
	row.init_seconds =
		std::chrono::duration<double>(init_end - init_start).count();
	row.assembly_seconds =
		std::chrono::duration<double>(assembly_end - assembly_start).count();
	row.solve_seconds =
		std::chrono::duration<double>(solve_end - solve_start).count();
	row.total_seconds =
		std::chrono::duration<double>(total_end - total_start).count();

	return row;
}

} // namespace

int main() {
	try {
		std::filesystem::create_directories("output/timing");

		const std::vector<int> cell_counts = {10, 20, 40, 80, 120};

		std::ofstream out("output/timing/graphene_spectral_complexity.csv");
		out << "n_cells_per_edge,n_dofs,init_seconds,assembly_seconds,"
				<< "solve_seconds,total_seconds\n";

		for (const int n_cells : cell_counts) {
			const TimingRow row = run_case(n_cells);

			out << row.n_cells << ","
					<< row.n_dofs << ","
					<< std::setprecision(16)
					<< row.init_seconds << ","
					<< row.assembly_seconds << ","
					<< row.solve_seconds << ","
					<< row.total_seconds << "\n";

			std::cout << "n_cells = " << row.n_cells
							<< ", dofs = " << row.n_dofs
							<< ", spectral solve = "
							<< row.solve_seconds << " s\n";
		}

		std::cout << "Wrote output/timing/graphene_spectral_complexity.csv\n";
		return 0;
	} catch (const std::exception &error) {
		std::cerr << "Error: " << error.what() << "\n";
		return 1;
	}
}