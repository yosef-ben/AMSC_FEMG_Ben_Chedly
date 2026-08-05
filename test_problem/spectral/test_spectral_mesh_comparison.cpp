#include "spectral_problem.hpp"

#include <algorithm>
#include <cmath>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <stdexcept>
#include <string>
#include <vector>

namespace {

struct EdgeData {
	int source = 0;
	int target = 0;
	double length = 1.0;
};

std::vector<EdgeData> nonuniform_tree_edges() {
	return {
		{0, 1, 1.0},
		{1, 2, 1.0}, {1, 3, 1.0},
		{2, 4, 0.5}, {2, 5, 0.5}, {3, 6, 0.5}, {3, 7, 0.5},
		{4, 8, 0.25}, {4, 9, 0.25}, {5, 10, 0.25}, {5, 11, 0.25},
		{6, 12, 0.25}, {6, 13, 0.25}, {7, 14, 0.25}, {7, 15, 0.25}
	};
}

std::string write_tree_file(
	const std::string &name,
	const std::vector<int> &n_cells) {
	std::filesystem::create_directories("output/mesh_comparison/graphs");

	const auto edges = nonuniform_tree_edges();
	if (n_cells.size() != edges.size()) {
		throw std::runtime_error("Invalid number of mesh counts.");
	}

	const std::string filename =
		"output/mesh_comparison/graphs/" + name + ".txt";
	std::ofstream out(filename);
	if (!out) {
		throw std::runtime_error("Unable to write graph file: " + filename);
	}

	out << "16 15\n";
	out << "0.0 3.0\n";
	out << "0.0 2.0\n";
	out << "-1.0 1.0\n";
	out << "1.0 1.0\n";
	out << "-1.5 0.0\n";
	out << "-0.5 0.0\n";
	out << "0.5 0.0\n";
	out << "1.5 0.0\n";
	out << "-1.75 -1.0\n";
	out << "-1.25 -1.0\n";
	out << "-0.75 -1.0\n";
	out << "-0.25 -1.0\n";
	out << "0.25 -1.0\n";
	out << "0.75 -1.0\n";
	out << "1.25 -1.0\n";
	out << "1.75 -1.0\n";

	for (std::size_t i = 0; i < edges.size(); ++i) {
		out << edges[i].source << " "
				<< edges[i].target << " "
				<< edges[i].length << " "
				<< n_cells[i] << "\n";
	}

	return filename;
}

std::vector<int> n_type_mesh(int n_cells_per_edge) {
	return std::vector<int>(nonuniform_tree_edges().size(), n_cells_per_edge);
}

std::vector<int> h_type_mesh(double h) {
	std::vector<int> n_cells;
	for (const auto &edge : nonuniform_tree_edges()) {
		n_cells.push_back(std::max(1, static_cast<int>(
			std::llround(edge.length / h))));
	}
	return n_cells;
}

struct SpectrumData {
	std::string name;
	std::string graph_file;
	std::size_t dofs = 0;
	std::vector<double> eigenvalues;
};

SpectrumData compute_spectrum(
	const std::string &name,
	const std::string &graph_file,
	std::size_t n_modes) {
	std::string executable_name = "test_spectral_mesh_comparison";
	std::string local_graph_file = graph_file;
	char *argv[] = {executable_name.data(), local_graph_file.data()};

	femg::spectral_problem problem(n_modes);
	problem.init(2, argv);
	problem.set_coefficients();
	problem.assembly();
	problem.solve();

	SpectrumData data;
	data.name = name;
	data.graph_file = graph_file;
	data.dofs = problem.dofs();
	data.eigenvalues.reserve(n_modes);
	for (std::size_t i = 0; i < n_modes; ++i) {
		data.eigenvalues.push_back(problem.eigenvalue(i));
	}
	return data;
}

} // namespace

int main() {
	try {
		const std::size_t n_modes = 30;
		const int n_type_cells = 60;

		const double total_length = 7.0;
		const int target_total_cells =
			static_cast<int>(nonuniform_tree_edges().size()) * n_type_cells;
		const double comparable_h =
			total_length / static_cast<double>(target_total_cells);
		const double minimum_n_type_h = 0.25 / static_cast<double>(n_type_cells);

		const std::string n_type_file =
			write_tree_file("tree_nonuniform_n_type", n_type_mesh(n_type_cells));
		const std::string h_type_file =
			write_tree_file("tree_nonuniform_h_type", h_type_mesh(comparable_h));
		const std::string h_min_file =
			write_tree_file("tree_nonuniform_h_min", h_type_mesh(minimum_n_type_h));

		const SpectrumData n_type =
			compute_spectrum("N-type", n_type_file, n_modes);
		const SpectrumData h_type =
			compute_spectrum("h-type", h_type_file, n_modes);
		const SpectrumData h_min =
			compute_spectrum("h-min", h_min_file, n_modes);

		std::filesystem::create_directories("output/mesh_comparison");
		std::ofstream out("output/mesh_comparison/tree_mesh_spectra.csv");
		out << "index,n_type,h_type,h_min\n";
		for (std::size_t i = 0; i < n_modes; ++i) {
			out << i << ","
					<< std::setprecision(16)
					<< n_type.eigenvalues[i] << ","
					<< h_type.eigenvalues[i] << ","
					<< h_min.eigenvalues[i] << "\n";
		}

		std::ofstream summary("output/mesh_comparison/tree_mesh_summary.csv");
		summary << "mesh,graph_file,dofs\n";
		summary << "N-type," << n_type.graph_file << "," << n_type.dofs << "\n";
		summary << "h-type," << h_type.graph_file << "," << h_type.dofs << "\n";
		summary << "h-min," << h_min.graph_file << "," << h_min.dofs << "\n";

		std::cout << "N-type dofs: " << n_type.dofs << "\n";
		std::cout << "h-type dofs: " << h_type.dofs << "\n";
		std::cout << "h-min dofs:  " << h_min.dofs << "\n";
		std::cout << "Wrote output/mesh_comparison/tree_mesh_spectra.csv\n";
		std::cout << "Wrote output/mesh_comparison/tree_mesh_summary.csv\n";

		return 0;
	} catch (const std::exception &error) {
		std::cerr << "Error: " << error.what() << "\n";
		return 1;
	}
}
