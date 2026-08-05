#include "spectral_problem.hpp"

#include <Eigen/Eigenvalues>

#include <algorithm>
#include <cmath>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <sstream>
#include <stdexcept>
#include <vector>

namespace femg {

spectral_problem::spectral_problem(std::size_t n_modes)
	: n_modes_(n_modes) {}

void spectral_problem::set_output_directory(std::string output_directory) {
	output_directory_ = std::move(output_directory);
}

void spectral_problem::set_eigenfunction_scale(double scale) {
	eigenfunction_scale_ = scale;
}

void spectral_problem::set_coefficients() {
	U = Vector::Zero(static_cast<Eigen::Index>(n_dofs()));
}

void spectral_problem::assembly() {
	assemble_matrices();
}

void spectral_problem::assemble_matrices() {
	std::vector<Eigen::Triplet<double>> mass_triplets;
	std::vector<Eigen::Triplet<double>> stiffness_triplets;

	mass_triplets.reserve(4 * n_dofs());
	stiffness_triplets.reserve(4 * n_dofs());

	auto [ei, ei_end] = boost::edges(graph_);
	for (auto edge_it = ei; edge_it != ei_end; ++edge_it) {
		const Edge edge = *edge_it;
		const int n_cells = graph_[edge].n_cells;
		const double h = graph_[edge].length / static_cast<double>(n_cells);

		const double m00 = h / 3.0;
		const double m01 = h / 6.0;
		const double m10 = h / 6.0;
		const double m11 = h / 3.0;

		const double k00 = 1.0 / h;
		const double k01 = -1.0 / h;
		const double k10 = -1.0 / h;
		const double k11 = 1.0 / h;

		for (int cell = 0; cell < n_cells; ++cell) {
			const std::size_t dof0 = edge_local_node_to_dof(edge, cell);
			const std::size_t dof1 = edge_local_node_to_dof(edge, cell + 1);

			const auto i0 = static_cast<Eigen::Index>(dof0);
			const auto i1 = static_cast<Eigen::Index>(dof1);

			mass_triplets.emplace_back(i0, i0, m00);
			mass_triplets.emplace_back(i0, i1, m01);
			mass_triplets.emplace_back(i1, i0, m10);
			mass_triplets.emplace_back(i1, i1, m11);

			stiffness_triplets.emplace_back(i0, i0, k00);
			stiffness_triplets.emplace_back(i0, i1, k01);
			stiffness_triplets.emplace_back(i1, i0, k10);
			stiffness_triplets.emplace_back(i1, i1, k11);
		}
	}

	const auto size = static_cast<Eigen::Index>(n_dofs());
	M.resize(size, size);
	H.resize(size, size);

	M.setFromTriplets(mass_triplets.begin(), mass_triplets.end());
	H.setFromTriplets(stiffness_triplets.begin(), stiffness_triplets.end());
}

bool spectral_problem::solve() {
	if (H.rows() == 0 || M.rows() == 0) {
		throw std::runtime_error("Matrices have not been assembled yet.");
	}

	const Eigen::MatrixXd dense_H = Eigen::MatrixXd(H);
	const Eigen::MatrixXd dense_M = Eigen::MatrixXd(M);

	Eigen::GeneralizedSelfAdjointEigenSolver<Eigen::MatrixXd> solver(
		dense_H, dense_M);

	if (solver.info() != Eigen::Success) {
		throw std::runtime_error("Generalized eigenvalue solve failed.");
	}

	const std::size_t available_modes =
		static_cast<std::size_t>(solver.eigenvalues().size());
	const std::size_t modes_to_store = std::min(n_modes_, available_modes);

	eigenvalues_ = solver.eigenvalues().head(
		static_cast<Eigen::Index>(modes_to_store));
	eigenvectors_.clear();
	eigenvectors_.reserve(modes_to_store);

	for (std::size_t mode = 0; mode < modes_to_store; ++mode) {
		eigenvectors_.push_back(solver.eigenvectors().col(
			static_cast<Eigen::Index>(mode)));
	}

	return true;
}

void spectral_problem::sol_export() {
	std::filesystem::create_directories(output_directory_);
	write_flat_domain();

	for (std::size_t mode = 0; mode < eigenvectors_.size(); ++mode) {
		write_eigenfunction(mode);
	}

	write_eigenvalues_csv(output_directory_ + "/eigenvalues.csv");
}

double spectral_problem::eigenvalue(std::size_t mode) const {
	if (mode >= static_cast<std::size_t>(eigenvalues_.size())) {
		throw std::out_of_range("Eigenvalue mode index out of range.");
	}

	return eigenvalues_(static_cast<Eigen::Index>(mode));
}

const Vector &spectral_problem::eigenvector(std::size_t mode) const {
	if (mode >= eigenvectors_.size()) {
		throw std::out_of_range("Eigenvector mode index out of range.");
	}

	return eigenvectors_[mode];
}

std::vector<std::vector<double>> spectral_problem::edge_nodal_values(
	std::size_t mode) const {
	const Vector &phi = eigenvector(mode);

	std::vector<std::vector<double>> values(n_edges_);

	auto [ei, ei_end] = boost::edges(graph_);
	for (auto edge_it = ei; edge_it != ei_end; ++edge_it) {
		const Edge edge = *edge_it;
		const std::size_t edge_index = graph_[edge].index;
		const int n_cells = graph_[edge].n_cells;

		values[edge_index].resize(static_cast<std::size_t>(n_cells + 1));
		for (int local_node = 0; local_node <= n_cells; ++local_node) {
			const std::size_t dof = edge_local_node_to_dof(edge, local_node);
			values[edge_index][static_cast<std::size_t>(local_node)] =
				phi(static_cast<Eigen::Index>(dof));
		}
	}

	return values;
}

std::vector<double> spectral_problem::edge_lengths() const {
	std::vector<double> lengths(n_edges_);

	auto [ei, ei_end] = boost::edges(graph_);
	for (auto edge_it = ei; edge_it != ei_end; ++edge_it) {
		const Edge edge = *edge_it;
		lengths[graph_[edge].index] = graph_[edge].length;
	}

	return lengths;
}

Vector spectral_problem::combinatorial_laplacian_eigenvalues() const {
	Eigen::MatrixXd laplacian =
		Eigen::MatrixXd::Zero(
			static_cast<Eigen::Index>(n_vertices_orig_),
			static_cast<Eigen::Index>(n_vertices_orig_));

	auto [ei, ei_end] = boost::edges(graph_);
	for (auto edge_it = ei; edge_it != ei_end; ++edge_it) {
		const Edge edge = *edge_it;
		const auto source = static_cast<Eigen::Index>(
			graph_[boost::source(edge, graph_)].index);
		const auto target = static_cast<Eigen::Index>(
			graph_[boost::target(edge, graph_)].index);

		laplacian(source, source) += 1.0;
		laplacian(target, target) += 1.0;
		laplacian(source, target) -= 1.0;
		laplacian(target, source) -= 1.0;
	}

	Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> solver(laplacian);
	if (solver.info() != Eigen::Success) {
		throw std::runtime_error("Combinatorial eigenvalue solve failed.");
	}

	return solver.eigenvalues();
}

void spectral_problem::write_eigenvalues_csv(
	const std::string &filename) const {
	std::ofstream out(filename);
	if (!out) {
		throw std::runtime_error("Unable to open output file: " + filename);
	}

	out << "mode,eigenvalue\n";
	for (std::size_t mode = 0; mode < eigenvectors_.size(); ++mode) {
		out << mode << "," << std::setprecision(16)
				<< eigenvalue(mode) << "\n";
	}
}

void spectral_problem::print_eigenvalues(std::ostream &out) const {
	out << "First " << eigenvectors_.size() << " eigenvalues\n";
	for (std::size_t mode = 0; mode < eigenvectors_.size(); ++mode) {
		out << "  lambda_" << mode << " = "
				<< std::setprecision(12) << eigenvalue(mode) << "\n";
	}
}

void spectral_problem::write_flat_domain() const {
	std::filesystem::create_directories(output_directory_);

	std::vector<Edge> edges;
	auto [ei, ei_end] = boost::edges(graph_);
	for (auto edge_it = ei; edge_it != ei_end; ++edge_it) {
		edges.push_back(*edge_it);
	}

	std::sort(edges.begin(), edges.end(),
						[this](const Edge &a, const Edge &b) {
							return graph_[a].index < graph_[b].index;
						});

	std::size_t n_points = 0;
	std::size_t n_line_cells = 0;
	for (const auto &edge : edges) {
		n_points += static_cast<std::size_t>(graph_[edge].n_cells + 1);
		n_line_cells += static_cast<std::size_t>(graph_[edge].n_cells);
	}

	const std::string filename = output_directory_ + "/domain.vtp";
	std::ofstream out(filename);
	if (!out) {
		throw std::runtime_error("Unable to open output file: " + filename);
	}

	out << "<?xml version=\"1.0\"?>\n";
	out << "<VTKFile type=\"PolyData\" version=\"0.1\" "
			<< "byte_order=\"LittleEndian\">\n";
	out << "  <PolyData>\n";
	out << "    <Piece NumberOfPoints=\"" << n_points
			<< "\" NumberOfLines=\"" << n_line_cells << "\">\n";
	out << "      <Points>\n";
	out << "        <DataArray type=\"Float64\" NumberOfComponents=\"3\" "
			<< "format=\"ascii\">\n";

	for (const auto &edge : edges) {
		const Vertex source = boost::source(edge, graph_);
		const Vertex target = boost::target(edge, graph_);
		const int n_cells = graph_[edge].n_cells;
		const auto source_coords = graph_[source].coords;
		const auto target_coords = graph_[target].coords;

		for (int local_node = 0; local_node <= n_cells; ++local_node) {
			const double t =
				static_cast<double>(local_node) / static_cast<double>(n_cells);
			const double x =
				(1.0 - t) * source_coords[0] + t * target_coords[0];
			const double y =
				(1.0 - t) * source_coords[1] + t * target_coords[1];
			out << "          " << x << " " << y << " 0\n";
		}
	}

	out << "        </DataArray>\n";
	out << "      </Points>\n";
	out << "      <Lines>\n";
	out << "        <DataArray type=\"Int64\" Name=\"connectivity\" "
			<< "format=\"ascii\">\n";

	std::size_t point_offset = 0;
	for (const auto &edge : edges) {
		const int n_cells = graph_[edge].n_cells;
		for (int cell = 0; cell < n_cells; ++cell) {
			out << "          " << point_offset + static_cast<std::size_t>(cell)
					<< " " << point_offset + static_cast<std::size_t>(cell + 1)
					<< "\n";
		}
		point_offset += static_cast<std::size_t>(n_cells + 1);
	}

	out << "        </DataArray>\n";
	out << "        <DataArray type=\"Int64\" Name=\"offsets\" "
			<< "format=\"ascii\">\n";
	for (std::size_t line = 0; line < n_line_cells; ++line) {
		out << "          " << 2 * (line + 1) << "\n";
	}
	out << "        </DataArray>\n";
	out << "      </Lines>\n";
	out << "    </Piece>\n";
	out << "  </PolyData>\n";
	out << "</VTKFile>\n";
}

void spectral_problem::write_eigenfunction(std::size_t mode) const {
	std::filesystem::create_directories(output_directory_);

	const Vector &phi = eigenvector(mode);

	std::vector<Edge> edges;
	auto [ei, ei_end] = boost::edges(graph_);
	for (auto edge_it = ei; edge_it != ei_end; ++edge_it) {
		edges.push_back(*edge_it);
	}

	std::sort(edges.begin(), edges.end(),
						[this](const Edge &a, const Edge &b) {
							return graph_[a].index < graph_[b].index;
						});

	std::size_t n_points = 0;
	std::size_t n_line_cells = 0;
	for (const auto &edge : edges) {
		n_points += static_cast<std::size_t>(graph_[edge].n_cells + 1);
		n_line_cells += static_cast<std::size_t>(graph_[edge].n_cells);
	}

	std::ostringstream filename;
	filename << output_directory_ << "/eigenmode_"
					 << std::setw(2) << std::setfill('0') << mode << ".vtp";

	std::ofstream out(filename.str());
	if (!out) {
		throw std::runtime_error("Unable to open output file: " + filename.str());
	}

	out << "<?xml version=\"1.0\"?>\n";
	out << "<VTKFile type=\"PolyData\" version=\"0.1\" "
			<< "byte_order=\"LittleEndian\">\n";
	out << "  <PolyData>\n";
	out << "    <Piece NumberOfPoints=\"" << n_points
			<< "\" NumberOfLines=\"" << n_line_cells << "\">\n";
	out << "      <Points>\n";
	out << "        <DataArray type=\"Float64\" NumberOfComponents=\"3\" "
			<< "format=\"ascii\">\n";

	for (const auto &edge : edges) {
		const Vertex source = boost::source(edge, graph_);
		const Vertex target = boost::target(edge, graph_);
		const int n_cells = graph_[edge].n_cells;
		const auto source_coords = graph_[source].coords;
		const auto target_coords = graph_[target].coords;

		for (int local_node = 0; local_node <= n_cells; ++local_node) {
			const double t =
				static_cast<double>(local_node) / static_cast<double>(n_cells);
			const double x =
				(1.0 - t) * source_coords[0] + t * target_coords[0];
			const double y =
				(1.0 - t) * source_coords[1] + t * target_coords[1];
			const std::size_t dof = edge_local_node_to_dof(edge, local_node);
			const double z = eigenfunction_scale_ *
				phi(static_cast<Eigen::Index>(dof));
			out << "          " << x << " " << y << " " << z << "\n";
		}
	}

	out << "        </DataArray>\n";
	out << "      </Points>\n";
	out << "      <Lines>\n";
	out << "        <DataArray type=\"Int64\" Name=\"connectivity\" "
			<< "format=\"ascii\">\n";

	std::size_t point_offset = 0;
	for (const auto &edge : edges) {
		const int n_cells = graph_[edge].n_cells;
		for (int cell = 0; cell < n_cells; ++cell) {
			out << "          " << point_offset + static_cast<std::size_t>(cell)
					<< " " << point_offset + static_cast<std::size_t>(cell + 1)
					<< "\n";
		}
		point_offset += static_cast<std::size_t>(n_cells + 1);
	}

	out << "        </DataArray>\n";
	out << "        <DataArray type=\"Int64\" Name=\"offsets\" "
			<< "format=\"ascii\">\n";
	for (std::size_t line = 0; line < n_line_cells; ++line) {
		out << "          " << 2 * (line + 1) << "\n";
	}
	out << "        </DataArray>\n";
	out << "      </Lines>\n";
	out << "      <PointData Scalars=\"phi\">\n";
	out << "        <DataArray type=\"Float64\" Name=\"phi\" "
			<< "NumberOfComponents=\"1\" format=\"ascii\">\n";

	for (const auto &edge : edges) {
		const int n_cells = graph_[edge].n_cells;
		for (int local_node = 0; local_node <= n_cells; ++local_node) {
			const std::size_t dof = edge_local_node_to_dof(edge, local_node);
			out << "          "
					<< phi(static_cast<Eigen::Index>(dof)) << "\n";
		}
	}

	out << "        </DataArray>\n";
	out << "        <DataArray type=\"Float64\" Name=\"lambda\" "
			<< "NumberOfComponents=\"1\" format=\"ascii\">\n";
	for (std::size_t point = 0; point < n_points; ++point) {
		out << "          " << eigenvalue(mode) << "\n";
	}
	out << "        </DataArray>\n";
	out << "      </PointData>\n";
	out << "    </Piece>\n";
	out << "  </PolyData>\n";
	out << "</VTKFile>\n";
}

} // namespace femg
