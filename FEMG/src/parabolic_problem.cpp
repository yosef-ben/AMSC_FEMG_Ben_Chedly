#include "parabolic_problem.hpp"

#include <Eigen/Sparse>
#include <Eigen/SparseLU>

#include <algorithm>
#include <chrono>
#include <cmath>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <sstream>
#include <stdexcept>
#include <utility>
#include <vector>

namespace femg {

parabolic_problem::parabolic_problem(
	double diffusion,
	double final_time,
	double time_step,
	double theta)
	: diffusion_(diffusion),
		final_time_(final_time),
		time_step_(time_step),
		theta_(theta),
		forcing_term_(std::make_shared<ForcingTerm>()),
		initial_condition_term_(std::make_shared<FunctionU0>()) {
	source_function_ = [this](std::size_t edge_index, double s, double time) {
		return forcing_term_->value(edge_index, s, time);
	};

	initial_condition_ = [this](std::size_t edge_index, double s, double time) {
		return initial_condition_term_->value(edge_index, s, time);
	};
}

void parabolic_problem::set_source_function(EdgeFunction source_function) {
	source_function_ = std::move(source_function);
}

void parabolic_problem::set_initial_condition(EdgeFunction initial_condition) {
	initial_condition_ = std::move(initial_condition);
}

void parabolic_problem::set_source_function(
	std::shared_ptr<const ForcingTerm> source_function) {
	forcing_term_ = std::move(source_function);
	source_function_ = [this](std::size_t edge_index, double s, double time) {
		return forcing_term_->value(edge_index, s, time);
	};
}

void parabolic_problem::set_initial_condition(
	std::shared_ptr<const FunctionU0> initial_condition) {
	initial_condition_term_ = std::move(initial_condition);
	initial_condition_ = [this](std::size_t edge_index, double s, double time) {
		return initial_condition_term_->value(edge_index, s, time);
	};
}

void parabolic_problem::set_reaction_coefficient(double reaction) {
	reaction_ = reaction;
}

void parabolic_problem::set_output_enabled(bool output_enabled) {
	output_enabled_ = output_enabled;
}

void parabolic_problem::set_output_directory(std::string output_directory) {
	output_directory_ = std::move(output_directory);
}

void parabolic_problem::set_verbose(bool verbose) {
	verbose_ = verbose;
}

void parabolic_problem::set_time_step_callback(TimeStepCallback callback) {
	time_step_callback_ = std::move(callback);
}

void parabolic_problem::set_dirichlet_vertices(
	std::vector<std::size_t> vertices,
	VertexFunction boundary_value) {
	dirichlet_vertices_ = std::move(vertices);
	dirichlet_value_ = std::move(boundary_value);

	if (lhs_matrix_.rows() > 0) {
		impose_dirichlet_on_lhs();
	}
}

void parabolic_problem::clear_dirichlet_vertices() {
	dirichlet_vertices_.clear();
	dirichlet_value_ = nullptr;
}

void parabolic_problem::set_coefficients() {
	U = Vector::Zero(static_cast<Eigen::Index>(n_dofs()));

	Vector counts = Vector::Zero(static_cast<Eigen::Index>(n_dofs()));

	auto [ei, ei_end] = boost::edges(graph_);
	for (auto edge_it = ei; edge_it != ei_end; ++edge_it) {
		const Edge edge = *edge_it;
		const std::size_t edge_index = graph_[edge].index;
		const int n_cells = graph_[edge].n_cells;
		const double h = graph_[edge].length / static_cast<double>(n_cells);

		for (int local_node = 0; local_node <= n_cells; ++local_node) {
			const std::size_t dof = edge_local_node_to_dof(edge, local_node);
			const double s = static_cast<double>(local_node) * h;

			U(static_cast<Eigen::Index>(dof)) +=
				initial_condition_(edge_index, s, 0.0);
			counts(static_cast<Eigen::Index>(dof)) += 1.0;
		}
	}

	for (Eigen::Index i = 0; i < U.size(); ++i) {
		if (counts(i) > 0.0) {
			U(i) /= counts(i);
		}
	}

	if (!dirichlet_vertices_.empty()) {
		for (const auto vertex_index : dirichlet_vertices_) {
			const auto dof = static_cast<Eigen::Index>(
				vertex_dof_offset_ + vertex_index);
			U(dof) = dirichlet_value_(vertex_index, 0.0);
		}
	}
}

void parabolic_problem::assembly() {
	assemble_matrices();
}

void parabolic_problem::assemble_matrices() {
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

		const double k00 = diffusion_ / h;
		const double k01 = -diffusion_ / h;
		const double k10 = -diffusion_ / h;
		const double k11 = diffusion_ / h;

		const double r00 = reaction_ * h / 3.0;
		const double r01 = reaction_ * h / 6.0;
		const double r10 = reaction_ * h / 6.0;
		const double r11 = reaction_ * h / 3.0;

		for (int cell = 0; cell < n_cells; ++cell) {
			const std::size_t dof0 = edge_local_node_to_dof(edge, cell);
			const std::size_t dof1 = edge_local_node_to_dof(edge, cell + 1);

			const auto i0 = static_cast<Eigen::Index>(dof0);
			const auto i1 = static_cast<Eigen::Index>(dof1);

			mass_triplets.emplace_back(i0, i0, m00);
			mass_triplets.emplace_back(i0, i1, m01);
			mass_triplets.emplace_back(i1, i0, m10);
			mass_triplets.emplace_back(i1, i1, m11);

			stiffness_triplets.emplace_back(i0, i0, k00 + r00);
			stiffness_triplets.emplace_back(i0, i1, k01 + r01);
			stiffness_triplets.emplace_back(i1, i0, k10 + r10);
			stiffness_triplets.emplace_back(i1, i1, k11 + r11);
		}
	}

	const auto size = static_cast<Eigen::Index>(n_dofs());
	M.resize(size, size);
	H.resize(size, size);

	M.setFromTriplets(mass_triplets.begin(), mass_triplets.end());
	H.setFromTriplets(stiffness_triplets.begin(), stiffness_triplets.end());

	build_theta_matrices();
	impose_dirichlet_on_lhs();
}

bool parabolic_problem::solve() {
	solve_with_timing();
	return true;
}

parabolic_problem::TimingData parabolic_problem::solve_with_timing() {
	if (lhs_matrix_.rows() == 0 || rhs_matrix_.rows() == 0) {
		throw std::runtime_error("Matrices have not been assembled yet.");
	}

	TimingData timing;
	const auto solve_start = std::chrono::steady_clock::now();
	const auto factorization_start = std::chrono::steady_clock::now();
	linear_solver_.compute(lhs_matrix_);
	const auto factorization_end = std::chrono::steady_clock::now();
	timing.factorization_seconds =
		std::chrono::duration<double>(
			factorization_end - factorization_start).count();

	if (linear_solver_.info() != Eigen::Success) {
		throw std::runtime_error("Unable to factorize the theta-method matrix.");
	}

	time_ = 0.0;
	std::size_t time_step_number = 0;

	if (verbose_) {
		std::cout << "Solving parabolic problem\n";
		std::cout << "  n = " << time_step_number
							<< ", t = " << time_
							<< ", ||u|| = " << U.norm() << "\n";
	}
	if (output_enabled_) {
		const auto output_start = std::chrono::steady_clock::now();
		output(time_step_number);
		const auto output_end = std::chrono::steady_clock::now();
		timing.output_seconds +=
			std::chrono::duration<double>(output_end - output_start).count();
	}
	if (time_step_callback_) {
		time_step_callback_(time_step_number, time_, U);
	}

	const auto time_stepping_start = std::chrono::steady_clock::now();
	while (time_ < final_time_ - 0.5 * time_step_) {
		time_ += time_step_;
		++time_step_number;

		assemble_rhs(time_);
		solve_time_step();

		if (verbose_) {
			std::cout << "  n = " << time_step_number
								<< ", t = " << time_
								<< ", ||u|| = " << U.norm() << "\n";
		}
		if (output_enabled_) {
			const auto output_start = std::chrono::steady_clock::now();
			output(time_step_number);
			const auto output_end = std::chrono::steady_clock::now();
			timing.output_seconds +=
				std::chrono::duration<double>(output_end - output_start).count();
		}
		if (time_step_callback_) {
			time_step_callback_(time_step_number, time_, U);
		}
	}
	const auto time_stepping_end = std::chrono::steady_clock::now();
	timing.time_stepping_seconds =
		std::chrono::duration<double>(
			time_stepping_end - time_stepping_start).count()
		- timing.output_seconds;

	if (output_enabled_) {
		const auto output_start = std::chrono::steady_clock::now();
		write_pvd_record(time_step_number);
		const auto output_end = std::chrono::steady_clock::now();
		timing.output_seconds +=
			std::chrono::duration<double>(output_end - output_start).count();
	}

	const auto solve_end = std::chrono::steady_clock::now();
	timing.total_solve_seconds =
		std::chrono::duration<double>(solve_end - solve_start).count();
	timing.n_time_steps = time_step_number;

	return timing;
}

void parabolic_problem::sol_export() {}

void parabolic_problem::build_theta_matrices() {
	lhs_matrix_ = M + theta_ * time_step_ * H;
	rhs_matrix_ = M - (1.0 - theta_) * time_step_ * H;
}

Vector parabolic_problem::assemble_load_vector(double time) const {
	Vector load = Vector::Zero(static_cast<Eigen::Index>(n_dofs()));

	const double reference_points[2] = {
		0.5 * (1.0 - 1.0 / std::sqrt(3.0)),
		0.5 * (1.0 + 1.0 / std::sqrt(3.0))
	};
	const double reference_weights[2] = {0.5, 0.5};

	auto [ei, ei_end] = boost::edges(graph_);
	for (auto edge_it = ei; edge_it != ei_end; ++edge_it) {
		const Edge edge = *edge_it;
		const std::size_t edge_index = graph_[edge].index;
		const int n_cells = graph_[edge].n_cells;
		const double h = graph_[edge].length / static_cast<double>(n_cells);

		for (int cell = 0; cell < n_cells; ++cell) {
			const std::size_t dof0 = edge_local_node_to_dof(edge, cell);
			const std::size_t dof1 = edge_local_node_to_dof(edge, cell + 1);

			double local_load_0 = 0.0;
			double local_load_1 = 0.0;

			for (int q = 0; q < 2; ++q) {
				const double xi = reference_points[q];
				const double w = reference_weights[q];
				const double s = (static_cast<double>(cell) + xi) * h;
				const double f = source_function_(edge_index, s, time);

				const double phi0 = 1.0 - xi;
				const double phi1 = xi;

				local_load_0 += f * phi0 * w * h;
				local_load_1 += f * phi1 * w * h;
			}

			load(static_cast<Eigen::Index>(dof0)) += local_load_0;
			load(static_cast<Eigen::Index>(dof1)) += local_load_1;
		}
	}

	return load;
}

double parabolic_problem::compute_l2_error(
	const ExactSolution &exact_solution) const {
	double error_squared = 0.0;

	const double reference_points[2] = {
		0.5 * (1.0 - 1.0 / std::sqrt(3.0)),
		0.5 * (1.0 + 1.0 / std::sqrt(3.0))
	};
	const double reference_weights[2] = {0.5, 0.5};

	auto [ei, ei_end] = boost::edges(graph_);
	for (auto edge_it = ei; edge_it != ei_end; ++edge_it) {
		const Edge edge = *edge_it;
		const std::size_t edge_index = graph_[edge].index;
		const int n_cells = graph_[edge].n_cells;
		const double h = graph_[edge].length / static_cast<double>(n_cells);

		for (int cell = 0; cell < n_cells; ++cell) {
			const std::size_t dof0 = edge_local_node_to_dof(edge, cell);
			const std::size_t dof1 = edge_local_node_to_dof(edge, cell + 1);

			const double u0 = U(static_cast<Eigen::Index>(dof0));
			const double u1 = U(static_cast<Eigen::Index>(dof1));

			for (int q = 0; q < 2; ++q) {
				const double xi = reference_points[q];
				const double w = reference_weights[q];
				const double s = (static_cast<double>(cell) + xi) * h;

				const double uh = (1.0 - xi) * u0 + xi * u1;
				const double ue = exact_solution.value(edge_index, s, time_);
				const double error = uh - ue;

				error_squared += error * error * w * h;
			}
		}
	}

	return std::sqrt(error_squared);
}

void parabolic_problem::assemble_rhs(double new_time) {
	const double old_time = new_time - time_step_;
	const Vector load_new = assemble_load_vector(new_time);
	const Vector load_old = assemble_load_vector(old_time);

	system_rhs_ = rhs_matrix_ * U
		+ time_step_ * (theta_ * load_new + (1.0 - theta_) * load_old);

	impose_dirichlet_on_rhs(new_time);
}

void parabolic_problem::solve_time_step() {
	const Vector new_solution = linear_solver_.solve(system_rhs_);

	if (linear_solver_.info() != Eigen::Success) {
		throw std::runtime_error("Linear solve failed during time stepping.");
	}

	U = new_solution;
}

bool parabolic_problem::is_dirichlet_dof(Eigen::Index dof) const {
	for (const auto vertex_index : dirichlet_vertices_) {
		const auto vertex_dof = static_cast<Eigen::Index>(
			vertex_dof_offset_ + vertex_index);
		if (dof == vertex_dof) {
			return true;
		}
	}

	return false;
}

void parabolic_problem::impose_dirichlet_on_lhs() {
	if (dirichlet_vertices_.empty()) {
		return;
	}
	if (!dirichlet_value_) {
		throw std::runtime_error("Missing Dirichlet boundary value function.");
	}

	std::vector<Eigen::Triplet<double>> triplets;
	triplets.reserve(static_cast<std::size_t>(lhs_matrix_.nonZeros())
		+ dirichlet_vertices_.size());

	for (int k = 0; k < lhs_matrix_.outerSize(); ++k) {
		for (SparseMatrix::InnerIterator it(lhs_matrix_, k); it; ++it) {
			if (!is_dirichlet_dof(it.row())) {
				triplets.emplace_back(it.row(), it.col(), it.value());
			}
		}
	}

	for (const auto vertex_index : dirichlet_vertices_) {
		const auto dof = static_cast<Eigen::Index>(
			vertex_dof_offset_ + vertex_index);
		triplets.emplace_back(dof, dof, 1.0);
	}

	const auto rows = lhs_matrix_.rows();
	const auto cols = lhs_matrix_.cols();
	lhs_matrix_.resize(rows, cols);
	lhs_matrix_.setFromTriplets(triplets.begin(), triplets.end());
	lhs_matrix_.makeCompressed();
}

void parabolic_problem::impose_dirichlet_on_rhs(double time) {
	if (dirichlet_vertices_.empty()) {
		return;
	}
	if (!dirichlet_value_) {
		throw std::runtime_error("Missing Dirichlet boundary value function.");
	}

	for (const auto vertex_index : dirichlet_vertices_) {
		const auto dof = static_cast<Eigen::Index>(
			vertex_dof_offset_ + vertex_index);
		system_rhs_(dof) = dirichlet_value_(vertex_index, time);
	}
}

void parabolic_problem::output(std::size_t time_step_number) const {
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

	std::ostringstream filename;
	filename << output_directory_ << "/solution_"
					 << std::setw(4) << std::setfill('0')
					 << time_step_number << ".vtp";

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
	out << "      <PointData Scalars=\"u\">\n";
	out << "        <DataArray type=\"Float64\" Name=\"u\" "
			<< "NumberOfComponents=\"1\" format=\"ascii\">\n";

	for (const auto &edge : edges) {
		const int n_cells = graph_[edge].n_cells;

		for (int local_node = 0; local_node <= n_cells; ++local_node) {
			const std::size_t dof = edge_local_node_to_dof(edge, local_node);
			out << "          " << U(static_cast<Eigen::Index>(dof)) << "\n";
		}
	}

	out << "        </DataArray>\n";
	out << "      </PointData>\n";
	out << "    </Piece>\n";
	out << "  </PolyData>\n";
	out << "</VTKFile>\n";
}

void parabolic_problem::write_pvd_record(
	std::size_t last_time_step_number) const {
	std::filesystem::create_directories(output_directory_);

	const std::string pvd_path = output_directory_ + "/solution.pvd";
	std::ofstream out(pvd_path);
	if (!out) {
		throw std::runtime_error("Unable to open output file: " + pvd_path);
	}

	out << "<?xml version=\"1.0\"?>\n";
	out << "<VTKFile type=\"Collection\" version=\"0.1\" "
			<< "byte_order=\"LittleEndian\">\n";
	out << "  <Collection>\n";

	for (std::size_t step = 0; step <= last_time_step_number; ++step) {
		std::ostringstream filename;
		filename << "solution_" << std::setw(4) << std::setfill('0')
						 << step << ".vtp";

		out << "    <DataSet timestep=\""
				<< static_cast<double>(step) * time_step_
				<< "\" group=\"\" part=\"0\" file=\""
				<< filename.str() << "\"/>\n";
	}

	out << "  </Collection>\n";
	out << "</VTKFile>\n";
}

void parabolic_problem::print_matrix_summary(std::ostream &out) const {
	out << "Parabolic problem summary\n";
	out << "  diffusion:  " << diffusion_ << "\n";
	out << "  final time: " << final_time_ << "\n";
	out << "  time step:  " << time_step_ << "\n";
	out << "  theta:      " << theta_ << "\n";
	out << "  reaction:   " << reaction_ << "\n";
	out << "  dofs:       " << n_dofs() << "\n";
	out << "  M:          " << M.rows() << " x " << M.cols()
			<< ", nnz = " << M.nonZeros() << "\n";
	out << "  H:          " << H.rows() << " x " << H.cols()
			<< ", nnz = " << H.nonZeros() << "\n";
	out << "  lhs:        " << lhs_matrix_.rows() << " x " << lhs_matrix_.cols()
			<< ", nnz = " << lhs_matrix_.nonZeros() << "\n";
	out << "  rhs:        " << rhs_matrix_.rows() << " x " << rhs_matrix_.cols()
			<< ", nnz = " << rhs_matrix_.nonZeros() << "\n";
}

} // namespace femg
