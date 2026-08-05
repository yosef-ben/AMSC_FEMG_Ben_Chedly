#include "fisher_kolmogorov_problem.hpp"

#include <Eigen/Sparse>

#include <algorithm>
#include <cmath>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <sstream>
#include <stdexcept>
#include <utility>

namespace femg {

fisher_kolmogorov_problem::fisher_kolmogorov_problem(
	double final_time,
	double time_step)
	: final_time_(final_time), time_step_(time_step) {
	if (final_time_ <= 0.0 || time_step_ <= 0.0) {
		throw std::invalid_argument("Final time and time step must be positive.");
	}
}

void fisher_kolmogorov_problem::set_edge_diffusion_coefficients(
	std::vector<double> coefficients) {
	edge_diffusion_ = std::move(coefficients);
}

void fisher_kolmogorov_problem::set_vertex_reaction_coefficients(
	std::vector<double> coefficients) {
	vertex_alpha_ = std::move(coefficients);
}

void fisher_kolmogorov_problem::set_vertex_initial_condition(
	std::vector<double> values) {
	vertex_initial_condition_ = std::move(values);
}

void fisher_kolmogorov_problem::set_vertex_visualization_coordinates(
	std::vector<std::array<double, 3>> coordinates) {
	visualization_coordinates_ = std::move(coordinates);
}

void fisher_kolmogorov_problem::set_output_enabled(bool enabled) {
	output_enabled_ = enabled;
}

void fisher_kolmogorov_problem::set_output_directory(std::string directory) {
	output_directory_ = std::move(directory);
}

void fisher_kolmogorov_problem::set_output_stride(std::size_t stride) {
	if (stride == 0) {
		throw std::invalid_argument("Output stride must be positive.");
	}
	output_stride_ = stride;
}

void fisher_kolmogorov_problem::set_verbose(bool verbose) {
	verbose_ = verbose;
}

void fisher_kolmogorov_problem::set_time_scheme(TimeScheme scheme) {
	time_scheme_ = scheme;
}

void fisher_kolmogorov_problem::set_newton_parameters(
	double tolerance,
	std::size_t max_iterations) {
	if (tolerance <= 0.0 || max_iterations == 0) {
		throw std::invalid_argument(
			"Newton tolerance and iteration count must be positive.");
	}
	newton_tolerance_ = tolerance;
	newton_max_iterations_ = max_iterations;
}

void fisher_kolmogorov_problem::set_time_step_callback(
	TimeStepCallback callback) {
	time_step_callback_ = std::move(callback);
}

Edge fisher_kolmogorov_problem::edge_with_index(std::size_t edge_index) const {
	if (edge_index >= n_edges_) {
		throw std::out_of_range("Invalid edge index.");
	}
	auto [edge_it, edge_end] = boost::edges(graph_);
	for (; edge_it != edge_end; ++edge_it) {
		if (graph_[*edge_it].index == edge_index) {
			return *edge_it;
		}
	}
	throw std::runtime_error("Edge index was not found in the graph.");
}

void fisher_kolmogorov_problem::set_edge_initial_values(
	std::size_t edge_index,
	const std::vector<double> &values) {
	if (U.size() != static_cast<Eigen::Index>(n_dofs())) {
		throw std::runtime_error(
			"Call set_coefficients before setting edge initial values.");
	}
	const Edge edge = edge_with_index(edge_index);
	const int n_cells = graph_[edge].n_cells;
	if (values.size() != static_cast<std::size_t>(n_cells + 1)) {
		throw std::invalid_argument(
			"One initial value is required for each local edge node.");
	}
	if (std::any_of(values.begin(), values.end(), [](double value) {
			return value < 0.0 || value > 1.0;
	})) {
		throw std::invalid_argument("Initial concentrations must lie in [0, 1].");
	}
	for (int node = 0; node <= n_cells; ++node) {
		const auto dof = static_cast<Eigen::Index>(
			edge_local_node_to_dof(edge, node));
		U(dof) = values[static_cast<std::size_t>(node)];
	}
}

std::vector<double> fisher_kolmogorov_problem::edge_values(
	std::size_t edge_index) const {
	const Edge edge = edge_with_index(edge_index);
	const int n_cells = graph_[edge].n_cells;
	std::vector<double> values(static_cast<std::size_t>(n_cells + 1));
	for (int node = 0; node <= n_cells; ++node) {
		const auto dof = static_cast<Eigen::Index>(
			edge_local_node_to_dof(edge, node));
		values[static_cast<std::size_t>(node)] = U(dof);
	}
	return values;
}

void fisher_kolmogorov_problem::validate_input_data() const {
	if (!edge_diffusion_.empty() && edge_diffusion_.size() != n_edges_) {
		throw std::invalid_argument(
			"One diffusion coefficient is required for each edge.");
	}
	if (!vertex_alpha_.empty() && vertex_alpha_.size() != n_vertices_orig_) {
		throw std::invalid_argument(
			"One reaction coefficient is required for each graph vertex.");
	}
	if (!vertex_initial_condition_.empty()
		&& vertex_initial_condition_.size() != n_vertices_orig_) {
		throw std::invalid_argument(
			"One initial value is required for each graph vertex.");
	}
	if (!visualization_coordinates_.empty()
		&& visualization_coordinates_.size() != n_vertices_orig_) {
		throw std::invalid_argument(
			"One visualization coordinate is required for each graph vertex.");
	}
	if (std::any_of(
			edge_diffusion_.begin(), edge_diffusion_.end(),
			[](double value) { return value < 0.0; })) {
		throw std::invalid_argument("Diffusion coefficients must be nonnegative.");
	}
	if (std::any_of(
			vertex_alpha_.begin(), vertex_alpha_.end(),
			[](double value) { return value < 0.0; })) {
		throw std::invalid_argument("Reaction coefficients must be nonnegative.");
	}
	if (std::any_of(
			vertex_initial_condition_.begin(), vertex_initial_condition_.end(),
			[](double value) { return value < 0.0 || value > 1.0; })) {
		throw std::invalid_argument("Initial concentrations must lie in [0, 1].");
	}
}

void fisher_kolmogorov_problem::interpolate_vertex_data_to_dofs(
	const std::vector<double> &vertex_values,
	Vector &dof_values) const {
	dof_values = Vector::Zero(static_cast<Eigen::Index>(n_dofs()));

	for (std::size_t vertex = 0; vertex < n_vertices_orig_; ++vertex) {
		const auto dof = static_cast<Eigen::Index>(vertex_dof_offset_ + vertex);
		dof_values(dof) = vertex_values[vertex];
	}

	auto [edge_it, edge_end] = boost::edges(graph_);
	for (; edge_it != edge_end; ++edge_it) {
		const Edge edge = *edge_it;
		const Vertex source = boost::source(edge, graph_);
		const Vertex target = boost::target(edge, graph_);
		const double source_value = vertex_values[graph_[source].index];
		const double target_value = vertex_values[graph_[target].index];
		const int n_cells = graph_[edge].n_cells;

		for (int local_node = 1; local_node < n_cells; ++local_node) {
			const double xi =
				static_cast<double>(local_node) / static_cast<double>(n_cells);
			const std::size_t dof =
				edge_local_node_to_dof(edge, local_node);
			dof_values(static_cast<Eigen::Index>(dof)) =
				(1.0 - xi) * source_value + xi * target_value;
		}
	}
}

void fisher_kolmogorov_problem::set_coefficients() {
	if (edge_diffusion_.empty()) {
		edge_diffusion_.assign(n_edges_, 1.0);
	}
	if (vertex_alpha_.empty()) {
		vertex_alpha_.assign(n_vertices_orig_, 0.0);
	}
	if (vertex_initial_condition_.empty()) {
		vertex_initial_condition_.assign(n_vertices_orig_, 0.0);
	}
	validate_input_data();

	interpolate_vertex_data_to_dofs(vertex_initial_condition_, U);
	interpolate_vertex_data_to_dofs(vertex_alpha_, alpha_dofs_);
}

void fisher_kolmogorov_problem::assembly() {
	assemble_matrices();
}

void fisher_kolmogorov_problem::assemble_matrices() {
	validate_input_data();
	if (U.size() != static_cast<Eigen::Index>(n_dofs())) {
		throw std::runtime_error("Call set_coefficients before matrix assembly.");
	}

	std::vector<Eigen::Triplet<double>> mass_triplets;
	std::vector<Eigen::Triplet<double>> diffusion_triplets;
	mass_triplets.reserve(4 * n_dofs());
	diffusion_triplets.reserve(4 * n_dofs());

	auto [edge_it, edge_end] = boost::edges(graph_);
	for (; edge_it != edge_end; ++edge_it) {
		const Edge edge = *edge_it;
		const std::size_t edge_index = graph_[edge].index;
		const int n_cells = graph_[edge].n_cells;
		const double h = graph_[edge].length / static_cast<double>(n_cells);
		const double diffusion = edge_diffusion_.at(edge_index);

		for (int cell = 0; cell < n_cells; ++cell) {
			const auto i0 = static_cast<Eigen::Index>(
				edge_local_node_to_dof(edge, cell));
			const auto i1 = static_cast<Eigen::Index>(
				edge_local_node_to_dof(edge, cell + 1));

			mass_triplets.emplace_back(i0, i0, h / 3.0);
			mass_triplets.emplace_back(i0, i1, h / 6.0);
			mass_triplets.emplace_back(i1, i0, h / 6.0);
			mass_triplets.emplace_back(i1, i1, h / 3.0);

			diffusion_triplets.emplace_back(i0, i0, diffusion / h);
			diffusion_triplets.emplace_back(i0, i1, -diffusion / h);
			diffusion_triplets.emplace_back(i1, i0, -diffusion / h);
			diffusion_triplets.emplace_back(i1, i1, diffusion / h);
		}
	}

	const auto size = static_cast<Eigen::Index>(n_dofs());
	M.resize(size, size);
	H.resize(size, size);
	M.setFromTriplets(mass_triplets.begin(), mass_triplets.end());
	H.setFromTriplets(diffusion_triplets.begin(), diffusion_triplets.end());
	M.makeCompressed();
	H.makeCompressed();

	lhs_matrix_ = M + 0.5 * time_step_ * H;
	rhs_matrix_ = M - 0.5 * time_step_ * H;
	lhs_matrix_.makeCompressed();
	rhs_matrix_.makeCompressed();
}

Vector fisher_kolmogorov_problem::assemble_reaction_vector(
	const Vector &state) const {
	if (state.size() != static_cast<Eigen::Index>(n_dofs())) {
		throw std::invalid_argument("Reaction state has the wrong size.");
	}

	Vector reaction = Vector::Zero(static_cast<Eigen::Index>(n_dofs()));
	const double root = std::sqrt(3.0 / 5.0);
	const double points[3] = {
		0.5 * (1.0 - root),
		0.5,
		0.5 * (1.0 + root)
	};
	const double weights[3] = {5.0 / 18.0, 4.0 / 9.0, 5.0 / 18.0};

	auto [edge_it, edge_end] = boost::edges(graph_);
	for (; edge_it != edge_end; ++edge_it) {
		const Edge edge = *edge_it;
		const int n_cells = graph_[edge].n_cells;
		const double h = graph_[edge].length / static_cast<double>(n_cells);

		for (int cell = 0; cell < n_cells; ++cell) {
			const auto i0 = static_cast<Eigen::Index>(
				edge_local_node_to_dof(edge, cell));
			const auto i1 = static_cast<Eigen::Index>(
				edge_local_node_to_dof(edge, cell + 1));
			const double c0 = state(i0);
			const double c1 = state(i1);
			const double alpha0 = alpha_dofs_(i0);
			const double alpha1 = alpha_dofs_(i1);

			for (int q = 0; q < 3; ++q) {
				const double phi0 = 1.0 - points[q];
				const double phi1 = points[q];
				const double concentration = phi0 * c0 + phi1 * c1;
				const double alpha = phi0 * alpha0 + phi1 * alpha1;
				const double value =
					alpha * concentration * (1.0 - concentration);
				reaction(i0) += value * phi0 * weights[q] * h;
				reaction(i1) += value * phi1 * weights[q] * h;
			}
		}
	}
	return reaction;
}

SparseMatrix fisher_kolmogorov_problem::assemble_reaction_weight_matrix(
	const Vector &state,
	double state_factor) const {
	if (state.size() != static_cast<Eigen::Index>(n_dofs())) {
		throw std::invalid_argument("Reaction state has the wrong size.");
	}

	std::vector<Eigen::Triplet<double>> triplets;
	triplets.reserve(4 * n_dofs());
	const double root = std::sqrt(3.0 / 5.0);
	const double points[3] = {
		0.5 * (1.0 - root),
		0.5,
		0.5 * (1.0 + root)
	};
	const double weights[3] = {5.0 / 18.0, 4.0 / 9.0, 5.0 / 18.0};

	auto [edge_it, edge_end] = boost::edges(graph_);
	for (; edge_it != edge_end; ++edge_it) {
		const Edge edge = *edge_it;
		const int n_cells = graph_[edge].n_cells;
		const double h = graph_[edge].length / static_cast<double>(n_cells);

		for (int cell = 0; cell < n_cells; ++cell) {
			const auto i0 = static_cast<Eigen::Index>(
				edge_local_node_to_dof(edge, cell));
			const auto i1 = static_cast<Eigen::Index>(
				edge_local_node_to_dof(edge, cell + 1));
			const double c0 = state(i0);
			const double c1 = state(i1);
			const double alpha0 = alpha_dofs_(i0);
			const double alpha1 = alpha_dofs_(i1);

			double local[2][2] = {};
			for (int q = 0; q < 3; ++q) {
				const double basis[2] = {1.0 - points[q], points[q]};
				const double concentration = basis[0] * c0 + basis[1] * c1;
				const double alpha = basis[0] * alpha0 + basis[1] * alpha1;
				const double coefficient =
					alpha * (1.0 - state_factor * concentration);
				for (int i = 0; i < 2; ++i) {
					for (int j = 0; j < 2; ++j) {
						local[i][j] += coefficient * basis[i] * basis[j]
							* weights[q] * h;
					}
				}
			}
			triplets.emplace_back(i0, i0, local[0][0]);
			triplets.emplace_back(i0, i1, local[0][1]);
			triplets.emplace_back(i1, i0, local[1][0]);
			triplets.emplace_back(i1, i1, local[1][1]);
		}
	}

	SparseMatrix matrix(
		static_cast<Eigen::Index>(n_dofs()),
		static_cast<Eigen::Index>(n_dofs()));
	matrix.setFromTriplets(triplets.begin(), triplets.end());
	matrix.makeCompressed();
	return matrix;
}

Vector fisher_kolmogorov_problem::solve_backward_euler_step(
	const Vector &old_solution) {
	Vector iterate = old_solution;
	const double residual_scale =
		std::max(1.0, (M * old_solution).norm());
	double last_residual = 0.0;

	for (std::size_t iteration = 0;
		iteration < newton_max_iterations_; ++iteration) {
		const Vector residual = M * (iterate - old_solution)
			+ time_step_ * (H * iterate)
			- time_step_ * assemble_reaction_vector(iterate);
		last_residual = residual.norm();
		if (last_residual <= newton_tolerance_ * residual_scale) {
			return iterate;
		}

		SparseMatrix jacobian = M + time_step_ * H
			- time_step_
				* assemble_reaction_weight_matrix(iterate, 2.0);
		jacobian.makeCompressed();
		Eigen::SparseLU<SparseMatrix> newton_solver;
		newton_solver.compute(jacobian);
		if (newton_solver.info() != Eigen::Success) {
			throw std::runtime_error(
				"Unable to factorize the Newton Jacobian.");
		}
		const Vector update = newton_solver.solve(-residual);
		if (newton_solver.info() != Eigen::Success) {
			throw std::runtime_error("Newton linear solve failed.");
		}

		// Backtracking line search. The undamped step can leave the physical
		// range whenever the reaction dominates the diffusion, and the
		// Jacobian is only guaranteed positive definite while
		// dt * max(alpha) * |1 - 2c| < 1. Outside newton_admissible_range_
		// that guarantee is lost and the iteration stalls, so a candidate is
		// accepted only if it both reduces the residual and stays admissible.
		double step = 1.0;
		Vector candidate = iterate + update;
		for (std::size_t trial = 0; trial < newton_max_backtracks_; ++trial) {
			const bool admissible =
				candidate.minCoeff() >= -newton_admissible_range_
				&& candidate.maxCoeff() <= 1.0 + newton_admissible_range_;
			if (admissible) {
				const double candidate_residual =
					(M * (candidate - old_solution)
						+ time_step_ * (H * candidate)
						- time_step_ * assemble_reaction_vector(candidate)).norm();
				if (candidate_residual < last_residual) {
					break;
				}
			}
			step *= 0.5;
			candidate = iterate + step * update;
		}
		iterate = candidate;
	}
	std::ostringstream message;
	message << "Newton method did not converge within "
		<< newton_max_iterations_ << " iterations at t = " << time_
		<< ": residual " << last_residual << " > tolerance "
		<< newton_tolerance_ * residual_scale
		<< " (min c = " << iterate.minCoeff()
		<< ", max c = " << iterate.maxCoeff() << ").";
	throw std::runtime_error(message.str());
}

bool fisher_kolmogorov_problem::solve() {
	if (M.rows() == 0 || H.rows() == 0) {
		throw std::runtime_error("Matrices have not been assembled yet.");
	}

	const double step_count = final_time_ / time_step_;
	const auto n_steps = static_cast<std::size_t>(std::llround(step_count));
	if (std::abs(step_count - static_cast<double>(n_steps)) > 1.0e-10) {
		throw std::invalid_argument("Final time must be a multiple of the time step.");
	}

	time_ = 0.0;
	Vector previous_solution = U;
	std::vector<std::size_t> output_steps;
	if (output_enabled_) {
		output(0);
		output_steps.push_back(0);
	}
	if (time_step_callback_) {
		time_step_callback_(0, time_, U);
	}

	for (std::size_t step = 1; step <= n_steps; ++step) {
		const Vector old_solution = U;
		if (time_scheme_ == TimeScheme::backward_euler) {
			U = solve_backward_euler_step(old_solution);
		} else {
			// Corti et al.: only the factor (1-c) is extrapolated.
			const Vector extrapolated = (step == 1)
				? U
				: 1.5 * U - 0.5 * previous_solution;
			const SparseMatrix reaction_matrix =
				assemble_reaction_weight_matrix(extrapolated, 1.0);
			lhs_matrix_ = M + 0.5 * time_step_ * H
				- 0.5 * time_step_ * reaction_matrix;
			rhs_matrix_ = M - 0.5 * time_step_ * H
				+ 0.5 * time_step_ * reaction_matrix;
			lhs_matrix_.makeCompressed();
			rhs_matrix_.makeCompressed();

			linear_solver_.compute(lhs_matrix_);
			if (linear_solver_.info() != Eigen::Success) {
				throw std::runtime_error(
					"Unable to factorize the Corti time-step matrix.");
			}
			U = linear_solver_.solve(rhs_matrix_ * old_solution);
			if (linear_solver_.info() != Eigen::Success) {
				throw std::runtime_error(
					"Fisher-Kolmogorov time step failed.");
			}
		}

		previous_solution = old_solution;
		time_ = static_cast<double>(step) * time_step_;

		if (verbose_) {
			std::cout << "  n = " << step << ", t = " << time_
				<< ", min(c) = " << U.minCoeff()
				<< ", max(c) = " << U.maxCoeff() << "\n";
		}
		if (time_step_callback_) {
			time_step_callback_(step, time_, U);
		}
		if (output_enabled_
			&& (step % output_stride_ == 0 || step == n_steps)) {
			output(step);
			output_steps.push_back(step);
		}
	}

	if (output_enabled_) {
		write_pvd_record(output_steps);
	}
	return true;
}

void fisher_kolmogorov_problem::sol_export() {}

double fisher_kolmogorov_problem::vertex_value(
	std::size_t vertex_index) const {
	if (vertex_index >= n_vertices_orig_) {
		throw std::out_of_range("Invalid graph vertex index.");
	}
	return U(static_cast<Eigen::Index>(vertex_dof_offset_ + vertex_index));
}

void fisher_kolmogorov_problem::output(std::size_t step) const {
	std::filesystem::create_directories(output_directory_);
	std::vector<Edge> edges;
	auto [edge_it, edge_end] = boost::edges(graph_);
	for (; edge_it != edge_end; ++edge_it) {
		edges.push_back(*edge_it);
	}
	std::sort(edges.begin(), edges.end(), [this](const Edge &a, const Edge &b) {
		return graph_[a].index < graph_[b].index;
	});

	std::size_t n_points = 0;
	std::size_t n_lines = 0;
	for (const Edge &edge : edges) {
		n_points += static_cast<std::size_t>(graph_[edge].n_cells + 1);
		n_lines += static_cast<std::size_t>(graph_[edge].n_cells);
	}

	std::ostringstream path;
	path << output_directory_ << "/solution_" << std::setw(4)
		<< std::setfill('0') << step << ".vtp";
	std::ofstream out(path.str());
	if (!out) {
		throw std::runtime_error("Unable to write " + path.str());
	}

	out << "<?xml version=\"1.0\"?>\n"
		<< "<VTKFile type=\"PolyData\" version=\"0.1\" "
		<< "byte_order=\"LittleEndian\">\n"
		<< "  <PolyData>\n"
		<< "    <Piece NumberOfPoints=\"" << n_points
		<< "\" NumberOfLines=\"" << n_lines << "\">\n"
		<< "      <Points>\n"
		<< "        <DataArray type=\"Float64\" NumberOfComponents=\"3\" "
		<< "format=\"ascii\">\n";
	for (const Edge &edge : edges) {
		const Vertex source_vertex = boost::source(edge, graph_);
		const Vertex target_vertex = boost::target(edge, graph_);
		std::array<double, 3> source = {
			graph_[source_vertex].coords[0], graph_[source_vertex].coords[1], 0.0};
		std::array<double, 3> target = {
			graph_[target_vertex].coords[0], graph_[target_vertex].coords[1], 0.0};
		if (!visualization_coordinates_.empty()) {
			source = visualization_coordinates_.at(graph_[source_vertex].index);
			target = visualization_coordinates_.at(graph_[target_vertex].index);
		}
		const int n_cells = graph_[edge].n_cells;
		for (int node = 0; node <= n_cells; ++node) {
			const double xi = static_cast<double>(node) / n_cells;
			out << "          " << (1.0 - xi) * source[0] + xi * target[0]
				<< " " << (1.0 - xi) * source[1] + xi * target[1]
				<< " " << (1.0 - xi) * source[2] + xi * target[2] << "\n";
		}
	}
	out << "        </DataArray>\n      </Points>\n      <Lines>\n"
		<< "        <DataArray type=\"Int64\" Name=\"connectivity\" "
		<< "format=\"ascii\">\n";
	std::size_t point_offset = 0;
	for (const Edge &edge : edges) {
		const int n_cells = graph_[edge].n_cells;
		for (int cell = 0; cell < n_cells; ++cell) {
			out << "          " << point_offset + cell << " "
				<< point_offset + cell + 1 << "\n";
		}
		point_offset += static_cast<std::size_t>(n_cells + 1);
	}
	out << "        </DataArray>\n"
		<< "        <DataArray type=\"Int64\" Name=\"offsets\" "
		<< "format=\"ascii\">\n";
	for (std::size_t line = 0; line < n_lines; ++line) {
		out << "          " << 2 * (line + 1) << "\n";
	}
	out << "        </DataArray>\n      </Lines>\n"
		<< "      <PointData Scalars=\"c\">\n"
		<< "        <DataArray type=\"Float64\" Name=\"c\" "
		<< "NumberOfComponents=\"1\" format=\"ascii\">\n";
	for (const Edge &edge : edges) {
		for (int node = 0; node <= graph_[edge].n_cells; ++node) {
			const auto dof = static_cast<Eigen::Index>(
				edge_local_node_to_dof(edge, node));
			out << "          " << U(dof) << "\n";
		}
	}
	out << "        </DataArray>\n"
		<< "        <DataArray type=\"Float64\" Name=\"alpha\" "
		<< "NumberOfComponents=\"1\" format=\"ascii\">\n";
	for (const Edge &edge : edges) {
		for (int node = 0; node <= graph_[edge].n_cells; ++node) {
			const auto dof = static_cast<Eigen::Index>(
				edge_local_node_to_dof(edge, node));
			out << "          " << alpha_dofs_(dof) << "\n";
		}
	}
	out << "        </DataArray>\n      </PointData>\n"
		<< "      <CellData Scalars=\"diffusion\">\n"
		<< "        <DataArray type=\"Float64\" Name=\"diffusion\" "
		<< "NumberOfComponents=\"1\" format=\"ascii\">\n";
	for (const Edge &edge : edges) {
		const double diffusion = edge_diffusion_[graph_[edge].index];
		for (int cell = 0; cell < graph_[edge].n_cells; ++cell) {
			out << "          " << diffusion << "\n";
		}
	}
	out << "        </DataArray>\n      </CellData>\n"
		<< "    </Piece>\n  </PolyData>\n</VTKFile>\n";
}

void fisher_kolmogorov_problem::write_pvd_record(
	const std::vector<std::size_t> &steps) const {
	const std::string path = output_directory_ + "/solution.pvd";
	std::ofstream out(path);
	if (!out) {
		throw std::runtime_error("Unable to write " + path);
	}
	out << "<?xml version=\"1.0\"?>\n"
		<< "<VTKFile type=\"Collection\" version=\"0.1\" "
		<< "byte_order=\"LittleEndian\">\n  <Collection>\n";
	for (const std::size_t step : steps) {
		out << "    <DataSet timestep=\"" << step * time_step_
			<< "\" group=\"\" part=\"0\" file=\"solution_"
			<< std::setw(4) << std::setfill('0') << step << ".vtp\"/>\n";
	}
	out << "  </Collection>\n</VTKFile>\n";
}

void fisher_kolmogorov_problem::print_matrix_summary(std::ostream &out) const {
	out << "Fisher-Kolmogorov problem summary\n"
		<< "  final time: " << final_time_ << "\n"
		<< "  time step:  " << time_step_ << "\n"
		<< "  vertices:   " << n_vertices_orig_ << "\n"
		<< "  edges:      " << n_edges_ << "\n"
		<< "  dofs:       " << n_dofs() << "\n"
		<< "  M:          " << M.rows() << " x " << M.cols()
		<< ", nnz = " << M.nonZeros() << "\n"
		<< "  K_D:        " << H.rows() << " x " << H.cols()
		<< ", nnz = " << H.nonZeros() << "\n";
}

} // namespace femg
