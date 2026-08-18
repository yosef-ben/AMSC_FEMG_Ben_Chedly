#ifndef FEMG_FISHER_KOLMOGOROV_PROBLEM_HPP
#define FEMG_FISHER_KOLMOGOROV_PROBLEM_HPP

#include "quantum_graph_problem.hpp"

#include <Eigen/SparseLU>

#include <array>
#include <cstddef>
#include <functional>
#include <iosfwd>
#include <string>
#include <vector>

namespace femg {

// FEM discretization of the Fisher-Kolmogorov equation on a metric graph.
class fisher_kolmogorov_problem : public quantum_graph_problem {
public:
	enum class TimeScheme {
		corti_semi_implicit,
		backward_euler
	};

	using TimeStepCallback =
		std::function<void(std::size_t, double, const Vector &)>;

	fisher_kolmogorov_problem(double final_time, double time_step);

	void set_edge_diffusion_coefficients(std::vector<double> coefficients);
	void set_vertex_reaction_coefficients(std::vector<double> coefficients);
	void set_vertex_initial_condition(std::vector<double> values);
	void set_vertex_visualization_coordinates(
		std::vector<std::array<double, 3>> coordinates);
	void set_output_enabled(bool enabled);
	void set_output_directory(std::string directory);
	void set_output_stride(std::size_t stride);
	void set_verbose(bool verbose);
	void set_time_scheme(TimeScheme scheme);
	// Row-sum mass lumping. The lumped mass, the reaction vector and the
	// reaction weight matrix are then evaluated with the vertex (trapezoidal)
	// rule, which removes the positive off-diagonal mass entries and restores
	// the discrete maximum principle in the reaction-dominated regime.
	void set_mass_lumping(bool enabled);
	bool mass_lumping() const { return mass_lumping_; }
	void set_newton_parameters(double tolerance, std::size_t max_iterations);
	void set_time_step_callback(TimeStepCallback callback);
	void set_edge_initial_values(
		std::size_t edge_index,
		const std::vector<double> &values);

	void set_coefficients() override;
	void assembly() override;
	void assemble_matrices();
	bool solve() override;
	void sol_export() override;

	const SparseMatrix &mass_matrix() const { return M; }
	const SparseMatrix &diffusion_matrix() const { return H; }
	const Vector &solution() const { return U; }
	const Vector &reaction_coefficients() const { return alpha_dofs_; }
	double vertex_value(std::size_t vertex_index) const;
	std::vector<double> edge_values(std::size_t edge_index) const;
	std::size_t number_of_vertices() const { return n_vertices_orig_; }
	std::size_t number_of_edges() const { return n_edges_; }
	std::size_t number_of_dofs() const { return n_dofs(); }
	double time() const { return time_; }

	Vector assemble_reaction_vector(const Vector &state) const;
	void print_matrix_summary(std::ostream &out) const;

private:
	void validate_input_data() const;
	void interpolate_vertex_data_to_dofs(
		const std::vector<double> &vertex_values,
		Vector &dof_values) const;
	Edge edge_with_index(std::size_t edge_index) const;
	SparseMatrix assemble_reaction_weight_matrix(
		const Vector &state,
		double state_factor) const;
	Vector solve_backward_euler_step(const Vector &old_solution);
	void output(std::size_t step) const;
	void write_pvd_record(const std::vector<std::size_t> &steps) const;

	double final_time_ = 1.0;
	double time_step_ = 0.01;
	double time_ = 0.0;
	bool output_enabled_ = true;
	bool verbose_ = true;
	TimeScheme time_scheme_ = TimeScheme::corti_semi_implicit;
	bool mass_lumping_ = false;
	double newton_tolerance_ = 1.0e-11;
	std::size_t newton_max_iterations_ = 30;
	std::size_t newton_max_backtracks_ = 20;
	double newton_admissible_range_ = 0.5;
	std::size_t output_stride_ = 1;
	std::string output_directory_ = "output/fisher_kolmogorov";

	std::vector<double> edge_diffusion_;
	std::vector<double> vertex_alpha_;
	std::vector<double> vertex_initial_condition_;
	std::vector<std::array<double, 3>> visualization_coordinates_;
	Vector alpha_dofs_;

	SparseMatrix lhs_matrix_;
	SparseMatrix rhs_matrix_;
	Eigen::SparseLU<SparseMatrix> linear_solver_;
	TimeStepCallback time_step_callback_;
};

} // namespace femg

#endif
