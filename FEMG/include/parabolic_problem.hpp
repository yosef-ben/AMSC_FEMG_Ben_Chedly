#ifndef FEMG_PARABOLIC_PROBLEM_HPP
#define FEMG_PARABOLIC_PROBLEM_HPP

#include "quantum_graph_problem.hpp"

#include <Eigen/SparseLU>

#include <cmath>
#include <functional>
#include <iosfwd>
#include <memory>
#include <string>
#include <vector>

namespace femg {

class parabolic_problem : public quantum_graph_problem {
public:
	using EdgeFunction = std::function<double(std::size_t, double, double)>;
	using VertexFunction = std::function<double(std::size_t, double)>;
	using TimeStepCallback =
		std::function<void(std::size_t, double, const Vector &)>;

	class ForcingTerm {
	public:
		virtual ~ForcingTerm() = default;

		virtual double value(
			std::size_t edge_index,
			double s,
			double time) const {
			return 0.0;
		}
	};

	class FunctionU0 {
	public:
		virtual ~FunctionU0() = default;

		virtual double value(
			std::size_t edge_index,
			double s,
			double time) const {
			return std::cos(M_PI * s);
		}
	};

	class ExactSolution {
	public:
		virtual ~ExactSolution() = default;

		virtual double value(
			std::size_t edge_index,
			double s,
			double time) const {
			return std::exp(-M_PI * M_PI * time) * std::cos(M_PI * s);
		}
	};

	struct TimingData {
		double factorization_seconds = 0.0;
		double time_stepping_seconds = 0.0;
		double output_seconds = 0.0;
		double total_solve_seconds = 0.0;
		std::size_t n_time_steps = 0;
	};

	parabolic_problem(
		double diffusion,
		double final_time,
		double time_step,
		double theta);

	void set_source_function(EdgeFunction source_function);
	void set_initial_condition(EdgeFunction initial_condition);
	void set_source_function(std::shared_ptr<const ForcingTerm> source_function);
	void set_initial_condition(std::shared_ptr<const FunctionU0> initial_condition);
	void set_reaction_coefficient(double reaction);
	void set_output_enabled(bool output_enabled);
	void set_output_directory(std::string output_directory);
	void set_verbose(bool verbose);
	void set_time_step_callback(TimeStepCallback callback);
	void set_dirichlet_vertices(
		std::vector<std::size_t> vertices,
		VertexFunction boundary_value);
	void clear_dirichlet_vertices();

	void set_coefficients() override;
	void assembly() override;
	void assemble_matrices();
	bool solve() override;
	TimingData solve_with_timing();
	void sol_export() override;

	void print_matrix_summary(std::ostream &out) const;

	const SparseMatrix &mass_matrix() const { return M; }
	const SparseMatrix &stiffness_matrix() const { return H; }
	const SparseMatrix &lhs_matrix() const { return lhs_matrix_; }
	const SparseMatrix &rhs_matrix() const { return rhs_matrix_; }
	const Vector &solution() const { return U; }
	Vector assemble_load_vector(double time) const;
	double compute_l2_error(const ExactSolution &exact_solution) const;

private:
	void build_theta_matrices();
	void assemble_rhs(double new_time);
	void solve_time_step();
	void output(std::size_t time_step_number) const;
	void write_pvd_record(std::size_t last_time_step_number) const;
	void impose_dirichlet_on_lhs();
	void impose_dirichlet_on_rhs(double time);
	bool is_dirichlet_dof(Eigen::Index dof) const;

	double diffusion_ = 1.0;
	double final_time_ = 1.0;
	double time_step_ = 0.01;
	double theta_ = 1.0;
	double reaction_ = 0.0;
	double time_ = 0.0;
	bool output_enabled_ = true;
	bool verbose_ = true;
	std::string output_directory_ = "output/visualization";

	SparseMatrix lhs_matrix_;
	SparseMatrix rhs_matrix_;
	Eigen::SparseLU<SparseMatrix> linear_solver_;
	Vector system_rhs_;
	EdgeFunction source_function_;
	EdgeFunction initial_condition_;
	std::shared_ptr<const ForcingTerm> forcing_term_;
	std::shared_ptr<const FunctionU0> initial_condition_term_;
	std::vector<std::size_t> dirichlet_vertices_;
	VertexFunction dirichlet_value_;
	TimeStepCallback time_step_callback_;
};

} // namespace femg

#endif
