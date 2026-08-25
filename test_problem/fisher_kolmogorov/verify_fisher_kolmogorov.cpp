// Verification suite for the Fisher-Kolmogorov metric-graph solver.
//
// Repository-internal correctness checks, independent of the connectome data
// and of the report: exact-solution convergence on an interval and on a
// symmetric star, pure diffusion against a generalized eigendecomposition,
// pure reaction against the logistic solution, a finite-difference check of
// the assembled Newton Jacobian, space/time refinement consistency on a
// nonlinear asymmetric graph, and the discrete invariants (conservation,
// stationary states, symmetry). See verification/README.md for the equations,
// the expectations and the tolerances. The program prints one PASS/FAIL line
// per check, writes machine-readable summaries into the directory given as
// its first argument (default verification/results) and returns a nonzero
// exit code if any check fails.

#include "fisher_kolmogorov_problem.hpp"

#include <Eigen/Dense>
#include <Eigen/Eigenvalues>

#include <cmath>
#include <cstdio>
#include <filesystem>
#include <fstream>
#include <functional>
#include <iostream>
#include <memory>
#include <sstream>
#include <string>
#include <vector>

namespace {

constexpr double kPi = 3.14159265358979323846;

std::string results_directory = "verification/results";
int failures = 0;
std::ofstream summary;

void record(const std::string &test, const std::string &check,
	double value, double bound, bool pass) {
	std::printf("%-4s %-22s %-42s value %-12.4g bound %-10.4g\n",
		pass ? "ok" : "FAIL", test.c_str(), check.c_str(), value, bound);
	summary << test << "," << check << "," << value << "," << bound << ","
		<< (pass ? "pass" : "fail") << "\n";
	if (!pass) {
		++failures;
	}
}

void check_below(const std::string &test, const std::string &check,
	double value, double bound) {
	record(test, check, value, bound, value <= bound);
}

void check_rate(const std::string &test, const std::string &check,
	double rate, double expected, double slack) {
	record(test, check, rate, expected,
		rate >= expected - slack && rate <= expected + slack);
}

// ---------------------------------------------------------------------------
// Graph files, written where the caller asks. Format of the library reader:
// "n_vertices n_edges" then one "u v length n_cells" record per edge.

std::string write_graph(const std::string &name,
	std::size_t n_vertices,
	const std::vector<std::array<double, 4>> &edges) {
	const std::string path = results_directory + "/" + name + ".txt";
	std::ofstream out(path);
	out << n_vertices << " " << edges.size() << "\n";
	for (const auto &edge : edges) {
		out << static_cast<std::size_t>(edge[0]) << " "
			<< static_cast<std::size_t>(edge[1]) << " " << edge[2] << " "
			<< static_cast<int>(edge[3]) << "\n";
	}
	return path;
}

std::string interval_graph(int n_cells) {
	return write_graph("graph_interval_" + std::to_string(n_cells), 2,
		{{0.0, 1.0, 1.0, static_cast<double>(n_cells)}});
}

std::string star_graph(int n_cells) {
	const auto cells = static_cast<double>(n_cells);
	return write_graph("graph_star_" + std::to_string(n_cells), 4,
		{{0.0, 1.0, 1.0, cells}, {0.0, 2.0, 1.0, cells},
		 {0.0, 3.0, 1.0, cells}});
}

// A small graph with a cycle, unequal lengths and a degree-3 vertex; the
// n_cells of every edge are scaled together by the refinement factor.
std::string asymmetric_graph(int refinement) {
	const auto scale = [refinement](double base) {
		return base * static_cast<double>(refinement);
	};
	return write_graph("graph_asym_" + std::to_string(refinement), 4,
		{{0.0, 1.0, 1.0, scale(2)}, {1.0, 2.0, 0.5, scale(1)},
		 {2.0, 3.0, 1.5, scale(3)}, {0.0, 2.0, 1.0, scale(2)},
		 {1.0, 3.0, 2.0, scale(4)}});
}

using Problem = femg::fisher_kolmogorov_problem;

struct Setup {
	std::string graph;
	double final_time = 1.0;
	double time_step = 0.1;
	std::vector<double> diffusion;
	std::vector<double> alpha;
	std::vector<double> vertex_initial;
	Problem::TimeScheme scheme = Problem::TimeScheme::backward_euler;
	bool lumped = false;
};

std::unique_ptr<Problem> build(const Setup &setup) {
	auto problem = std::make_unique<Problem>(setup.final_time, setup.time_step);
	std::string path = setup.graph;
	char program[] = "verify_fisher_kolmogorov";
	char *arguments[] = {program, path.data()};
	problem->init(2, arguments);
	if (!setup.diffusion.empty()) {
		problem->set_edge_diffusion_coefficients(setup.diffusion);
	}
	if (!setup.alpha.empty()) {
		problem->set_vertex_reaction_coefficients(setup.alpha);
	}
	if (!setup.vertex_initial.empty()) {
		problem->set_vertex_initial_condition(setup.vertex_initial);
	}
	problem->set_time_scheme(setup.scheme);
	problem->set_mass_lumping(setup.lumped);
	problem->set_output_enabled(false);
	problem->set_verbose(false);
	problem->set_coefficients();
	return problem;
}

double least_squares_rate(const std::vector<double> &sizes,
	const std::vector<double> &errors) {
	// Slope of log(error) against log(size) over all the levels.
	double sx = 0.0, sy = 0.0, sxx = 0.0, sxy = 0.0;
	const auto n = static_cast<double>(sizes.size());
	for (std::size_t k = 0; k < sizes.size(); ++k) {
		const double x = std::log(sizes[k]);
		const double y = std::log(errors[k]);
		sx += x; sy += y; sxx += x * x; sxy += x * y;
	}
	return (n * sxy - sx * sy) / (n * sxx - sx * sx);
}

// ---------------------------------------------------------------------------
// Test 1: exact cosine mode of the heat equation on the interval and on the
// symmetric star. c(xi, t) = 1/2 + (1/4) exp(-D pi^2 t / L^2) cos(pi xi / L)
// has zero flux at both ends of every unit edge, so it solves the metric
// graph problem with alpha = 0 exactly, on the star as well because the three
// identical edges give zero Kirchhoff sum at the centre.

void set_cosine_initial(Problem &problem, std::size_t n_edges, int n_cells) {
	for (std::size_t edge = 0; edge < n_edges; ++edge) {
		std::vector<double> values(static_cast<std::size_t>(n_cells) + 1);
		for (int node = 0; node <= n_cells; ++node) {
			const double xi = static_cast<double>(node) / n_cells;
			values[static_cast<std::size_t>(node)] =
				0.5 + 0.25 * std::cos(kPi * xi);
		}
		problem.set_edge_initial_values(edge, values);
	}
}

// L2 and H1-seminorm errors of the P1 solution against the exact cosine mode
// at the final time, by 3-point Gauss quadrature on every cell.
std::array<double, 2> cosine_errors(Problem &problem, std::size_t n_edges,
	int n_cells, double time) {
	const double decay = std::exp(-kPi * kPi * time);
	const double root = std::sqrt(3.0 / 5.0);
	const double points[3] = {0.5 * (1.0 - root), 0.5, 0.5 * (1.0 + root)};
	const double weights[3] = {5.0 / 18.0, 4.0 / 9.0, 5.0 / 18.0};
	const double h = 1.0 / n_cells;
	double l2 = 0.0;
	double h1 = 0.0;
	for (std::size_t edge = 0; edge < n_edges; ++edge) {
		const std::vector<double> values = problem.edge_values(edge);
		for (int cell = 0; cell < n_cells; ++cell) {
			const double c0 = values[static_cast<std::size_t>(cell)];
			const double c1 = values[static_cast<std::size_t>(cell) + 1];
			const double slope = (c1 - c0) / h;
			for (int q = 0; q < 3; ++q) {
				const double xi = (cell + points[q]) * h;
				const double exact = 0.5 + 0.25 * decay * std::cos(kPi * xi);
				const double exact_slope =
					-0.25 * kPi * decay * std::sin(kPi * xi);
				const double numeric = c0 + (c1 - c0) * points[q];
				l2 += weights[q] * h * (numeric - exact) * (numeric - exact);
				h1 += weights[q] * h
					* (slope - exact_slope) * (slope - exact_slope);
			}
		}
	}
	return {std::sqrt(l2), std::sqrt(h1)};
}

void test_exact_solution() {
	for (const bool star : {false, true}) {
		const std::string label =
			star ? "1 exact star" : "1 exact interval";
		const std::size_t n_edges = star ? 3 : 1;
		std::vector<double> l2_errors;
		std::vector<double> h1_errors;
		std::vector<double> sizes;
		for (const int n_cells : {4, 8, 16, 32}) {
			Setup setup;
			setup.graph = star ? star_graph(n_cells) : interval_graph(n_cells);
			setup.final_time = 0.02;
			setup.time_step = 1.0e-4;   // temporal error negligible (CN)
			setup.scheme = Problem::TimeScheme::corti_semi_implicit;
			auto problem = build(setup);
			set_cosine_initial(*problem, n_edges, n_cells);
			problem->assemble_matrices();
			problem->solve();
			const auto errors =
				cosine_errors(*problem, n_edges, n_cells, setup.final_time);
			l2_errors.push_back(errors[0]);
			h1_errors.push_back(errors[1]);
			sizes.push_back(1.0 / n_cells);
		}
		check_rate(label, "L2 rate in h (expect 2)",
			least_squares_rate(sizes, l2_errors), 2.0, 0.25);
		check_rate(label, "H1 rate in h (expect 1)",
			least_squares_rate(sizes, h1_errors), 1.0, 0.25);
		// The sharp P1 interpolation bound is (h^2/pi^2) |u|_H2 per edge,
		// with |u|_H2 = 0.25 exp(-pi^2 T) pi^2 sqrt(1/2) = 1.432 at T = 0.02,
		// i.e. 1.42e-4 at h = 1/32, times sqrt(3) on the three-edge star; the
		// finite element solution sits a few percent above its interpolant.
		check_below(label, "L2 error at h = 1/32 (theory 1.4e-4 per edge)",
			l2_errors.back(), (star ? std::sqrt(3.0) : 1.0) * 1.8e-4);
	}

	// Temporal convergence at fixed fine mesh, against the exact solution.
	for (const auto scheme : {Problem::TimeScheme::backward_euler,
			Problem::TimeScheme::corti_semi_implicit}) {
		const bool euler = scheme == Problem::TimeScheme::backward_euler;
		const std::string label = euler
			? "1 time backward-euler" : "1 time crank-nicolson";
		const int n_cells = 64;
		std::vector<double> errors;
		std::vector<double> steps;
		// The spatial error at h = 1/64 is ~1e-5, so the Euler errors stay
		// well above it for these steps; Crank-Nicolson is compared against
		// a small-step run of itself to remove the spatial floor.
		Eigen::VectorXd reference;
		if (!euler) {
			Setup setup;
			setup.graph = interval_graph(n_cells);
			setup.final_time = 0.5;
			setup.time_step = 0.5 / 4096.0;
			setup.scheme = scheme;
			auto problem = build(setup);
			set_cosine_initial(*problem, 1, n_cells);
			problem->assemble_matrices();
			problem->solve();
			reference = problem->solution();
		}
		for (const int divisions : {8, 16, 32, 64}) {
			Setup setup;
			setup.graph = interval_graph(n_cells);
			setup.final_time = 0.5;
			setup.time_step = 0.5 / divisions;
			setup.scheme = scheme;
			auto problem = build(setup);
			set_cosine_initial(*problem, 1, n_cells);
			problem->assemble_matrices();
			problem->solve();
			if (euler) {
				const auto err =
					cosine_errors(*problem, 1, n_cells, setup.final_time);
				errors.push_back(err[0]);
			} else {
				const Eigen::VectorXd difference =
					problem->solution() - reference;
				errors.push_back(std::sqrt(
					difference.dot(problem->mass_matrix() * difference)));
			}
			steps.push_back(setup.time_step);
		}
		check_rate(label, euler ? "L2 rate in dt (expect 1)"
			: "M-norm rate in dt (expect 2)",
			least_squares_rate(steps, errors), euler ? 1.0 : 2.0, 0.25);
	}
}

// ---------------------------------------------------------------------------
// Test 2: pure diffusion on the asymmetric graph against the generalized
// eigendecomposition H v = lambda M v, which gives the exact solution of the
// semi-discrete system M c' + H c = 0.

void test_pure_diffusion() {
	for (const bool lumped : {false, true}) {
		for (const auto scheme : {Problem::TimeScheme::backward_euler,
				Problem::TimeScheme::corti_semi_implicit}) {
			const bool euler = scheme == Problem::TimeScheme::backward_euler;
			std::ostringstream name;
			name << "2 diffusion " << (lumped ? "lumped" : "consistent")
				<< (euler ? " be" : " cn");
			std::vector<double> errors;
			std::vector<double> steps;
			for (const int divisions : {16, 32, 64, 128}) {
				Setup setup;
				setup.graph = asymmetric_graph(2);
				setup.final_time = 1.0;
				setup.time_step = 1.0 / divisions;
				setup.diffusion = {1.0, 0.3, 2.0, 0.7, 1.5};
				setup.vertex_initial = {0.9, 0.1, 0.5, 0.2};
				setup.scheme = scheme;
				setup.lumped = lumped;
				auto problem = build(setup);
				problem->assemble_matrices();

				const Eigen::MatrixXd mass =
					Eigen::MatrixXd(problem->mass_matrix());
				const Eigen::MatrixXd stiffness =
					Eigen::MatrixXd(problem->diffusion_matrix());
				Eigen::GeneralizedSelfAdjointEigenSolver<Eigen::MatrixXd>
					eigen(stiffness, mass);
				const Eigen::VectorXd coefficients =
					eigen.eigenvectors().transpose() * mass
						* problem->solution();
				Eigen::VectorXd reference =
					Eigen::VectorXd::Zero(problem->solution().size());
				for (Eigen::Index k = 0; k < coefficients.size(); ++k) {
					reference += coefficients(k)
						* std::exp(-eigen.eigenvalues()(k) * setup.final_time)
						* eigen.eigenvectors().col(k);
				}

				problem->solve();
				const Eigen::VectorXd difference =
					problem->solution() - reference;
				errors.push_back(std::sqrt(difference.dot(mass * difference)));
				steps.push_back(setup.time_step);
			}
			check_rate(name.str(),
				euler ? "M-norm rate in dt (expect 1)"
					: "M-norm rate in dt (expect 2)",
				least_squares_rate(steps, errors), euler ? 1.0 : 2.0, 0.25);
			check_below(name.str(), "M-norm error at the smallest step",
				errors.back(), euler ? 2.0e-3 : 5.0e-6);
		}
	}
}

// ---------------------------------------------------------------------------
// Test 3: pure reaction. With the lumped mass and the vertex rule the system
// decouples into one logistic equation per degree of freedom, with the exact
// solution c(t) = c0 e^{at} / (1 + c0 (e^{at} - 1)). With the consistent
// mass the reaction couples neighbouring nodes through M^{-1}, so the
// semi-discrete system is nodewise logistic only when the state is uniform;
// the non-uniform consistent case is therefore checked against a small-step
// Runge-Kutta integration of M c' = R(c) built from the same assembled
// reaction vector.

double logistic(double c0, double alpha, double time) {
	const double growth = std::exp(alpha * time);
	return c0 * growth / (1.0 + c0 * (growth - 1.0));
}

void test_pure_reaction() {
	const double alpha = 0.7;
	const double final_time = 4.0;

	const auto run = [&](bool lumped, const std::vector<double> &initial,
			int divisions) {
		Setup setup;
		setup.graph = asymmetric_graph(1);
		setup.final_time = final_time;
		setup.time_step = final_time / divisions;
		setup.diffusion = {0.0, 0.0, 0.0, 0.0, 0.0};
		setup.alpha = {alpha, alpha, alpha, alpha};
		setup.vertex_initial = initial;
		setup.scheme = Problem::TimeScheme::backward_euler;
		setup.lumped = lumped;
		auto problem = build(setup);
		problem->assemble_matrices();
		problem->solve();
		return problem;
	};

	// (a) lumped, non-uniform state: exact nodewise logistic.
	{
		std::vector<double> errors;
		std::vector<double> steps;
		const std::vector<double> initial = {0.05, 0.3, 0.6, 0.85};
		for (const int divisions : {40, 80, 160, 320}) {
			auto problem = run(true, initial, divisions);
			double error = 0.0;
			for (std::size_t vertex = 0; vertex < 4; ++vertex) {
				error = std::max(error, std::abs(problem->vertex_value(vertex)
					- logistic(initial[vertex], alpha, final_time)));
			}
			errors.push_back(error);
			steps.push_back(final_time / divisions);
		}
		check_rate("3 reaction lumped", "max rate in dt (expect 1)",
			least_squares_rate(steps, errors), 1.0, 0.2);
		check_below("3 reaction lumped", "max error at the smallest step",
			errors.back(), 2.0e-3);
	}

	// (b) consistent, uniform state: the reaction vector is alpha c (1-c)
	// times the mass row sums, so the semi-discrete solution is the scalar
	// logistic exactly.
	{
		std::vector<double> errors;
		std::vector<double> steps;
		for (const int divisions : {40, 80, 160, 320}) {
			auto problem = run(false, {0.2, 0.2, 0.2, 0.2}, divisions);
			double error = 0.0;
			for (std::size_t vertex = 0; vertex < 4; ++vertex) {
				error = std::max(error, std::abs(problem->vertex_value(vertex)
					- logistic(0.2, alpha, final_time)));
			}
			errors.push_back(error);
			steps.push_back(final_time / divisions);
		}
		check_rate("3 reaction uniform", "max rate in dt (expect 1)",
			least_squares_rate(steps, errors), 1.0, 0.2);
		check_below("3 reaction uniform", "max error at the smallest step",
			errors.back(), 2.0e-3);
	}

	// (c) consistent, non-uniform state, against RK4 on M c' = R(c).
	{
		const std::vector<double> initial = {0.05, 0.3, 0.6, 0.85};
		auto assembler = run(false, initial, 40);   // reused for M and R
		const Eigen::MatrixXd mass =
			Eigen::MatrixXd(assembler->mass_matrix());
		const Eigen::PartialPivLU<Eigen::MatrixXd> factor(mass);
		Eigen::VectorXd state = [&] {
			auto fresh = run(false, initial, 1);
			return Eigen::VectorXd(fresh->solution());
		}();
		// Rebuild the initial state: run(..., 1) already advanced one step,
		// so take the initial condition from a fresh unsolved problem.
		{
			Setup setup;
			setup.graph = asymmetric_graph(1);
			setup.final_time = final_time;
			setup.time_step = final_time;
			setup.diffusion = {0.0, 0.0, 0.0, 0.0, 0.0};
			setup.alpha = {alpha, alpha, alpha, alpha};
			setup.vertex_initial = initial;
			auto fresh = build(setup);
			fresh->assemble_matrices();
			state = fresh->solution();
		}
		const int reference_steps = 40000;
		const double dt = final_time / reference_steps;
		const auto slope = [&](const Eigen::VectorXd &value) {
			return Eigen::VectorXd(
				factor.solve(assembler->assemble_reaction_vector(value)));
		};
		for (int step = 0; step < reference_steps; ++step) {
			const Eigen::VectorXd k1 = slope(state);
			const Eigen::VectorXd k2 = slope(state + 0.5 * dt * k1);
			const Eigen::VectorXd k3 = slope(state + 0.5 * dt * k2);
			const Eigen::VectorXd k4 = slope(state + dt * k3);
			state += dt / 6.0 * (k1 + 2.0 * k2 + 2.0 * k3 + k4);
		}
		std::vector<double> errors;
		std::vector<double> steps;
		for (const int divisions : {40, 80, 160, 320}) {
			auto problem = run(false, initial, divisions);
			const Eigen::VectorXd difference = problem->solution() - state;
			errors.push_back(difference.cwiseAbs().maxCoeff());
			steps.push_back(final_time / divisions);
		}
		check_rate("3 reaction consistent", "max rate in dt (expect 1)",
			least_squares_rate(steps, errors), 1.0, 0.2);
		check_below("3 reaction consistent", "max error at the smallest step",
			errors.back(), 2.0e-3);
		// The consistent semi-discrete flow is genuinely different from the
		// nodewise logistic: the gap between the two references must be far
		// above the discretization errors, or test (c) would be vacuous.
		double gap = 0.0;
		for (std::size_t vertex = 0; vertex < 4; ++vertex) {
			const auto dof = state.size() - 4 + static_cast<Eigen::Index>(vertex);
			gap = std::max(gap, std::abs(state(dof)
				- logistic(initial[vertex], alpha, final_time)));
		}
		record("3 reaction consistent",
			"consistent flow differs from nodewise logistic", gap, 1.0e-3,
			gap > 1.0e-3);
	}
}

// ---------------------------------------------------------------------------
// Test 4: the assembled Newton Jacobian of the backward Euler residual
// F(c) = M (c - c_old) + dt H c - dt R(c) against central finite
// differences of F, over a range of increments.

void test_jacobian() {
	for (const bool lumped : {false, true}) {
		const std::string label = lumped
			? "4 jacobian lumped" : "4 jacobian consistent";
		Setup setup;
		setup.graph = asymmetric_graph(2);
		setup.final_time = 1.0;
		setup.time_step = 0.25;
		setup.diffusion = {1.0, 0.3, 2.0, 0.7, 1.5};
		setup.alpha = {0.8, 0.5, 0.3, 0.6};
		setup.vertex_initial = {0.9, 0.1, 0.5, 0.2};
		setup.lumped = lumped;
		auto problem = build(setup);
		problem->assemble_matrices();
		const double dt = setup.time_step;
		const auto size = problem->solution().size();
		const Eigen::VectorXd old_state = problem->solution();

		const auto residual = [&](const Eigen::VectorXd &state) {
			return Eigen::VectorXd(
				problem->mass_matrix() * (state - old_state)
				+ dt * (problem->diffusion_matrix() * state)
				- dt * problem->assemble_reaction_vector(state));
		};

		unsigned seed = 12345;
		const auto uniform = [&seed]() {
			seed = 1664525u * seed + 1013904223u;
			return static_cast<double>(seed) / 4294967296.0;
		};

		double worst_best = 0.0;
		for (int trial = 0; trial < 5; ++trial) {
			Eigen::VectorXd state(size);
			Eigen::VectorXd direction(size);
			for (Eigen::Index i = 0; i < size; ++i) {
				state(i) = 0.05 + 0.9 * uniform();
				direction(i) = 2.0 * uniform() - 1.0;
			}
			direction /= direction.norm();
			const femg::SparseMatrix jacobian =
				femg::SparseMatrix(problem->mass_matrix())
				+ dt * femg::SparseMatrix(problem->diffusion_matrix())
				- dt * problem->assemble_reaction_weight_matrix(state, 2.0);
			const Eigen::VectorXd analytic = jacobian * direction;
			double best = 1.0e30;
			for (double eps = 1.0e-3; eps >= 1.0e-9; eps *= 0.1) {
				const Eigen::VectorXd numeric =
					(residual(state + eps * direction)
					- residual(state - eps * direction)) / (2.0 * eps);
				best = std::min(best,
					(numeric - analytic).norm()
						/ std::max(analytic.norm(), 1.0e-30));
			}
			worst_best = std::max(worst_best, best);
		}
		check_below(label,
			"worst directional relative error over 5 trials",
			worst_best, 1.0e-9);
	}
}

// ---------------------------------------------------------------------------
// Test 5: space/time refinement consistency on a nonlinear problem coupling
// diffusion and reaction on the asymmetric graph.

void test_refinement_consistency() {
	const std::vector<double> alpha = {0.8, 0.5, 0.3, 0.6};
	const std::vector<double> diffusion = {0.4, 0.12, 0.8, 0.28, 0.6};
	const std::vector<double> initial = {0.4, 0.02, 0.02, 0.02};

	const auto vertex_state = [&](int refinement, int divisions) {
		Setup setup;
		setup.graph = asymmetric_graph(refinement);
		setup.final_time = 8.0;
		setup.time_step = 8.0 / divisions;
		setup.diffusion = diffusion;
		setup.alpha = alpha;
		setup.vertex_initial = initial;
		auto problem = build(setup);
		problem->assemble_matrices();
		problem->solve();
		Eigen::Vector4d values;
		for (std::size_t vertex = 0; vertex < 4; ++vertex) {
			values(static_cast<Eigen::Index>(vertex)) =
				problem->vertex_value(vertex);
		}
		return values;
	};

	// Refine h at small dt: successive differences shrink at the P1 rate.
	{
		std::vector<Eigen::Vector4d> states;
		for (const int refinement : {1, 2, 4, 8}) {
			states.push_back(vertex_state(refinement, 3200));
		}
		std::vector<double> differences;
		std::vector<double> sizes;
		for (std::size_t level = 0; level + 1 < states.size(); ++level) {
			differences.push_back(
				(states[level] - states[level + 1]).cwiseAbs().maxCoeff());
			sizes.push_back(1.0 / static_cast<double>(1 << level));
		}
		check_rate("5 h refinement", "vertex rate in h (expect 2)",
			least_squares_rate(sizes, differences), 2.0, 0.4);
		check_below("5 h refinement", "difference at the finest pair",
			differences.back(), 1.0e-4);
	}

	// Refine dt at the finest mesh: backward Euler rate one.
	{
		std::vector<Eigen::Vector4d> states;
		std::vector<double> steps;
		for (const int divisions : {50, 100, 200, 400, 800}) {
			states.push_back(vertex_state(8, divisions));
			steps.push_back(8.0 / divisions);
		}
		std::vector<double> differences;
		for (std::size_t level = 0; level + 1 < states.size(); ++level) {
			differences.push_back(
				(states[level] - states[level + 1]).cwiseAbs().maxCoeff());
		}
		steps.pop_back();
		check_rate("5 dt refinement", "vertex rate in dt (expect 1)",
			least_squares_rate(steps, differences), 1.0, 0.2);
	}
}

// ---------------------------------------------------------------------------
// Test 6: discrete invariants.

void test_invariants() {
	// (a) alpha = 0: the total metric mass 1^T M c is conserved by both
	// schemes, because the row sums of H vanish.
	for (const bool lumped : {false, true}) {
		for (const auto scheme : {Problem::TimeScheme::backward_euler,
				Problem::TimeScheme::corti_semi_implicit}) {
			Setup setup;
			setup.graph = asymmetric_graph(2);
			setup.final_time = 2.0;
			setup.time_step = 0.05;
			setup.diffusion = {1.0, 0.3, 2.0, 0.7, 1.5};
			setup.vertex_initial = {0.9, 0.1, 0.5, 0.2};
			setup.scheme = scheme;
			setup.lumped = lumped;
			auto problem = build(setup);
			problem->assemble_matrices();
			const Eigen::VectorXd ones =
				Eigen::VectorXd::Ones(problem->solution().size());
			const double before =
				ones.dot(problem->mass_matrix() * problem->solution());
			problem->solve();
			const double after =
				ones.dot(problem->mass_matrix() * problem->solution());
			std::ostringstream name;
			name << "6 conservation " << (lumped ? "lumped" : "consistent")
				<< (scheme == Problem::TimeScheme::backward_euler
					? " be" : " cn");
			check_below(name.str(), "relative drift of 1'Mc",
				std::abs(after - before) / std::abs(before), 1.0e-12);
		}
	}

	// (b) stationary states: a uniform state with alpha = 0, and the two
	// equilibria c = 0 and c = 1 of the reaction.
	const auto stationary = [&](const std::string &check, double level,
			double alpha) {
		Setup setup;
		setup.graph = asymmetric_graph(2);
		setup.final_time = 2.0;
		setup.time_step = 0.1;
		setup.diffusion = {1.0, 0.3, 2.0, 0.7, 1.5};
		setup.alpha = {alpha, alpha, alpha, alpha};
		setup.vertex_initial = {level, level, level, level};
		auto problem = build(setup);
		problem->assemble_matrices();
		problem->solve();
		const Eigen::VectorXd target = Eigen::VectorXd::Constant(
			problem->solution().size(), level);
		check_below("6 stationary", check,
			(problem->solution() - target).cwiseAbs().maxCoeff(), 1.0e-12);
	};
	stationary("uniform state with alpha = 0 stays", 0.37, 0.0);
	stationary("c = 0 equilibrium stays", 0.0, 0.9);
	stationary("c = 1 equilibrium stays", 1.0, 0.9);

	// (c) symmetry: on the star with symmetric data the three edges carry
	// the same solution at all times.
	{
		Setup setup;
		setup.graph = star_graph(8);
		setup.final_time = 3.0;
		setup.time_step = 0.05;
		setup.diffusion = {0.5, 0.5, 0.5};
		setup.alpha = {0.6, 0.6, 0.6, 0.6};
		setup.vertex_initial = {0.5, 0.05, 0.05, 0.05};
		auto problem = build(setup);
		problem->assemble_matrices();
		problem->solve();
		const std::vector<double> first = problem->edge_values(0);
		double asymmetry = 0.0;
		for (std::size_t edge = 1; edge < 3; ++edge) {
			const std::vector<double> other = problem->edge_values(edge);
			for (std::size_t node = 0; node < first.size(); ++node) {
				asymmetry = std::max(asymmetry,
					std::abs(other[node] - first[node]));
			}
		}
		check_below("6 symmetry", "max deviation between the star edges",
			asymmetry, 1.0e-12);
	}
}

} // namespace

int main(int argc, char *argv[]) {
	if (argc >= 2) {
		results_directory = argv[1];
	}
	std::filesystem::create_directories(results_directory);
	summary.open(results_directory + "/verification_summary.csv");
	summary << "test,check,value,bound,status\n";

	try {
		test_exact_solution();
		test_pure_diffusion();
		test_pure_reaction();
		test_jacobian();
		test_refinement_consistency();
		test_invariants();
	} catch (const std::exception &error) {
		std::cerr << "Error: " << error.what() << "\n";
		summary << "abort,exception,0,0,fail\n";
		return 2;
	}

	std::printf("\n%d check(s) failed.\n", failures);
	return failures == 0 ? 0 : 1;
}
