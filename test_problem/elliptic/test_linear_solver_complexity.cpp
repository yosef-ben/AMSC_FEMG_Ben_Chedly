#include "parabolic_problem.hpp"

#include <Eigen/IterativeLinearSolvers>
#include <Eigen/SparseLU>
#include <unsupported/Eigen/IterativeSolvers>

#include <chrono>
#include <cmath>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <stdexcept>
#include <string>
#include <vector>

namespace {

using femg::SparseMatrix;
using femg::Vector;

std::string write_graphene_file(int n_cells) {
	std::filesystem::create_directories("output/timing/graphs");

	const std::string filename =
		"output/timing/graphs/graphene_linear_" + std::to_string(n_cells) + ".txt";

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

SparseMatrix assemble_elliptic_matrix(int n_cells, int &n_dofs) {
	const std::string graph_file = write_graphene_file(n_cells);
	std::string executable_name = "test_linear_solver_complexity";
	std::string local_graph_file = graph_file;
	char *argv[] = {executable_name.data(), local_graph_file.data()};

	femg::parabolic_problem problem(1.0, 0.1, 0.01, 0.5);
	problem.set_output_enabled(false);
	problem.set_verbose(false);
	problem.init(2, argv);
	problem.set_coefficients();
	problem.assemble_matrices();

	n_dofs = static_cast<int>(problem.solution().size());
	SparseMatrix matrix = problem.stiffness_matrix() + problem.mass_matrix();
	matrix.makeCompressed();
	return matrix;
}

Vector build_rhs(Eigen::Index size) {
	Vector rhs(size);
	for (Eigen::Index i = 0; i < size; ++i) {
		const double x = static_cast<double>(i + 1);
		rhs(i) = 1.0 + std::sin(0.013 * x) + 0.5 * std::cos(0.031 * x);
	}
	return rhs;
}

struct SolverTiming {
	double seconds = 0.0;
	int iterations = 0;
	double error = 0.0;
	int info = 0;
};

template <typename Solver>
SolverTiming run_iterative_solver(
	Solver &solver,
	const SparseMatrix &matrix,
	const Vector &rhs) {
	const auto start = std::chrono::steady_clock::now();
	solver.setTolerance(1.0e-8);
	solver.setMaxIterations(static_cast<int>(std::max<Eigen::Index>(1000, 2 * matrix.rows())));
	solver.compute(matrix);
	const Vector solution = solver.solve(rhs);
	const auto end = std::chrono::steady_clock::now();

	SolverTiming timing;
	timing.seconds = std::chrono::duration<double>(end - start).count();
	timing.iterations = static_cast<int>(solver.iterations());
	timing.error = solver.error();
	timing.info = static_cast<int>(solver.info());

	const double residual =
		(matrix * solution - rhs).norm() / std::max(rhs.norm(), 1.0e-16);
	if (residual > timing.error) {
		timing.error = residual;
	}

	return timing;
}

SolverTiming run_lu_solver(const SparseMatrix &matrix, const Vector &rhs) {
	const auto start = std::chrono::steady_clock::now();
	Eigen::SparseLU<SparseMatrix> solver;
	solver.compute(matrix);
	const Vector solution = solver.solve(rhs);
	const auto end = std::chrono::steady_clock::now();

	SolverTiming timing;
	timing.seconds = std::chrono::duration<double>(end - start).count();
	timing.iterations = 1;
	timing.info = static_cast<int>(solver.info());
	timing.error =
		(matrix * solution - rhs).norm() / std::max(rhs.norm(), 1.0e-16);
	return timing;
}

struct TimingRow {
	int n_cells = 0;
	int n_dofs = 0;
	int repetitions = 0;
	SolverTiming cg;
	SolverTiming gmres;
	SolverTiming bicgstab;
	SolverTiming lu;
};

void add_timing(SolverTiming &average, const SolverTiming &current) {
	average.seconds += current.seconds;
	average.iterations += current.iterations;
	average.error += current.error;
	average.info = std::max(average.info, current.info);
}

void scale_timing(SolverTiming &timing, double scale) {
	timing.seconds *= scale;
	timing.iterations = static_cast<int>(std::round(timing.iterations * scale));
	timing.error *= scale;
}

TimingRow run_case(int n_cells, int repetitions) {
	int n_dofs = 0;
	const SparseMatrix matrix = assemble_elliptic_matrix(n_cells, n_dofs);
	const Vector rhs = build_rhs(matrix.rows());

	TimingRow row;
	row.n_cells = n_cells;
	row.n_dofs = n_dofs;
	row.repetitions = repetitions;

	for (int rep = 0; rep < repetitions; ++rep) {
		Eigen::ConjugateGradient<SparseMatrix, Eigen::Lower | Eigen::Upper, Eigen::IncompleteCholesky<double>> cg;
		Eigen::GMRES<SparseMatrix, Eigen::IncompleteLUT<double>> gmres;
		Eigen::BiCGSTAB<SparseMatrix, Eigen::IncompleteLUT<double>> bicgstab;

		gmres.set_restart(50);

		add_timing(row.cg, run_iterative_solver(cg, matrix, rhs));
		add_timing(row.gmres, run_iterative_solver(gmres, matrix, rhs));
		add_timing(row.bicgstab, run_iterative_solver(bicgstab, matrix, rhs));
		add_timing(row.lu, run_lu_solver(matrix, rhs));
	}

	const double scale = 1.0 / static_cast<double>(repetitions);
	scale_timing(row.cg, scale);
	scale_timing(row.gmres, scale);
	scale_timing(row.bicgstab, scale);
	scale_timing(row.lu, scale);

	return row;
}

void write_timing(
	std::ofstream &out,
	const TimingRow &row,
	const std::string &method,
	const SolverTiming &timing) {
	out << row.n_cells << ","
			<< row.n_dofs << ","
			<< row.repetitions << ","
			<< method << ","
			<< std::setprecision(16)
			<< timing.seconds << ","
			<< timing.iterations << ","
			<< timing.error << ","
			<< timing.info << "\n";
}

} // namespace

int main() {
	try {
		std::filesystem::create_directories("output/timing");

		const std::vector<int> cell_counts = {20, 40, 80, 160, 320};
		const int repetitions = 1;

		std::ofstream out("output/timing/graphene_linear_solver_complexity.csv");
		out << "n_cells_per_edge,n_dofs,repetitions,method,seconds,"
				<< "iterations,relative_residual,info\n";

		for (const int n_cells : cell_counts) {
			const TimingRow row = run_case(n_cells, repetitions);
			write_timing(out, row, "CG", row.cg);
			write_timing(out, row, "GMRES", row.gmres);
			write_timing(out, row, "BiCGSTAB", row.bicgstab);
			write_timing(out, row, "LU", row.lu);

			std::cout << "n_cells = " << row.n_cells
								<< ", dofs = " << row.n_dofs
								<< ", CG = " << row.cg.seconds
								<< " s, GMRES = " << row.gmres.seconds
								<< " s, BiCGSTAB = " << row.bicgstab.seconds
								<< " s, LU = " << row.lu.seconds << " s\n";
		}

		std::cout << "Wrote output/timing/graphene_linear_solver_complexity.csv\n";
		return 0;
	} catch (const std::exception &error) {
		std::cerr << "Error: " << error.what() << "\n";
		return 1;
	}
}
