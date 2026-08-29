// Sequential static condensation of the semi-implicit time step.
//
// The DoF numbering of the library makes the time-step matrix a bordered
// arrow matrix: one tridiagonal block per edge chain (the interior nodes)
// bordered by the original graph vertices. This stepper advances the
// semi-implicit scheme of the production solver,
//
//     [M + (dt/2) H - (dt/2) W(c_hat)] c_new
//         = [M - (dt/2) H + (dt/2) W(c_hat)] c_old,
//     c_hat = 1.5 c_n - 0.5 c_{n-1}  (c_hat = c_0 at the first step),
//
// without ever assembling a global matrix. Every cell contributes a 2x2
// block (consistent mass, 3-point Gauss on the reaction weight, exactly the
// quadrature of the library); the blocks are accumulated per edge into the
// local tridiagonal system, the interior unknowns are eliminated by a
// Thomas factorization, and each edge leaves a 2x2 contribution to the
// Schur complement on the original vertices plus a condensed right-hand
// side. The vertex system (83 unknowns on the connectome, adjacency
// sparsity, factorized dense) is solved once, and the interiors are
// recovered by one saxpy per edge from the vectors stored during the
// elimination:
//
//     u_int = T^{-1} r_int - u_source T^{-1}(c_s e_1) - u_target T^{-1}(c_t e_m).
//
// The per-edge loop is the natural parallel unit of the later stages: it
// reads the state, writes only edge-local buffers and accumulates into the
// vertex system, which a thread-local copy makes race-free.

#ifndef FEMG_TEST_CONDENSED_STEPPER_HPP
#define FEMG_TEST_CONDENSED_STEPPER_HPP

#include <Eigen/Cholesky>
#include <Eigen/Dense>

#include <chrono>

#include <cmath>
#include <cstddef>
#include <stdexcept>
#include <vector>

namespace condensed {

struct EdgeData {
	int source = 0;             // original vertex index
	int target = 0;
	int cells = 1;
	double h = 0.0;             // cell length
	double diffusion = 0.0;     // edge diffusivity D_e
	std::size_t offset = 0;     // first interior DoF of the edge
};

struct StepTimes {
	double local = 0.0;         // per-edge assembly, elimination, Schur
	double interface = 0.0;     // vertex system assembly close + solve
	double back = 0.0;          // interior back-substitution
};

class Stepper {
public:
	Stepper(std::vector<EdgeData> edges, std::size_t n_vertices,
		std::vector<double> alpha_dofs, double time_step)
		: edges_(std::move(edges)), n_vertices_(n_vertices),
		  alpha_(std::move(alpha_dofs)), dt_(time_step) {
		std::size_t interiors = 0;
		for (const EdgeData &edge : edges_) {
			if (edge.offset != interiors) {
				throw std::invalid_argument("Edge offsets are not packed.");
			}
			interiors += static_cast<std::size_t>(edge.cells - 1);
		}
		vertex_offset_ = interiors;
		n_dofs_ = interiors + n_vertices_;
		if (alpha_.size() != n_dofs_) {
			throw std::invalid_argument("One alpha per DoF is required.");
		}
		diagonal_.resize(interiors);
		lower_.resize(interiors);
		forward_.resize(interiors);
		source_column_.resize(interiors);
		target_column_.resize(interiors);
		schur_.resize(n_vertices_, n_vertices_);
		condensed_rhs_.resize(n_vertices_);
	}

	std::size_t n_dofs() const { return n_dofs_; }
	std::size_t vertex_offset() const { return vertex_offset_; }

	// One semi-implicit step: state holds every DoF and is advanced in
	// place; extrapolated is the reaction state c_hat.
	StepTimes step(const std::vector<double> &extrapolated,
		std::vector<double> &state);

private:
	// The 2x2 cell blocks of lhs = m + (dt/2) k - (dt/2) w and
	// rhs = m - (dt/2) k + (dt/2) w, with the 3-point Gauss reaction
	// weight of the library on the extrapolated state.
	void cell_blocks(double h, double diffusion,
		double alpha0, double alpha1, double chat0, double chat1,
		double lhs[2][2], double rhs[2][2]) const;

	std::vector<EdgeData> edges_;
	std::size_t n_vertices_ = 0;
	std::vector<double> alpha_;
	double dt_ = 0.0;
	std::size_t vertex_offset_ = 0;
	std::size_t n_dofs_ = 0;

	// Per-interior work vectors, packed like the DoFs.
	std::vector<double> diagonal_;       // Thomas pivots
	std::vector<double> lower_;          // subdiagonal of the chain
	std::vector<double> forward_;        // T^{-1} r_int
	std::vector<double> source_column_;  // T^{-1} (c_s e_1)
	std::vector<double> target_column_;  // T^{-1} (c_t e_m)
	Eigen::MatrixXd schur_;
	Eigen::VectorXd condensed_rhs_;
};

inline void Stepper::cell_blocks(double h, double diffusion,
	double alpha0, double alpha1, double chat0, double chat1,
	double lhs[2][2], double rhs[2][2]) const {
	static const double root = std::sqrt(3.0 / 5.0);
	static const double points[3] = {
		0.5 * (1.0 - root), 0.5, 0.5 * (1.0 + root)};
	static const double weights[3] = {5.0 / 18.0, 4.0 / 9.0, 5.0 / 18.0};

	const double mass[2][2] = {{h / 3.0, h / 6.0}, {h / 6.0, h / 3.0}};
	const double k = diffusion / h;
	double reaction[2][2] = {};
	for (int q = 0; q < 3; ++q) {
		const double basis[2] = {1.0 - points[q], points[q]};
		const double concentration = basis[0] * chat0 + basis[1] * chat1;
		const double alpha = basis[0] * alpha0 + basis[1] * alpha1;
		const double coefficient = alpha * (1.0 - concentration);
		for (int i = 0; i < 2; ++i) {
			for (int j = 0; j < 2; ++j) {
				reaction[i][j] +=
					coefficient * basis[i] * basis[j] * weights[q] * h;
			}
		}
	}
	const double stiffness[2][2] = {{k, -k}, {-k, k}};
	for (int i = 0; i < 2; ++i) {
		for (int j = 0; j < 2; ++j) {
			lhs[i][j] = mass[i][j] + 0.5 * dt_ * stiffness[i][j]
				- 0.5 * dt_ * reaction[i][j];
			rhs[i][j] = mass[i][j] - 0.5 * dt_ * stiffness[i][j]
				+ 0.5 * dt_ * reaction[i][j];
		}
	}
}

inline StepTimes Stepper::step(const std::vector<double> &extrapolated,
	std::vector<double> &state) {
	StepTimes times;
	const auto tick = []() {
		return std::chrono::steady_clock::now();
	};
	const auto seconds = [](auto start) {
		return std::chrono::duration<double>(
			std::chrono::steady_clock::now() - start).count();
	};

	auto start = tick();
	schur_.setZero();
	condensed_rhs_.setZero();

	for (const EdgeData &edge : edges_) {
		const std::size_t m = static_cast<std::size_t>(edge.cells - 1);
		const std::size_t source_dof = vertex_offset_
			+ static_cast<std::size_t>(edge.source);
		const std::size_t target_dof = vertex_offset_
			+ static_cast<std::size_t>(edge.target);

		// Local assembly: chain of edge.cells 2x2 blocks over the DoFs
		// [source, interior 0 .. m-1, target].
		double diag_source = 0.0;
		double diag_target = 0.0;
		double couple_source = 0.0;   // source to first interior
		double couple_target = 0.0;   // last interior to target
		double direct = 0.0;          // source to target (one-cell edge)
		double rhs_source = 0.0;
		double rhs_target = 0.0;
		double *diagonal = diagonal_.data() + edge.offset;
		double *lower = lower_.data() + edge.offset;
		double *rhs_int = forward_.data() + edge.offset;
		for (std::size_t k = 0; k < m; ++k) {
			diagonal[k] = 0.0;
			lower[k] = 0.0;
			rhs_int[k] = 0.0;
		}
		for (int cell = 0; cell < edge.cells; ++cell) {
			const std::size_t dof0 = (cell == 0)
				? source_dof
				: edge.offset + static_cast<std::size_t>(cell - 1);
			const std::size_t dof1 = (cell == edge.cells - 1)
				? target_dof
				: edge.offset + static_cast<std::size_t>(cell);
			double lhs[2][2];
			double rhs[2][2];
			cell_blocks(edge.h, edge.diffusion,
				alpha_[dof0], alpha_[dof1],
				extrapolated[dof0], extrapolated[dof1], lhs, rhs);
			const double r0 = rhs[0][0] * state[dof0]
				+ rhs[0][1] * state[dof1];
			const double r1 = rhs[1][0] * state[dof0]
				+ rhs[1][1] * state[dof1];

			if (cell == 0) {
				diag_source += lhs[0][0];
				rhs_source += r0;
			} else {
				diagonal[cell - 1] += lhs[0][0];
				rhs_int[cell - 1] += r0;
			}
			if (cell == edge.cells - 1) {
				diag_target += lhs[1][1];
				rhs_target += r1;
			} else {
				diagonal[cell] += lhs[1][1];
				rhs_int[cell] += r1;
			}
			if (edge.cells == 1) {
				direct = lhs[0][1];
			} else if (cell == 0) {
				couple_source = lhs[0][1];
			} else if (cell == edge.cells - 1) {
				couple_target = lhs[0][1];
			} else {
				lower[cell - 1] = lhs[0][1];
			}
		}

		double s_ss = diag_source;
		double s_tt = diag_target;
		double s_st = direct;
		double r_s = rhs_source;
		double r_t = rhs_target;
		if (m > 0) {
			// Thomas factorization of the interior chain and the three
			// forward solves (rhs, source column, target column) fused.
			double *y = rhs_int;
			double *zs = source_column_.data() + edge.offset;
			double *zt = target_column_.data() + edge.offset;
			zs[0] = couple_source;
			zt[0] = 0.0;
			for (std::size_t k = 1; k < m; ++k) {
				zs[k] = 0.0;
				zt[k] = 0.0;
			}
			zt[m - 1] = couple_target;
			for (std::size_t k = 1; k < m; ++k) {
				const double factor = lower[k - 1] / diagonal[k - 1];
				diagonal[k] -= factor * lower[k - 1];
				y[k] -= factor * y[k - 1];
				zs[k] -= factor * zs[k - 1];
				zt[k] -= factor * zt[k - 1];
			}
			y[m - 1] /= diagonal[m - 1];
			zs[m - 1] /= diagonal[m - 1];
			zt[m - 1] /= diagonal[m - 1];
			for (std::size_t k = m - 1; k-- > 0;) {
				const double factor = lower[k];
				y[k] = (y[k] - factor * y[k + 1]) / diagonal[k];
				zs[k] = (zs[k] - factor * zs[k + 1]) / diagonal[k];
				zt[k] = (zt[k] - factor * zt[k + 1]) / diagonal[k];
			}
			s_ss -= couple_source * zs[0];
			s_st = -couple_source * zt[0];
			s_tt -= couple_target * zt[m - 1];
			r_s -= couple_source * y[0];
			r_t -= couple_target * y[m - 1];
		}
		schur_(edge.source, edge.source) += s_ss;
		schur_(edge.target, edge.target) += s_tt;
		schur_(edge.source, edge.target) += s_st;
		schur_(edge.target, edge.source) += s_st;
		condensed_rhs_(edge.source) += r_s;
		condensed_rhs_(edge.target) += r_t;
	}
	times.local = seconds(start);

	start = tick();
	const Eigen::VectorXd vertices =
		Eigen::LDLT<Eigen::MatrixXd>(schur_).solve(condensed_rhs_);
	for (std::size_t vertex = 0; vertex < n_vertices_; ++vertex) {
		state[vertex_offset_ + vertex] =
			vertices(static_cast<Eigen::Index>(vertex));
	}
	times.interface = seconds(start);

	start = tick();
	for (const EdgeData &edge : edges_) {
		const std::size_t m = static_cast<std::size_t>(edge.cells - 1);
		if (m == 0) {
			continue;
		}
		const double u_source =
			vertices(static_cast<Eigen::Index>(edge.source));
		const double u_target =
			vertices(static_cast<Eigen::Index>(edge.target));
		const double *y = forward_.data() + edge.offset;
		const double *zs = source_column_.data() + edge.offset;
		const double *zt = target_column_.data() + edge.offset;
		double *interior = state.data() + edge.offset;
		for (std::size_t k = 0; k < m; ++k) {
			interior[k] = y[k] - u_source * zs[k] - u_target * zt[k];
		}
	}
	times.back = seconds(start);
	return times;
}

} // namespace condensed

#endif
