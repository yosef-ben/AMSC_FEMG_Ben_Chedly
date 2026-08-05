#ifndef FEMG_SPECTRAL_PROBLEM_HPP
#define FEMG_SPECTRAL_PROBLEM_HPP

#include "quantum_graph_problem.hpp"

#include <cstddef>
#include <iosfwd>
#include <string>
#include <vector>

namespace femg {

class spectral_problem : public quantum_graph_problem {
public:
	explicit spectral_problem(std::size_t n_modes);

	void set_output_directory(std::string output_directory);
	void set_eigenfunction_scale(double scale);

	void set_coefficients() override;
	void assembly() override;
	bool solve() override;
	void sol_export() override;

	double eigenvalue(std::size_t mode) const;
	const Vector &eigenvector(std::size_t mode) const;
	std::size_t dofs() const { return n_dofs(); }
	std::vector<std::vector<double>> edge_nodal_values(
		std::size_t mode) const;
	std::vector<double> edge_lengths() const;
	Vector combinatorial_laplacian_eigenvalues() const;

	void write_eigenvalues_csv(const std::string &filename) const;
	void print_eigenvalues(std::ostream &out) const;

private:
	void assemble_matrices();
	void write_flat_domain() const;
	void write_eigenfunction(std::size_t mode) const;

	std::size_t n_modes_ = 6;
	double eigenfunction_scale_ = 1.0;
	std::string output_directory_ = "output/spectral";

	Vector eigenvalues_;
	std::vector<Vector> eigenvectors_;
};

} // namespace femg

#endif
