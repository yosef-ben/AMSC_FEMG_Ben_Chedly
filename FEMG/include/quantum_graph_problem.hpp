// Core base class for quantum-graph PDE problems.
#ifndef FEMG_QUANTUM_GRAPH_PROBLEM_HPP
#define FEMG_QUANTUM_GRAPH_PROBLEM_HPP

#include <boost/graph/adjacency_list.hpp> // Boost Graph Library

#include <Eigen/Sparse> // Eigen for sparse linear algebra
#include <Eigen/Dense> // Eigen for dense linear algebra

#include <array> 
#include <cstddef>
#include <string>
#include <vector>

namespace femg {

// Vertex and edge properties for the metric graph.
struct VertexProperty {
	std::array<double, 2> coords{0.0, 0.0};
	std::size_t index = 0; // original vertex index in input (up to N-1)
};

struct EdgeProperty {
	double length = 0.0; // length of the edge l_e
	int n_cells = 0;       // number of sub-intervals on the edge (up to n_e-1)
	std::size_t index = 0; // edge index in input order (up to M-1)
};

using BoostGraph = boost::adjacency_list< 
	boost::vecS, // edges of the graph are stored in a std::vector
	boost::vecS, // vertices of the graph are stored in a std::vector
	boost::undirectedS, // the graph is undirected (edges have no direction)
	VertexProperty,
	EdgeProperty
>;

using Vertex = boost::graph_traits<BoostGraph>::vertex_descriptor;
using Edge = boost::graph_traits<BoostGraph>::edge_descriptor;

using SparseMatrix = Eigen::SparseMatrix<double>;
using Vector = Eigen::VectorXd;

// Abstract base class defining the common interface.
class quantum_graph_problem {
public:
	virtual ~quantum_graph_problem() = default;

	// Default init: read input, build graph, set defaults.
	virtual void init(int argc, char *argv[]);

	// Evaluate known coefficients and store them in derived class.
	virtual void set_coefficients() = 0;

	// Assemble discrete matrices/vectors.
	virtual void assembly() = 0;

	// Solve the discrete problem.
	virtual bool solve() = 0;

	// Export solution for post-processing.
	virtual void sol_export() = 0;

protected:
	quantum_graph_problem() = default;

	// Read graph from a simple text file (list of edges).
	void read_graph_txt(const std::string &path);

	// Build global DoF numbering: internal nodes first, vertices last.
	void build_dof_numbering();

	// Return the global DoF associated with a local node on an edge.
	std::size_t edge_local_node_to_dof(const Edge &edge, int local_node) const;

	std::size_t n_dofs() const { return n_dofs_; }

	// Graph data.
	BoostGraph graph_;
	std::size_t n_vertices_orig_ = 0;
	std::size_t n_edges_ = 0;

	// DoF layout.
	std::size_t n_dofs_ = 0;
	std::size_t vertex_dof_offset_ = 0;
	std::vector<std::size_t> edge_dof_offset_;

	// Discrete problem data (common to all problems).
	SparseMatrix H;
	SparseMatrix M;
	Vector U;

	// Input path (graph file).
	std::string graph_file_;
};

} // end namespace femg

#endif
