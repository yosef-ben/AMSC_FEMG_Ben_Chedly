#include "quantum_graph_problem.hpp"

#include <algorithm>
#include <fstream>
#include <stdexcept>

namespace femg {

void quantum_graph_problem::init(int argc, char *argv[]) {
	if (argc < 2 || argv == nullptr) {
		throw std::invalid_argument("Missing graph input file path.");
	}

	graph_file_ = argv[1];
	read_graph_txt(graph_file_);
	build_dof_numbering();
}

void quantum_graph_problem::read_graph_txt(const std::string &path) {
	std::ifstream in(path);
	if (!in) {
		throw std::runtime_error("Unable to open graph file: " + path);
	}

	std::size_t n_vertices = 0;
	std::size_t n_edges = 0;
	if (!(in >> n_vertices >> n_edges)) {
		throw std::runtime_error("Invalid header in graph file: " + path);
	}

	n_vertices_orig_ = n_vertices;
	n_edges_ = n_edges;

	graph_.clear();
	edge_dof_offset_.clear();

	for (std::size_t i = 0; i < n_vertices_orig_; ++i) {
		Vertex v = boost::add_vertex(graph_);
		graph_[v].index = i;
		graph_[v].coords = {static_cast<double>(i), 0.0};
	}

	std::vector<double> data;
	double value = 0.0;
	while (in >> value) {
		data.push_back(value);
	}

	const std::size_t old_format_size = 4 * n_edges_;
	const std::size_t coordinate_format_size =
		2 * n_vertices_orig_ + 4 * n_edges_;

	std::size_t offset = 0;
	if (data.size() == coordinate_format_size) {
		for (std::size_t i = 0; i < n_vertices_orig_; ++i) {
			graph_[i].coords = {data[offset], data[offset + 1]};
			offset += 2;
		}
	} else if (data.size() != old_format_size) {
		throw std::runtime_error("Invalid graph data size in file: " + path);
	}

	for (std::size_t e = 0; e < n_edges_; ++e) {
		const std::size_t u = static_cast<std::size_t>(data[offset]);
		const std::size_t v = static_cast<std::size_t>(data[offset + 1]);
		const double length = data[offset + 2];
		const int n_cells = static_cast<int>(data[offset + 3]);
		offset += 4;

		if (u >= n_vertices_orig_ || v >= n_vertices_orig_) {
			throw std::runtime_error("Edge references invalid vertex index.");
		}

		if (length <= 0.0 || n_cells < 1) {
			throw std::runtime_error("Edge length or n_cells invalid.");
		}

		auto [edge, ok] = boost::add_edge(u, v, graph_);
		if (!ok) {
			throw std::runtime_error("Failed to add edge to graph.");
		}

		graph_[edge].length = length;
		graph_[edge].n_cells = n_cells;
		graph_[edge].index = e;
	}
}

void quantum_graph_problem::build_dof_numbering() {
	edge_dof_offset_.assign(n_edges_, 0);

	std::vector<Edge> edges;
	edges.reserve(n_edges_);
	auto [ei, ei_end] = boost::edges(graph_);
	for (auto it = ei; it != ei_end; ++it) {
		edges.push_back(*it);
	}

	std::sort(edges.begin(), edges.end(),
						[this](const Edge &a, const Edge &b) {
							return graph_[a].index < graph_[b].index;
						});

	std::size_t offset = 0;
	for (const auto &edge : edges) {
		const std::size_t e_idx = graph_[edge].index;
		edge_dof_offset_[e_idx] = offset;
		const int n_cells = graph_[edge].n_cells;
		offset += static_cast<std::size_t>(n_cells - 1);
	}

	vertex_dof_offset_ = offset;
	n_dofs_ = vertex_dof_offset_ + n_vertices_orig_;
}

std::size_t quantum_graph_problem::edge_local_node_to_dof(
	const Edge &edge,
	int local_node) const {
	const int n_cells = graph_[edge].n_cells;

	if (local_node < 0 || local_node > n_cells) {
		throw std::out_of_range("Local node index is outside the edge.");
	}

	if (local_node == 0) {
		const Vertex source = boost::source(edge, graph_);
		return vertex_dof_offset_ + graph_[source].index;
	}

	if (local_node == n_cells) {
		const Vertex target = boost::target(edge, graph_);
		return vertex_dof_offset_ + graph_[target].index;
	}

	const std::size_t edge_index = graph_[edge].index;
	return edge_dof_offset_[edge_index]
		+ static_cast<std::size_t>(local_node - 1);
}

} // namespace femg
