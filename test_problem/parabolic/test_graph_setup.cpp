#include "parabolic_problem.hpp"

#include <algorithm>
#include <iostream>
#include <string>
#include <vector>

class graph_setup_debug_problem : public femg::parabolic_problem {
public:
	graph_setup_debug_problem()
		: femg::parabolic_problem(1.0, 1.0, 0.05, 1.0) {}

	void print_summary() const {
		std::cout << "Graph setup summary\n";
		std::cout << "  vertices: " << n_vertices_orig_ << "\n";
		std::cout << "  edges:    " << n_edges_ << "\n";
		std::cout << "  dofs:     " << n_dofs() << "\n\n";

		std::vector<femg::Edge> edges;
		auto [ei, ei_end] = boost::edges(graph_);
		for (auto it = ei; it != ei_end; ++it) {
			edges.push_back(*it);
		}

		std::sort(edges.begin(), edges.end(),
							[this](const femg::Edge &a, const femg::Edge &b) {
								return graph_[a].index < graph_[b].index;
							});

		for (const auto &edge : edges) {
			const auto source = boost::source(edge, graph_);
			const auto target = boost::target(edge, graph_);
			const auto &data = graph_[edge];

			std::cout << "Edge " << data.index << ": "
								<< graph_[source].index << " -- "
								<< graph_[target].index
								<< ", length = " << data.length
								<< ", cells = " << data.n_cells << "\n";

			std::cout << "  local node -> global dof\n";
			for (int local_node = 0; local_node <= data.n_cells; ++local_node) {
				const std::size_t global_dof =
					edge_local_node_to_dof(edge, local_node);
				std::cout << "    " << local_node << " -> " << global_dof << "\n";
			}
		}
	}
};

int main(int argc, char *argv[]) {
	std::string graph_file = "data/interval_1d.txt";
	char *local_argv[] = {argv[0], graph_file.data()};

	graph_setup_debug_problem problem;
	if (argc >= 2) {
		problem.init(argc, argv);
	} else {
		problem.init(2, local_argv);
	}

	problem.print_summary();
	return 0;
}
