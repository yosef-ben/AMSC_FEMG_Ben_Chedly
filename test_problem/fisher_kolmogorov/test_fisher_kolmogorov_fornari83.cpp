#include "fisher_kolmogorov_problem.hpp"

#include <Eigen/SparseLU>

#include <algorithm>
#include <array>
#include <cctype>
#include <cmath>
#include <filesystem>
#include <fstream>
#include <functional>
#include <iomanip>
#include <iostream>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>

namespace {

struct Node {
	std::size_t id;
	std::string name;
	int lobe;
	std::array<double, 3> coordinates;
};

struct Edge {
	std::size_t source;
	std::size_t target;
	double weight;
};

std::vector<std::string> split(const std::string &line) {
	std::vector<std::string> fields;
	std::stringstream stream(line);
	std::string field;
	while (std::getline(stream, field, ',')) {
		if (!field.empty() && field.back() == '\r') field.pop_back();
		fields.push_back(field);
	}
	return fields;
}

std::string lower(std::string value) {
	std::transform(value.begin(), value.end(), value.begin(),
		[](unsigned char c) { return static_cast<char>(std::tolower(c)); });
	return value;
}

bool has(const std::string &text, const std::string &token) {
	return text.find(token) != std::string::npos;
}

// 0: temporal, 1: frontal, 2: parietal, 3: occipital, -1: other.
int classify(const std::string &raw) {
	const std::string name = lower(raw);
	if (has(name, "temporal") || has(name, "bankssts")
		|| has(name, "entorhinal") || has(name, "fusiform")
		|| has(name, "parahippocampal")) return 0;
	if (has(name, "frontal") || has(name, "orbitofrontal")
		|| has(name, "parsopercularis") || has(name, "parsorbitalis")
		|| has(name, "parstriangularis") || has(name, "precentral")) return 1;
	if (has(name, "parietal") || has(name, "postcentral")
		|| has(name, "precuneus") || has(name, "supramarginal")
		|| has(name, "paracentral")) return 2;
	if (has(name, "cuneus") || has(name, "occipital")
		|| has(name, "lingual") || has(name, "pericalcarine")) return 3;
	return -1;
}

std::vector<Node> read_nodes(const std::string &path) {
	std::ifstream input(path);
	if (!input) throw std::runtime_error("Unable to open " + path);
	std::string line;
	std::getline(input, line);
	std::vector<Node> nodes;
	while (std::getline(input, line)) {
		const auto fields = split(line);
		if (fields.size() != 8) throw std::runtime_error("Invalid node CSV row.");
		Node node{
			static_cast<std::size_t>(std::stoul(fields[0])),
			fields[2],
			classify(fields[2]),
			{std::stod(fields[5]), std::stod(fields[6]), std::stod(fields[7])}
		};
		if (node.id != nodes.size()) throw std::runtime_error("Non-contiguous node IDs.");
		nodes.push_back(std::move(node));
	}
	return nodes;
}

std::vector<Edge> read_edges(const std::string &path) {
	std::ifstream input(path);
	if (!input) throw std::runtime_error("Unable to open " + path);
	std::string line;
	std::getline(input, line);
	std::vector<Edge> edges;
	while (std::getline(input, line)) {
		const auto fields = split(line);
		if (fields.size() != 7) throw std::runtime_error("Invalid edge CSV row.");
		if (static_cast<std::size_t>(std::stoul(fields[0])) != edges.size())
			throw std::runtime_error("Non-contiguous edge IDs.");
		edges.push_back({
			static_cast<std::size_t>(std::stoul(fields[1])),
			static_cast<std::size_t>(std::stoul(fields[2])),
			std::stod(fields[6])
		});
	}
	return edges;
}

femg::SparseMatrix make_laplacian(
	std::size_t n_nodes, const std::vector<Edge> &edges) {
	std::vector<Eigen::Triplet<double>> triplets;
	triplets.reserve(4 * edges.size());
	for (const Edge &edge : edges) {
		const auto i = static_cast<Eigen::Index>(edge.source);
		const auto j = static_cast<Eigen::Index>(edge.target);
		triplets.emplace_back(i, i, edge.weight);
		triplets.emplace_back(j, j, edge.weight);
		triplets.emplace_back(i, j, -edge.weight);
		triplets.emplace_back(j, i, -edge.weight);
	}
	femg::SparseMatrix matrix(
		static_cast<Eigen::Index>(n_nodes), static_cast<Eigen::Index>(n_nodes));
	matrix.setFromTriplets(triplets.begin(), triplets.end());
	matrix.makeCompressed();
	return matrix;
}

femg::Vector nodal_step(
	const femg::Vector &old,
	const femg::SparseMatrix &laplacian,
	double alpha,
	double dt) {
	femg::Vector solution = old;
	for (std::size_t iteration = 0; iteration < 30; ++iteration) {
		const femg::Vector residual = solution - old + dt * (laplacian * solution)
			- dt * alpha
				* (solution.array() * (1.0 - solution.array())).matrix();
		if (residual.norm() < 1.0e-12) return solution;

		femg::SparseMatrix jacobian = dt * laplacian;
		for (Eigen::Index i = 0; i < solution.size(); ++i)
			jacobian.coeffRef(i, i) +=
				1.0 - dt * alpha * (1.0 - 2.0 * solution(i));
		jacobian.makeCompressed();
		Eigen::SparseLU<femg::SparseMatrix> solver;
		solver.compute(jacobian);
		if (solver.info() != Eigen::Success)
			throw std::runtime_error("Nodal Newton factorization failed.");
		const femg::Vector update = solver.solve(-residual);
		if (solver.info() != Eigen::Success)
			throw std::runtime_error("Nodal Newton solve failed.");
		solution += update;
		if (update.norm() < 1.0e-12 * (1.0 + solution.norm())) return solution;
	}
	throw std::runtime_error("Nodal Newton method did not converge.");
}

void write_header(std::ofstream &output) {
	output << "time,global,temporal,frontal,parietal,occipital,min,max\n";
}

void write_profile_header(std::ofstream &output, std::size_t n_nodes) {
	output << "time";
	for (std::size_t node = 0; node < n_nodes; ++node)
		output << ",node_" << node;
	output << "\n";
}

void write_profile(
	std::ofstream &output,
	double time,
	std::size_t n_nodes,
	const std::function<double(std::size_t)> &value) {
	output << std::setprecision(16) << time;
	for (std::size_t node = 0; node < n_nodes; ++node)
		output << "," << value(node);
	output << "\n";
}

void write_biomarkers(
	std::ofstream &output,
	double time,
	const std::vector<Node> &nodes,
	const std::function<double(std::size_t)> &value) {
	std::array<double, 4> sums{};
	std::array<std::size_t, 4> counts{};
	double global = 0.0;
	double minimum = 1.0;
	double maximum = 0.0;
	for (const Node &node : nodes) {
		const double concentration = value(node.id);
		global += concentration;
		minimum = std::min(minimum, concentration);
		maximum = std::max(maximum, concentration);
		if (node.lobe >= 0) {
			sums[static_cast<std::size_t>(node.lobe)] += concentration;
			++counts[static_cast<std::size_t>(node.lobe)];
		}
	}
	output << std::setprecision(16) << time << ","
		<< 100.0 * global / nodes.size();
	for (std::size_t lobe = 0; lobe < sums.size(); ++lobe)
		output << "," << 100.0 * sums[lobe] / counts[lobe];
	output << "," << minimum << "," << maximum << "\n";
}

} // namespace

int main(int argc, char *argv[]) {
	try {
		const std::string data_dir = "data/connectome/fornari83";
		const std::size_t cells_per_edge =
			(argc >= 2) ? static_cast<std::size_t>(std::stoul(argv[1])) : 1;
		const double dt = (argc >= 3) ? std::stod(argv[2]) : 0.4;
		const std::string output_dir = (argc >= 4)
			? argv[3]
			: "output/fisher_kolmogorov/fornari83";
		const bool write_vtk = (argc < 4);
		const double alpha = (argc >= 5) ? std::stod(argv[4]) : 0.5;
		const double final_time = (argc >= 6) ? std::stod(argv[5]) : 40.0;
		// Uniform scaling of the connectivity weights. rho = 1 is the literal
		// Laplacian of Fornari et al.; benchmark 21 uses rho = 1/max(w).
		const double diffusion_scaling = (argc >= 7) ? std::stod(argv[6]) : 1.0;
		if (diffusion_scaling < 0.0)
			throw std::invalid_argument(
				"The diffusion scaling must be nonnegative.");
		// "be":    Backward Euler with Newton, as in Fornari et al.
		// "cn":    the semi-implicit Crank-Nicolson scheme of Corti et al.
		// "nodal": nodal reference only, without the metric-graph FEM.
		// A "_lumped" suffix runs the FEM with the row-sum lumped mass and
		// the vertex-rule reaction, the stabilized variant of benchmark 23.
		std::string scheme = (argc >= 8) ? argv[7] : "be";
		bool lumped = false;
		if (scheme.size() > 7 && scheme.substr(scheme.size() - 7) == "_lumped") {
			lumped = true;
			scheme = scheme.substr(0, scheme.size() - 7);
		}
		if (scheme != "be" && scheme != "cn" && scheme != "nodal")
			throw std::invalid_argument(
				"The time scheme must be be, cn, nodal, be_lumped or cn_lumped.");
		// Seeding: absent or negative keeps the entorhinal seeding of the
		// paper (the tau pattern); a vertex index seeds that vertex alone;
		// the keyword "neocortex" seeds every vertex of the four cortical
		// lobes, which is where amyloid-beta deposits first appear
		// (Weickenmeier et al., section 2.7).
		const std::string seed_argument = (argc >= 9) ? argv[8] : "-1";
		const bool seed_neocortex = seed_argument == "neocortex";
		const int seed_vertex =
			seed_neocortex ? -1 : std::stoi(seed_argument);
		const std::size_t n_steps =
			static_cast<std::size_t>(std::llround(final_time / dt));
		if (std::abs(n_steps * dt - final_time) > 1.0e-10)
			throw std::invalid_argument(
				"The final time must be a multiple of dt.");

		const auto nodes = read_nodes(data_dir + "/nodes.csv");
		const auto edges = read_edges(data_dir + "/edges.csv");
		if (nodes.size() != 83 || edges.size() != 1130)
			throw std::runtime_error("Unexpected FreeSurfer graph dimensions.");
		std::filesystem::create_directories(output_dir);

		femg::Vector initial =
			femg::Vector::Zero(static_cast<Eigen::Index>(nodes.size()));
		std::size_t seeds = 0;
		if (seed_neocortex) {
			// Amyloid-beta seeding: every vertex of the four cortical lobes,
			// at the same level the entorhinal seeding uses.
			for (const Node &node : nodes) {
				if (node.lobe >= 0) {
					initial(static_cast<Eigen::Index>(node.id)) = 0.1;
					++seeds;
				}
			}
			if (seeds != 58)
				throw std::runtime_error("Expected 58 neocortical seeds.");
		} else if (seed_vertex < 0) {
			// Default seeding of Fornari et al.: the entorhinal cortex.
			for (const Node &node : nodes) {
				if (has(lower(node.name), "entorhinal")) {
					initial(static_cast<Eigen::Index>(node.id)) = 0.1;
					++seeds;
				}
			}
			if (seeds != 2)
				throw std::runtime_error("Expected two entorhinal seeds.");
		} else {
			// Single-region seeding, for the regional vulnerability study.
			if (static_cast<std::size_t>(seed_vertex) >= nodes.size())
				throw std::invalid_argument("Seed vertex out of range.");
			initial(seed_vertex) = 0.1;
			seeds = 1;
		}

		const femg::SparseMatrix laplacian =
			diffusion_scaling * make_laplacian(nodes.size(), edges);
		femg::Vector nodal = initial;
		std::ofstream nodal_output(output_dir + "/nodal_biomarkers.csv");
		std::ofstream nodal_profiles(output_dir + "/nodal_profiles.csv");
		write_header(nodal_output);
		write_profile_header(nodal_profiles, nodes.size());
		write_biomarkers(nodal_output, 0.0, nodes, [&](std::size_t i) {
			return nodal(static_cast<Eigen::Index>(i));
		});
		write_profile(nodal_profiles, 0.0, nodes.size(), [&](std::size_t i) {
			return nodal(static_cast<Eigen::Index>(i));
		});
		for (std::size_t step = 1; step <= n_steps; ++step) {
			nodal = nodal_step(nodal, laplacian, alpha, dt);
			write_biomarkers(nodal_output, step * dt, nodes, [&](std::size_t i) {
				return nodal(static_cast<Eigen::Index>(i));
			});
			write_profile(nodal_profiles, step * dt, nodes.size(),
				[&](std::size_t i) {
					return nodal(static_cast<Eigen::Index>(i));
				});
		}

		if (scheme == "nodal") {
			std::cout << "Fornari 83-region nodal reference only\n"
				<< "  alpha: " << alpha << "\n"
				<< "  diffusion scaling: " << diffusion_scaling << "\n"
				<< "  dt: " << dt << ", T: " << final_time << "\n"
				<< "  nodal final range: [" << nodal.minCoeff()
				<< ", " << nodal.maxCoeff() << "]\n";
			return 0;
		}

		std::string graph_file = data_dir + "/graph_fem_"
			+ std::to_string(cells_per_edge) + ".txt";
		char program[] = "test_fisher_kolmogorov_fornari83";
		char *arguments[] = {program, graph_file.data()};
		femg::fisher_kolmogorov_problem fem(final_time, dt);
		fem.init(2, arguments);

		std::vector<double> diffusion;
		std::vector<double> reaction(nodes.size(), alpha);
		std::vector<double> initial_values(nodes.size());
		std::vector<std::array<double, 3>> coordinates;
		diffusion.reserve(edges.size());
		coordinates.reserve(nodes.size());
		for (const Edge &edge : edges)
			diffusion.push_back(diffusion_scaling * edge.weight);
		for (const Node &node : nodes) {
			initial_values[node.id] = initial(static_cast<Eigen::Index>(node.id));
			coordinates.push_back(node.coordinates);
		}

		fem.set_edge_diffusion_coefficients(diffusion);
		fem.set_vertex_reaction_coefficients(reaction);
		fem.set_vertex_initial_condition(initial_values);
		fem.set_vertex_visualization_coordinates(coordinates);
		fem.set_time_scheme(scheme == "cn"
			? femg::fisher_kolmogorov_problem::TimeScheme::corti_semi_implicit
			: femg::fisher_kolmogorov_problem::TimeScheme::backward_euler);
		fem.set_mass_lumping(lumped);
		// Weakly diffusive graphs need a damped Newton solve with more than
		// the default number of iterations; see benchmark 23.
		fem.set_newton_parameters(1.0e-11, 200);
		fem.set_output_directory(output_dir + "/fem_solution");
		fem.set_output_stride(10);
		fem.set_output_enabled(write_vtk);
		fem.set_verbose(false);
		fem.set_coefficients();
		fem.assemble_matrices();

		// At one element per connection and unit lengths the assembled
		// diffusion matrix must be the graph Laplacian the nodal reference
		// uses, entry by entry; the check ties the two models to the same
		// operator and is reported with the run.
		double laplacian_gap = -1.0;
		if (cells_per_edge == 1) {
			const femg::SparseMatrix difference =
				fem.diffusion_matrix() - laplacian;
			laplacian_gap = 0.0;
			for (int k = 0; k < difference.outerSize(); ++k)
				for (femg::SparseMatrix::InnerIterator it(difference, k); it;
					++it)
					laplacian_gap = std::max(laplacian_gap, std::abs(it.value()));
			std::cout << "  max |K_D - L| at one element per connection: "
				<< std::setprecision(3) << laplacian_gap << "\n";
		}

		std::ofstream fem_output(output_dir + "/fem_biomarkers.csv");
		std::ofstream fem_profiles(output_dir + "/fem_profiles.csv");
		std::ofstream fem_metric_mass(output_dir + "/fem_metric_mass.csv");
		write_header(fem_output);
		write_profile_header(fem_profiles, nodes.size());
		fem_metric_mass << "time,metric_global\n";
		const femg::Vector mass_weights = fem.mass_matrix()
			* femg::Vector::Ones(
				static_cast<Eigen::Index>(fem.number_of_dofs()));
		const double domain_measure = mass_weights.sum();
		fem.set_time_step_callback(
			[&](std::size_t, double time, const femg::Vector &state) {
				write_biomarkers(fem_output, time, nodes, [&](std::size_t i) {
					return fem.vertex_value(i);
				});
				write_profile(fem_profiles, time, nodes.size(),
					[&](std::size_t i) { return fem.vertex_value(i); });
				fem_metric_mass << std::setprecision(16) << time << ","
					<< 100.0 * mass_weights.dot(state) / domain_measure
					<< "\n";
			});
		fem.solve();

		if (nodal.minCoeff() < -1.0e-8 || nodal.maxCoeff() > 1.0 + 1.0e-8
			|| fem.solution().minCoeff() < -1.0e-8
			|| fem.solution().maxCoeff() > 1.0 + 1.0e-8)
		{
			std::ostringstream message;
			message << "A concentration left [0,1]: nodal ["
				<< nodal.minCoeff() << ", " << nodal.maxCoeff() << "], FEM ["
				<< fem.solution().minCoeff() << ", "
				<< fem.solution().maxCoeff() << "].";
			throw std::runtime_error(message.str());
		}

		std::cout << "Fornari 83-region Fisher-Kolmogorov comparison\n"
			<< "  vertices: " << nodes.size() << "\n"
			<< "  edges: " << edges.size() << "\n"
			<< "  entorhinal seeds: " << seeds << "\n"
			<< "  alpha: " << alpha << "\n"
			<< "  diffusion scaling: " << diffusion_scaling << "\n"
			<< "  mass matrix: " << (lumped ? "lumped" : "consistent") << "\n"
			<< "  cells per edge: " << cells_per_edge
			<< ", DoFs: " << fem.number_of_dofs() << "\n"
			<< "  dt: " << dt << ", T: " << final_time << "\n"
			<< "  nodal final range: [" << nodal.minCoeff()
			<< ", " << nodal.maxCoeff() << "]\n"
			<< "  FEM final range: [" << fem.solution().minCoeff()
			<< ", " << fem.solution().maxCoeff() << "]\n";
		return 0;
	} catch (const std::exception &error) {
		std::cerr << "Error: " << error.what() << "\n";
		return 1;
	}
}
