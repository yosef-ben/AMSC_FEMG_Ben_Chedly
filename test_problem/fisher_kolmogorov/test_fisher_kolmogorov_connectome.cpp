#include "fisher_kolmogorov_problem.hpp"

#include <algorithm>
#include <array>
#include <cctype>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <numeric>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>

namespace {

enum class BrainRegion : std::size_t {
	frontal = 0,
	temporal,
	parietal,
	insular,
	limbic,
	occipital,
	subcortical,
	count
};

constexpr std::size_t n_regions =
	static_cast<std::size_t>(BrainRegion::count);

const std::array<std::string, n_regions> region_names = {
	"frontal", "temporal", "parietal", "insular", "limbic",
	"occipital", "subcortical"
};

// Mean values estimated in Table 1 of Corti et al.
const std::array<double, n_regions> mean_reaction_coefficients = {
	0.1801, 0.1421, 0.0627, 0.1005, 0.1351, 0.0545, 0.1147
};

struct NodeRecord {
	std::size_t local_id = 0;
	std::string parent_name;
	BrainRegion region = BrainRegion::subcortical;
};

std::vector<std::string> split_csv_line(const std::string &line) {
	std::vector<std::string> fields;
	std::stringstream stream(line);
	std::string field;
	while (std::getline(stream, field, ',')) {
		if (!field.empty() && field.back() == '\r') {
			field.pop_back();
		}
		fields.push_back(field);
	}
	return fields;
}

std::string lower_case(std::string text) {
	std::transform(text.begin(), text.end(), text.begin(),
		[](unsigned char character) {
			return static_cast<char>(std::tolower(character));
		});
	return text;
}

bool contains(const std::string &text, const std::string &token) {
	return text.find(token) != std::string::npos;
}

BrainRegion classify_region(const std::string &parent_name) {
	const std::string name = lower_case(parent_name);
	if (contains(name, "insula")) {
		return BrainRegion::insular;
	}
	if (contains(name, "cingulate") || contains(name, "entorhinal")
		|| contains(name, "parahippocampal") || contains(name, "hippocampus")
		|| contains(name, "amygdala")) {
		return BrainRegion::limbic;
	}
	if (contains(name, "cuneus") || contains(name, "occipital")
		|| contains(name, "lingual") || contains(name, "pericalcarine")) {
		return BrainRegion::occipital;
	}
	if (contains(name, "parietal") || contains(name, "supramarginal")
		|| contains(name, "precuneus") || contains(name, "postcentral")
		|| contains(name, "paracentral")) {
		return BrainRegion::parietal;
	}
	if (contains(name, "temporal") || contains(name, "bankssts")
		|| contains(name, "fusiform")) {
		return BrainRegion::temporal;
	}
	if (contains(name, "frontal") || contains(name, "orbitofrontal")
		|| contains(name, "parsopercularis")
		|| contains(name, "parstriangularis") || contains(name, "precentral")) {
		return BrainRegion::frontal;
	}
	if (contains(name, "caudate") || contains(name, "putamen")
		|| contains(name, "pallidum") || contains(name, "thalamus")
		|| contains(name, "accumbens") || contains(name, "brain-stem")) {
		return BrainRegion::subcortical;
	}
	throw std::runtime_error("Unable to classify anatomical region: " + parent_name);
}

std::vector<NodeRecord> read_nodes(const std::string &path) {
	std::ifstream input(path);
	if (!input) {
		throw std::runtime_error("Unable to open node metadata: " + path);
	}
	std::string line;
	std::getline(input, line);
	std::vector<NodeRecord> nodes;
	while (std::getline(input, line)) {
		const auto fields = split_csv_line(line);
		if (fields.size() != 10) {
			throw std::runtime_error("Unexpected node CSV row.");
		}
		NodeRecord node;
		node.local_id = static_cast<std::size_t>(std::stoul(fields[0]));
		node.parent_name = fields[7];
		node.region = classify_region(node.parent_name);
		nodes.push_back(std::move(node));
	}
	std::sort(nodes.begin(), nodes.end(), [](const NodeRecord &a, const NodeRecord &b) {
		return a.local_id < b.local_id;
	});
	for (std::size_t index = 0; index < nodes.size(); ++index) {
		if (nodes[index].local_id != index) {
			throw std::runtime_error("Node IDs are not contiguous.");
		}
	}
	return nodes;
}

std::vector<std::array<double, 3>> read_visualization_coordinates(
	const std::string &path) {
	std::ifstream input(path);
	if (!input) {
		throw std::runtime_error(
			"Unable to open anatomical coordinates: " + path);
	}
	std::string line;
	std::getline(input, line);
	std::vector<std::array<double, 3>> coordinates;
	while (std::getline(input, line)) {
		const auto fields = split_csv_line(line);
		if (fields.size() != 5) {
			throw std::runtime_error(
				"Unexpected anatomical coordinate CSV row.");
		}
		const std::size_t local_id =
			static_cast<std::size_t>(std::stoul(fields[0]));
		if (local_id != coordinates.size()) {
			throw std::runtime_error(
				"Anatomical coordinate IDs are not contiguous.");
		}
		coordinates.push_back(
			{std::stod(fields[2]), std::stod(fields[3]), std::stod(fields[4])});
	}
	return coordinates;
}

std::vector<double> read_normalized_edge_weights(const std::string &path) {
	std::ifstream input(path);
	if (!input) {
		throw std::runtime_error("Unable to open edge metadata: " + path);
	}
	std::string line;
	std::getline(input, line);
	std::vector<double> weights;
	while (std::getline(input, line)) {
		const auto fields = split_csv_line(line);
		if (fields.size() != 8) {
			throw std::runtime_error("Unexpected edge CSV row.");
		}
		const std::size_t edge_id =
			static_cast<std::size_t>(std::stoul(fields[0]));
		if (edge_id != weights.size()) {
			throw std::runtime_error("Edge IDs are not contiguous.");
		}
		weights.push_back(std::stod(fields[5]));
	}
	const double maximum = *std::max_element(weights.begin(), weights.end());
	for (double &weight : weights) {
		weight /= maximum;
	}
	return weights;
}

bool is_seed_region(const std::string &parent_name) {
	const std::string name = lower_case(parent_name);
	return contains(name, "entorhinal") || contains(name, "hippocampus");
}

} // namespace

int main(int argc, char *argv[]) {
	try {
		const double final_time = 20.0;
		const double time_step = 0.2;
		const double background_concentration = 0.01;
		const double seed_concentration = 0.10;
		const std::size_t output_stride = 5;

		std::string graph_file = "data/connectome/budapest_lcc_fem.txt";
		char *local_argv[] = {argv[0], graph_file.data()};
		const int graph_argc = (argc >= 2) ? argc : 2;
		char **graph_argv = (argc >= 2) ? argv : local_argv;

		const auto nodes = read_nodes("data/connectome/budapest_lcc_nodes.csv");
		const auto diffusion = read_normalized_edge_weights(
			"data/connectome/budapest_lcc_edges.csv");
		const auto visualization_coordinates = read_visualization_coordinates(
			"data/connectome/budapest_lcc_anatomical_coordinates.csv");
		std::vector<double> alpha(nodes.size());
		std::vector<double> initial_condition(nodes.size(), background_concentration);
		std::array<std::vector<std::size_t>, n_regions> vertices_by_region;
		std::size_t seed_count = 0;
		for (const NodeRecord &node : nodes) {
			const auto region = static_cast<std::size_t>(node.region);
			alpha[node.local_id] = mean_reaction_coefficients[region];
			vertices_by_region[region].push_back(node.local_id);
			if (is_seed_region(node.parent_name)) {
				initial_condition[node.local_id] = seed_concentration;
				++seed_count;
			}
		}

		femg::fisher_kolmogorov_problem problem(final_time, time_step);
		problem.init(graph_argc, graph_argv);
		if (problem.number_of_vertices() != nodes.size()
			|| problem.number_of_edges() != diffusion.size()) {
			throw std::runtime_error("Connectome metadata does not match FEM graph.");
		}
		problem.set_edge_diffusion_coefficients(diffusion);
		problem.set_vertex_reaction_coefficients(alpha);
		problem.set_vertex_initial_condition(initial_condition);
		problem.set_vertex_visualization_coordinates(visualization_coordinates);
		problem.set_time_scheme(
			femg::fisher_kolmogorov_problem::TimeScheme::corti_semi_implicit);
		problem.set_output_directory("output/fisher_kolmogorov/connectome");
		problem.set_output_stride(output_stride);
		problem.set_verbose(false);
		problem.set_coefficients();
		problem.assemble_matrices();
		problem.print_matrix_summary(std::cout);

		std::filesystem::create_directories("output/fisher_kolmogorov/connectome");
		std::ofstream averages(
			"output/fisher_kolmogorov/connectome/regional_averages.csv");
		averages << "time,global,min,max";
		for (const std::string &name : region_names) {
			averages << "," << name;
		}
		averages << "\n";

		problem.set_time_step_callback(
			[&](std::size_t, double time, const femg::Vector &) {
				double global_sum = 0.0;
				double minimum = 1.0;
				double maximum = 0.0;
				std::array<double, n_regions> region_sum{};
				for (std::size_t vertex = 0; vertex < nodes.size(); ++vertex) {
					const double value = problem.vertex_value(vertex);
					global_sum += value;
					minimum = std::min(minimum, value);
					maximum = std::max(maximum, value);
					region_sum[static_cast<std::size_t>(nodes[vertex].region)] += value;
				}
				averages << std::setprecision(16) << time << ","
					<< global_sum / nodes.size() << "," << minimum << "," << maximum;
				for (std::size_t region = 0; region < n_regions; ++region) {
					averages << "," << region_sum[region]
						/ vertices_by_region[region].size();
				}
				averages << "\n";
			});

		std::cout << "  seed vertices: " << seed_count << "\n"
			<< "  diffusion scaling: D_e = w_e / max(w)\n"
			<< "  reaction coefficients: deterministic means from Corti et al.\n";
		problem.solve();

		const double final_minimum = problem.solution().minCoeff();
		const double final_maximum = problem.solution().maxCoeff();
		if (final_minimum < -1.0e-8 || final_maximum > 1.0 + 1.0e-8) {
			throw std::runtime_error("Concentration left the physical interval [0, 1].");
		}
		std::cout << "Final concentration range: [" << final_minimum
			<< ", " << final_maximum << "]\n";
		return 0;
	} catch (const std::exception &error) {
		std::cerr << "Error: " << error.what() << "\n";
		return 1;
	}
}
