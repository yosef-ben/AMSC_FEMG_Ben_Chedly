#include "fisher_kolmogorov_problem.hpp"

#include <algorithm>
#include <array>
#include <cctype>
#include <chrono>
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
		|| contains(name, "parstriangularis") || contains(name, "parsorbitalis")
		|| contains(name, "precentral")) {
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
		if (fields.size() != 8) {
			throw std::runtime_error("Unexpected node CSV row.");
		}
		NodeRecord node;
		node.local_id = static_cast<std::size_t>(std::stoul(fields[0]));
		node.parent_name = fields[2];
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
		if (fields.size() != 8) {
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
			{std::stod(fields[5]), std::stod(fields[6]), std::stod(fields[7])});
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
		if (fields.size() != 7) {
			throw std::runtime_error("Unexpected edge CSV row.");
		}
		const std::size_t edge_id =
			static_cast<std::size_t>(std::stoul(fields[0]));
		if (edge_id != weights.size()) {
			throw std::runtime_error("Edge IDs are not contiguous.");
		}
		weights.push_back(std::stod(fields[6]));
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
		// Twenty years is the horizon of Corti et al.; a longer one is
		// accepted so that the spreading can be followed to saturation.
		double final_time = 20.0;
		const double time_step = 0.2;
		const double background_concentration = 0.01;
		const double seed_concentration = 0.10;
		const std::size_t output_stride = 5;

		std::string graph_file = "data/connectome/fornari83/graph_fem_4.txt";
		if (argc >= 2) graph_file = argv[1];
		const bool performance_mode =
			argc >= 3 && std::string(argv[2]) == "--performance";
		char *graph_argv[] = {argv[0], graph_file.data()};
		// An explicit output directory keeps refined runs from overwriting the
		// four-element run that benchmark 21 reports.
		const std::string output_directory = performance_mode
			? "output/fisher_kolmogorov/corti83/performance_tmp"
			: (argc >= 4 ? std::string(argv[3])
			             : std::string("output/fisher_kolmogorov/corti83"));
		if (argc >= 5) final_time = std::stod(argv[4]);
		const bool uniform_rates =
			argc >= 6 && std::string(argv[5]) == "uniform";
		using clock = std::chrono::steady_clock;
		const auto total_start = clock::now();

		const auto nodes = read_nodes("data/connectome/fornari83/nodes.csv");
		const auto diffusion = read_normalized_edge_weights(
			"data/connectome/fornari83/edges.csv");
		const auto visualization_coordinates = read_visualization_coordinates(
			"data/connectome/fornari83/nodes.csv");
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
		if (uniform_rates) {
			// Control variant for the activation-order comparison: the seven
			// regional rates are replaced by their vertex mean, so the only
			// anatomy left in the model is the connectivity.
			const double mean_rate =
				std::accumulate(alpha.begin(), alpha.end(), 0.0)
				/ static_cast<double>(alpha.size());
			std::fill(alpha.begin(), alpha.end(), mean_rate);
		}

		femg::fisher_kolmogorov_problem problem(final_time, time_step);
		const auto init_start = clock::now();
		problem.init(2, graph_argv);
		const auto init_end = clock::now();
		if (problem.number_of_vertices() != nodes.size()
			|| problem.number_of_edges() != diffusion.size()) {
			throw std::runtime_error("Connectome metadata does not match FEM graph.");
		}
		const auto assembly_start = clock::now();
		problem.set_edge_diffusion_coefficients(diffusion);
		problem.set_vertex_reaction_coefficients(alpha);
		problem.set_vertex_initial_condition(initial_condition);
		problem.set_vertex_visualization_coordinates(visualization_coordinates);
		problem.set_time_scheme(
			femg::fisher_kolmogorov_problem::TimeScheme::corti_semi_implicit);
		problem.set_output_enabled(!performance_mode);
		problem.set_output_directory(output_directory);
		problem.set_output_stride(output_stride);
		problem.set_verbose(false);
		problem.set_coefficients();
		problem.assemble_matrices();
		const auto assembly_end = clock::now();
		if (!performance_mode) problem.print_matrix_summary(std::cout);

		std::filesystem::create_directories(output_directory);
		std::ofstream coefficients(
			output_directory + "/reaction_coefficients.csv");
		coefficients << "node_id,name,region,alpha\n";
		for (const NodeRecord &node : nodes) {
			const std::size_t region = static_cast<std::size_t>(node.region);
			coefficients << node.local_id << "," << node.parent_name << ","
				<< region_names[region] << "," << std::setprecision(16)
				<< alpha[node.local_id] << "\n";
		}
		std::ofstream averages(
			output_directory + "/regional_averages.csv");
		averages << "time,global,min,max";
		for (const std::string &name : region_names) {
			averages << "," << name;
		}
		averages << "\n";

		if (!performance_mode) problem.set_time_step_callback(
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

		if (!performance_mode) {
			std::cout << "  seed vertices: " << seed_count << "\n"
				<< "  diffusion scaling: D_e = w_e / max(w)\n"
				<< (uniform_rates
					? "  reaction coefficients: uniform vertex mean of the regional rates\n"
					: "  reaction coefficients: deterministic means from Corti et al.\n");
		}
		const auto solve_start = clock::now();
		problem.solve();
		const auto solve_end = clock::now();

		const double final_minimum = problem.solution().minCoeff();
		const double final_maximum = problem.solution().maxCoeff();
		if (final_minimum < -1.0e-8 || final_maximum > 1.0 + 1.0e-8) {
			throw std::runtime_error("Concentration left the physical interval [0, 1].");
		}
		if (performance_mode) {
			const auto seconds = [](auto begin, auto end) {
				return std::chrono::duration<double>(end - begin).count();
			};
			std::cout << std::setprecision(16)
				<< "PERFORMANCE," << problem.number_of_vertices()
				<< "," << problem.number_of_edges()
				<< "," << problem.number_of_dofs()
				<< "," << problem.mass_matrix().nonZeros()
				<< "," << seconds(init_start, init_end)
				<< "," << seconds(assembly_start, assembly_end)
				<< "," << seconds(solve_start, solve_end)
				<< "," << seconds(total_start, solve_end) << "\n";
		} else {
			std::cout << "Final concentration range: [" << final_minimum
				<< ", " << final_maximum << "]\n";
		}
		return 0;
	} catch (const std::exception &error) {
		std::cerr << "Error: " << error.what() << "\n";
		return 1;
	}
}
