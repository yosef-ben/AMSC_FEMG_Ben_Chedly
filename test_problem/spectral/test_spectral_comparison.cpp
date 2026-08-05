#include "spectral_problem.hpp"

#include <array>
#include <cmath>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <stdexcept>
#include <string>
#include <vector>

namespace {

constexpr double pi = 3.141592653589793238462643383279502884;

struct CaseData {
	std::string name;
	std::string graph_file;
	std::size_t n_modes;
};

enum class AngleLaw {
	Fixed,
	InverseLevel
};

enum class LengthLaw {
	Fixed,
	InverseLevel,
	InverseLevelSquared
};

struct TreeVariant {
	std::string name;
	AngleLaw angle_law;
	double fixed_angle;
	LengthLaw length_law;
};

double level_length(int level, LengthLaw law) {
	if (level <= 0) {
		return 1.0;
	}

	const int bifurcation_level = std::max(level - 1, 1);

	switch (law) {
	case LengthLaw::Fixed:
		return 1.0;
	case LengthLaw::InverseLevel:
		return 1.0 / static_cast<double>(bifurcation_level);
	case LengthLaw::InverseLevelSquared:
		return 1.0 / static_cast<double>(
			bifurcation_level * bifurcation_level);
	}

	return 1.0;
}

double level_angle(int level, AngleLaw law, double fixed_angle) {
	if (law == AngleLaw::Fixed) {
		return fixed_angle;
	}

	const int safe_level = std::max(level, 1);
	return (pi / 4.0) / static_cast<double>(safe_level);
}

int cells_for_length(double length) {
	(void) length;
	return 100;
}

std::string write_tree_variant(const TreeVariant &variant) {
	std::filesystem::create_directories("output/spectral/comparison/graphs");

	const std::string filename =
		"output/spectral/comparison/graphs/" + variant.name + ".txt";

	const std::array<std::array<int, 2>, 15> edges = {{
		{{0, 1}}, {{1, 2}}, {{1, 3}}, {{2, 4}}, {{2, 5}},
		{{3, 6}}, {{3, 7}}, {{4, 8}}, {{4, 9}}, {{5, 10}},
		{{5, 11}}, {{6, 12}}, {{6, 13}}, {{7, 14}}, {{7, 15}}
	}};

	const std::array<int, 16> levels = {{
		0, 1, 2, 2, 3, 3, 3, 3,
		4, 4, 4, 4, 4, 4, 4, 4
	}};

	std::array<std::array<double, 2>, 16> coords{};
	std::array<double, 16> directions{};
	coords[0] = {{0.0, 0.0}};
	directions[0] = pi / 2.0;

	for (const auto &edge : edges) {
		const int parent = edge[0];
		const int child = edge[1];
		const int level = levels[child];
		const double length = level_length(level, variant.length_law);

		double direction = pi / 2.0;
		if (parent == 0) {
			direction = pi / 2.0;
		} else {
			const double angle = level_angle(level, variant.angle_law, variant.fixed_angle);
			const double sign = (child % 2 == 0) ? -1.0 : 1.0;
			direction = directions[parent] + sign * angle;
		}

		directions[child] = direction;
		coords[child] = {{
			coords[parent][0] + length * std::cos(direction),
			coords[parent][1] + length * std::sin(direction)
		}};
	}

	std::ofstream out(filename);
	if (!out) {
		throw std::runtime_error("Unable to write tree graph: " + filename);
	}

	out << "16 15\n";
	for (const auto &coord : coords) {
		out << std::setprecision(16) << coord[0] << " " << coord[1] << "\n";
	}

	for (const auto &edge : edges) {
		const int child = edge[1];
		const int level = levels[child];
		const double length = level_length(level, variant.length_law);
		out << edge[0] << " " << edge[1] << " "
				<< std::setprecision(16) << length << " "
				<< cells_for_length(length) << "\n";
	}

	return filename;
}

void write_comparison(const CaseData &case_data) {
	std::string executable_name = "test_spectral_comparison";
	std::string graph_file = case_data.graph_file;
	char *argv[] = {executable_name.data(), graph_file.data()};

	femg::spectral_problem problem(case_data.n_modes);
	problem.init(2, argv);
	problem.set_coefficients();
	problem.assembly();
	problem.solve();

	const femg::Vector combinatorial =
		problem.combinatorial_laplacian_eigenvalues();

	const std::string output_dir =
		"output/spectral/comparison/" + case_data.name;
	std::filesystem::create_directories(output_dir);

	std::ofstream out(output_dir + "/spectral_comparison.csv");
	if (!out) {
		throw std::runtime_error("Unable to open comparison output file.");
	}

	out << "index,combinatorial_laplacian,metric_laplacian\n";

	const std::size_t rows = std::min(
		static_cast<std::size_t>(combinatorial.size()),
		case_data.n_modes);

	for (std::size_t i = 0; i < rows; ++i) {
		out << i << ","
				<< std::setprecision(16)
				<< combinatorial(static_cast<Eigen::Index>(i)) << ","
				<< problem.eigenvalue(i) << "\n";
	}

	std::cout << "Wrote " << output_dir
					<< "/spectral_comparison.csv\n";
}

} // namespace

int main() {
	try {
		std::vector<CaseData> cases = {
			{"star", "data/star_4.txt", 5},
			{"graphene", "data/graphene_13.txt", 12}
		};

		const std::vector<TreeVariant> tree_variants = {
			{"tree_fixed_length_varying_angle",
			 AngleLaw::InverseLevel, pi / 4.0, LengthLaw::Fixed},
			{"tree_angle_pi4_length_inv",
			 AngleLaw::Fixed, pi / 4.0, LengthLaw::InverseLevel},
			{"tree_angle_pi6_length_inv",
			 AngleLaw::Fixed, pi / 6.0, LengthLaw::InverseLevel},
			{"tree_angle_pi4_length_inv2",
			 AngleLaw::Fixed, pi / 4.0, LengthLaw::InverseLevelSquared},
			{"tree_angle_pi6_length_inv2",
			 AngleLaw::Fixed, pi / 6.0, LengthLaw::InverseLevelSquared},
			{"tree_varying_angle_length_inv",
			 AngleLaw::InverseLevel, pi / 4.0, LengthLaw::InverseLevel},
			{"tree_varying_angle_length_inv2",
			 AngleLaw::InverseLevel, pi / 4.0, LengthLaw::InverseLevelSquared}
		};

		for (const auto &variant : tree_variants) {
			cases.push_back({variant.name, write_tree_variant(variant), 16});
		}

		for (const auto &case_data : cases) {
			write_comparison(case_data);
		}

		return 0;
	} catch (const std::exception &error) {
		std::cerr << "Error: " << error.what() << "\n";
		return 1;
	}
}