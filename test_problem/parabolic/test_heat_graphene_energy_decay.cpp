#include "parabolic_problem.hpp"
#include "spectral_problem.hpp"

#include <algorithm>
#include <cmath>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

namespace {

struct ModalComponent {
	std::size_t mode = 0;
	double coefficient = 0.0;
	double eigenvalue = 0.0;
	std::vector<std::vector<double>> edge_values;
};

double interpolate_on_edge(
	const std::vector<double> &values,
	double length,
	double s) {
	if (values.empty()) {
		throw std::runtime_error("Missing modal values on an edge.");
	}
	if (values.size() == 1) {
		return values.front();
	}

	const double h = length / static_cast<double>(values.size() - 1);
	const double clamped_s = std::clamp(s, 0.0, length);
	const double position = clamped_s / h;
	const auto cell = static_cast<std::size_t>(
		std::min(
			std::floor(position),
			static_cast<double>(values.size() - 2)));
	const double xi = position - static_cast<double>(cell);

	return (1.0 - xi) * values[cell] + xi * values[cell + 1];
}

class ModalCombinationInitialCondition
	: public femg::parabolic_problem::FunctionU0 {
public:
	ModalCombinationInitialCondition(
		std::vector<ModalComponent> components,
		std::vector<double> edge_lengths)
		: components_(std::move(components)),
			edge_lengths_(std::move(edge_lengths)) {}

	double value(std::size_t edge_index, double s, double) const override {
		double result = 0.0;
		for (const auto &component : components_) {
			result += component.coefficient *
				interpolate_on_edge(
					component.edge_values.at(edge_index),
					edge_lengths_.at(edge_index),
					s);
		}
		return result;
	}

private:
	std::vector<ModalComponent> components_;
	std::vector<double> edge_lengths_;
};

double exact_l2_squared(
	const std::vector<ModalComponent> &components,
	double diffusion,
	double time) {
	double value = 0.0;
	for (const auto &component : components) {
		const double decay =
			std::exp(-2.0 * diffusion * component.eigenvalue * time);
		value += component.coefficient * component.coefficient * decay;
	}
	return value;
}

double exact_energy(
	const std::vector<ModalComponent> &components,
	double diffusion,
	double time) {
	double value = 0.0;
	for (const auto &component : components) {
		const double decay =
			std::exp(-2.0 * diffusion * component.eigenvalue * time);
		value += diffusion * component.eigenvalue *
			component.coefficient * component.coefficient * decay;
	}
	return value;
}

} // namespace

int main(int argc, char *argv[]) {
	try {
		const double diffusion = 1.0;
		const double T = 1.0;
		const double deltat = 0.01;
		const double theta = 0.5;

		std::string graph_file = "data/graphene_13.txt";
		char *local_argv[] = {argv[0], graph_file.data()};
		const int graph_argc = (argc >= 2) ? argc : 2;
		char **graph_argv = (argc >= 2) ? argv : local_argv;

		femg::spectral_problem eigen_problem(10);
		eigen_problem.init(graph_argc, graph_argv);
		eigen_problem.set_coefficients();
		eigen_problem.assembly();
		eigen_problem.solve();

		const std::vector<std::pair<std::size_t, double>> modal_coefficients = {
			{6, 1.0},
			{7, -0.75},
			{8, 0.50},
			{9, -0.35}
		};

		std::vector<ModalComponent> components;
		components.reserve(modal_coefficients.size());
		for (const auto &[mode, coefficient] : modal_coefficients) {
			ModalComponent component;
			component.mode = mode;
			component.coefficient = coefficient;
			component.eigenvalue = eigen_problem.eigenvalue(mode);
			component.edge_values = eigen_problem.edge_nodal_values(mode);
			components.push_back(std::move(component));
		}

		const auto edge_lengths = eigen_problem.edge_lengths();

		femg::parabolic_problem problem(diffusion, T, deltat, theta);
		problem.set_initial_condition(
			std::make_shared<ModalCombinationInitialCondition>(
				components,
				edge_lengths));
		problem.set_output_directory("output/visualization/graphene_energy_decay");
		problem.set_verbose(false);

		problem.init(graph_argc, graph_argv);
		problem.set_coefficients();
		problem.assemble_matrices();
		problem.print_matrix_summary(std::cout);

		std::cout << "Graphene energy-decay benchmark\n";
		for (const auto &component : components) {
			std::cout << "  mode " << component.mode
								<< ", coefficient = " << component.coefficient
								<< ", lambda = " << component.eigenvalue << "\n";
		}

		std::filesystem::create_directories("output/energy");
		std::ofstream energy_out("output/energy/graphene_energy_decay.csv");
		energy_out
			<< "time,l2_squared,energy,exact_l2_squared,exact_energy\n";

		problem.set_time_step_callback(
			[&](std::size_t, double time, const femg::Vector &solution) {
				const double l2_squared =
					solution.dot(problem.mass_matrix() * solution);
				const double energy =
					solution.dot(problem.stiffness_matrix() * solution);
				energy_out << std::setprecision(16)
									 << time << ","
									 << l2_squared << ","
									 << energy << ","
									 << exact_l2_squared(components, diffusion, time) << ","
									 << exact_energy(components, diffusion, time)
									 << "\n";
			});

		problem.solve();

		std::cout << "Energy data written to "
							<< "output/energy/graphene_energy_decay.csv\n";
		return 0;
	} catch (const std::exception &error) {
		std::cerr << "Error: " << error.what() << "\n";
		return 1;
	}
}
