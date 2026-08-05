#ifndef FEMG_STAR_LOCALIZED_HEAT_FUNCTIONS_HPP
#define FEMG_STAR_LOCALIZED_HEAT_FUNCTIONS_HPP

#include "parabolic_problem.hpp"

#include <cmath>
#include <cstddef>

constexpr double star_localized_heat_reaction = 0.0;
constexpr double star_localized_constant = 4.0;
constexpr double star_localized_amplitude = 2.0;

inline double star_localized_lambda() {
	return 0.25 * M_PI * M_PI;
}

inline double star_localized_weight(std::size_t edge_index) {
	if (edge_index == 0) {
		return 3.0;
	}
	return -1.0;
}

inline double star_localized_heat_value(
	std::size_t edge_index,
	double s,
	double time,
	double diffusion) {
	const double decay = std::exp(
		-diffusion * star_localized_lambda() * time);
	return star_localized_constant
		+ star_localized_amplitude
			* star_localized_weight(edge_index)
			* decay
			* std::sin(0.5 * M_PI * s);
}

class StarLocalizedHeatInitialCondition
	: public femg::parabolic_problem::FunctionU0 {
public:
	explicit StarLocalizedHeatInitialCondition(double diffusion)
		: diffusion_(diffusion) {}

	double value(std::size_t edge_index, double s, double) const override {
		return star_localized_heat_value(edge_index, s, 0.0, diffusion_);
	}

private:
	double diffusion_ = 1.0;
};

class StarLocalizedHeatExactSolution
	: public femg::parabolic_problem::ExactSolution {
public:
	explicit StarLocalizedHeatExactSolution(double diffusion)
		: diffusion_(diffusion) {}

	double value(std::size_t edge_index, double s, double time) const override {
		return star_localized_heat_value(edge_index, s, time, diffusion_);
	}

private:
	double diffusion_ = 1.0;
};

inline double star_localized_boundary_value(std::size_t, double) {
	return 0.0;
}

#endif
