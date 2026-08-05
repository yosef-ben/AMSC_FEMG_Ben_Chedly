#ifndef FEMG_STAR_RADIAL_DECAY_FUNCTIONS_HPP
#define FEMG_STAR_RADIAL_DECAY_FUNCTIONS_HPP

#include "parabolic_problem.hpp"

#include <cmath>
#include <cstddef>

constexpr double star_radial_decay_reaction = 1.0;

inline double star_radial_decay_lambda() {
	return 0.25 * M_PI * M_PI;
}

inline double star_radial_decay_value(double s, double time) {
	const double decay =
		std::exp(-(star_radial_decay_lambda()
			+ star_radial_decay_reaction) * time);
	return decay * std::cos(0.5 * M_PI * s);
}

class StarRadialDecayInitialCondition
	: public femg::parabolic_problem::FunctionU0 {
public:
	double value(std::size_t, double s, double time) const override {
		return star_radial_decay_value(s, time);
	}
};

class StarRadialDecayExactSolution
	: public femg::parabolic_problem::ExactSolution {
public:
	double value(std::size_t, double s, double time) const override {
		return star_radial_decay_value(s, time);
	}
};

#endif
