#ifndef FEMG_STAR_SINE_FUNCTIONS_HPP
#define FEMG_STAR_SINE_FUNCTIONS_HPP

#include "parabolic_problem.hpp"

#include <cmath>
#include <cstddef>

inline double star_sine_value(std::size_t edge_index, double s) {
	if (edge_index == 0) {
		return std::sin(2.0 * M_PI * s);
	}
	if (edge_index == 2) {
		return std::sin(-2.0 * M_PI * s);
	}
	return 0.0;
}

class StarSineInitialCondition
	: public femg::parabolic_problem::FunctionU0 {
public:
	double value(std::size_t edge_index, double s, double) const override {
		return star_sine_value(edge_index, s);
	}
};

class StarSineForcingTerm
	: public femg::parabolic_problem::ForcingTerm {
public:
	double value(std::size_t edge_index, double s, double) const override {
		return 4.0 * M_PI * M_PI * star_sine_value(edge_index, s);
	}
};

class StarSineExactSolution
	: public femg::parabolic_problem::ExactSolution {
public:
	double value(std::size_t edge_index, double s, double) const override {
		return star_sine_value(edge_index, s);
	}
};

#endif
