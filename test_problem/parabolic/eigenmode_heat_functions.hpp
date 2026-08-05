#ifndef FEMG_EIGENMODE_HEAT_FUNCTIONS_HPP
#define FEMG_EIGENMODE_HEAT_FUNCTIONS_HPP

#include "parabolic_problem.hpp"

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <stdexcept>
#include <utility>
#include <vector>

class EigenmodeHeatData {
public:
	EigenmodeHeatData(
		std::vector<std::vector<double>> edge_values,
		std::vector<double> edge_lengths,
		double eigenvalue,
		double diffusion)
		: edge_values_(std::move(edge_values)),
			edge_lengths_(std::move(edge_lengths)),
			eigenvalue_(eigenvalue),
			diffusion_(diffusion) {}

	double eigenvalue() const { return eigenvalue_; }
	double diffusion() const { return diffusion_; }

	double phi(std::size_t edge_index, double s) const {
		if (edge_index >= edge_values_.size()) {
			throw std::out_of_range("Invalid edge index in eigenmode data.");
		}

		const auto &values = edge_values_[edge_index];
		if (values.empty()) {
			throw std::runtime_error("Missing edge values in eigenmode data.");
		}
		if (values.size() == 1) {
			return values.front();
		}

		const double length = edge_lengths_[edge_index];
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

	double exact_value(std::size_t edge_index, double s, double time) const {
		return std::exp(-diffusion_ * eigenvalue_ * time) * phi(edge_index, s);
	}

private:
	std::vector<std::vector<double>> edge_values_;
	std::vector<double> edge_lengths_;
	double eigenvalue_ = 0.0;
	double diffusion_ = 1.0;
};

class EigenmodeHeatInitialCondition
	: public femg::parabolic_problem::FunctionU0 {
public:
	explicit EigenmodeHeatInitialCondition(EigenmodeHeatData data)
		: data_(std::move(data)) {}

	double value(std::size_t edge_index, double s, double) const override {
		return data_.phi(edge_index, s);
	}

private:
	EigenmodeHeatData data_;
};

class EigenmodeHeatExactSolution
	: public femg::parabolic_problem::ExactSolution {
public:
	explicit EigenmodeHeatExactSolution(EigenmodeHeatData data)
		: data_(std::move(data)) {}

	double value(std::size_t edge_index, double s, double time) const override {
		return data_.exact_value(edge_index, s, time);
	}

private:
	EigenmodeHeatData data_;
};

#endif
