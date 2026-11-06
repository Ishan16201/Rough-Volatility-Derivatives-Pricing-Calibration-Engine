#pragma once

#include "rvre/core/types.hpp"
#include "rvre/pricing/option.hpp"

#include <cstdint>

namespace rvre::models {

struct HestonParameters {
    core::AssetPrice initial_asset;
    core::Variance initial_variance;
    double mean_reversion_speed;
    core::Variance long_run_variance;
    double volatility_of_variance;
    double correlation;
    core::InterestRate risk_free_rate;
    core::InterestRate dividend_yield;
};

struct HestonCharacteristicFunctionConfig {
    double integration_upper_bound{200.0};
    double absolute_tolerance{1.0e-9};
    std::uint32_t maximum_recursion_depth{20U};
};

struct HestonCharacteristicFunctionDiagnostics {
    double estimated_quadrature_error;
    std::uint64_t integrand_evaluation_count;
    double integration_upper_bound;
    bool converged;
};

struct HestonCharacteristicFunctionResult {
    core::OptionPrice price;
    HestonCharacteristicFunctionDiagnostics diagnostics;
};

enum class HestonVarianceScheme {
    full_truncation_euler
};

struct HestonMonteCarloConfig {
    std::uint64_t path_count;
    std::uint64_t time_step_count;
    std::uint64_t seed;
    double confidence_level{0.95};
    HestonVarianceScheme variance_scheme{
        HestonVarianceScheme::full_truncation_euler};
};

struct HestonMonteCarloDiagnostics {
    core::AssetPrice discounted_forward_mean;
    core::OptionPrice discounted_forward_standard_error;
    double fraction_of_variance_updates_at_zero;
    double minimum_variance;
    double maximum_variance;
};

struct HestonMonteCarloResult {
    core::OptionPrice price;
    core::OptionPrice standard_error;
    core::OptionPrice confidence_interval_lower;
    core::OptionPrice confidence_interval_upper;
    std::uint64_t path_count;
    std::uint64_t seed;
    std::uint64_t time_step_count;
    HestonMonteCarloDiagnostics diagnostics;
};

void validate_heston_parameters(const HestonParameters& parameters);

[[nodiscard]] HestonCharacteristicFunctionResult price_heston_characteristic_function(
    const HestonParameters& parameters,
    double strike,
    core::Time maturity,
    pricing::OptionRight right,
    const HestonCharacteristicFunctionConfig& config = {});

[[nodiscard]] HestonMonteCarloResult price_heston_monte_carlo(
    const HestonParameters& parameters,
    double strike,
    core::Time maturity,
    pricing::OptionRight right,
    const HestonMonteCarloConfig& config);

}  // namespace rvre::models
