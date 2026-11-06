#pragma once

#include "rvre/core/types.hpp"
#include "rvre/models/rough_bergomi.hpp"
#include "rvre/pricing/option.hpp"

#include <cstddef>
#include <cstdint>
#include <vector>

namespace rvre::pricing {

enum class VarianceReduction {
    standard,
    antithetic,
    stock_control_variate
};

struct EuropeanOption {
    double strike;
    core::Time maturity;
    OptionRight right;
};

struct MonteCarloPricingConfig {
    std::uint64_t path_count;
    std::uint64_t batch_size;
    std::uint64_t time_step_count;
    std::uint64_t pilot_path_count;
    std::uint64_t seed;
    double confidence_level;
    models::RoughBergomiConfig model;
};

struct ConfidenceInterval {
    double lower;
    double upper;
};

struct PriceEstimate {
    EuropeanOption option;
    VarianceReduction method;
    core::OptionPrice mean;
    core::OptionPrice standard_error;
    ConfidenceInterval confidence_interval;
    std::uint64_t path_count;
    std::uint64_t independent_observation_count;
    double control_coefficient;
};

struct PricingResult {
    std::vector<PriceEstimate> estimates;
    double confidence_level;
    std::uint64_t seed;
    std::uint64_t path_count;
    std::uint64_t batch_size;
    std::uint64_t pilot_path_count;
    std::uint64_t time_step_count;
};

[[nodiscard]] double black_scholes_price(core::AssetPrice initial_asset,
                                         double strike,
                                         core::Time maturity,
                                         double total_variance,
                                         core::DiscountFactor risk_free_discount,
                                         core::DiscountFactor dividend_discount,
                                         OptionRight right);

[[nodiscard]] PricingResult price_rough_bergomi(
    const models::RoughBergomiParameters& parameters,
    const MonteCarloPricingConfig& config,
    const std::vector<double>& strikes,
    const std::vector<core::Time>& maturities,
    const std::vector<VarianceReduction>& methods);

[[nodiscard]] const char* to_string(OptionRight right) noexcept;
[[nodiscard]] const char* to_string(VarianceReduction method) noexcept;

}  // namespace rvre::pricing
