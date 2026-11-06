#pragma once

#include "rvre/core/types.hpp"
#include "rvre/pricing/option.hpp"

#include <cstdint>

namespace rvre::pricing {

struct BlackScholesGreeks {
    double delta;
    double gamma;
    double vega;
    double theta;
};

struct ImpliedVolatilityResult {
    core::Volatility volatility;
    std::uint32_t iterations;
    double price_residual;
};

struct BlackScholesMonteCarloResult {
    core::OptionPrice price;
    core::OptionPrice standard_error;
    core::OptionPrice confidence_interval_lower;
    core::OptionPrice confidence_interval_upper;
    std::uint64_t path_count;
    std::uint64_t seed;
};

[[nodiscard]] double black_scholes_price(core::AssetPrice initial_asset,
                                         double strike,
                                         core::Time maturity,
                                         double total_variance,
                                         core::DiscountFactor risk_free_discount,
                                         core::DiscountFactor dividend_discount,
                                         OptionRight right);

[[nodiscard]] BlackScholesGreeks black_scholes_greeks(
    core::AssetPrice initial_asset,
    double strike,
    core::Time maturity,
    core::Volatility volatility,
    core::InterestRate risk_free_rate,
    core::InterestRate dividend_yield,
    OptionRight right);

[[nodiscard]] ImpliedVolatilityResult black_scholes_implied_volatility(
    double option_price,
    core::AssetPrice initial_asset,
    double strike,
    core::Time maturity,
    core::DiscountFactor risk_free_discount,
    core::DiscountFactor dividend_discount,
    OptionRight right,
    double price_tolerance = 1.0e-10,
    std::uint32_t maximum_iterations = 200U);

[[nodiscard]] BlackScholesMonteCarloResult price_black_scholes_monte_carlo(
    const core::BlackScholesParameters& parameters,
    double strike,
    OptionRight right,
    std::uint64_t path_count,
    std::uint64_t seed,
    double confidence_level = 0.95);

}  // namespace rvre::pricing
