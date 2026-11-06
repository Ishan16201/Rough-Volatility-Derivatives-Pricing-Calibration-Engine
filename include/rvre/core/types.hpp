#pragma once

#include <cstdint>
#include <optional>
#include <string>

namespace rvre::core {

struct AssetPrice {
    explicit constexpr AssetPrice(double value_in_currency) noexcept
        : value(value_in_currency) {}

    double value;
};

struct OptionPrice {
    explicit constexpr OptionPrice(double value_in_currency) noexcept
        : value(value_in_currency) {}

    double value;
};

struct DiscountFactor {
    explicit constexpr DiscountFactor(double value) noexcept : value{value} {}

    double value;
};

struct Variance {
    explicit constexpr Variance(double value_per_year) noexcept
        : value(value_per_year) {}

    double value;
};

struct Volatility {
    explicit constexpr Volatility(double value_per_sqrt_year) noexcept
        : value(value_per_sqrt_year) {}

    double value;
};

struct Time {
    explicit constexpr Time(double value_in_years) noexcept
        : value(value_in_years) {}

    double value;
};

struct InterestRate {
    explicit constexpr InterestRate(double value_per_year) noexcept
        : value(value_per_year) {}

    double value;
};

struct BlackScholesParameters {
    AssetPrice initial_asset;
    Volatility volatility;
    InterestRate risk_free_rate;
    InterestRate dividend_yield;
    Time maturity;
};

struct ModelParameters {
    BlackScholesParameters black_scholes;
};

struct MonteCarloConfig {
    std::uint64_t seed;
    std::uint64_t path_count;
    std::uint64_t time_step_count;
};

struct ExperimentConfig {
    std::string name;
    ModelParameters model;
    MonteCarloConfig monte_carlo;
};

struct PriceEstimate {
    OptionPrice mean;
    OptionPrice standard_error;
};

enum class SimulationStatus {
    configured_only,
    completed
};

struct SimulationResult {
    SimulationStatus status;
    std::uint64_t seed;
    std::uint64_t path_count;
    std::uint64_t time_step_count;
    std::optional<PriceEstimate> price_estimate;
};

}  // namespace rvre::core
