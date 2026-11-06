#pragma once

#include "rvre/core/types.hpp"
#include "rvre/models/heston.hpp"
#include "rvre/models/rough_bergomi.hpp"
#include "rvre/pricing/option.hpp"

#include <cstddef>
#include <cstdint>
#include <optional>
#include <string>
#include <variant>
#include <vector>

namespace rvre::hedging {

struct BlackScholesModel {
    core::AssetPrice initial_asset;
    core::Volatility volatility;
    core::InterestRate risk_free_rate;
    core::InterestRate dividend_yield;
};

struct RoughBergomiModel {
    models::RoughBergomiParameters parameters;
};

using MarketModel =
    std::variant<BlackScholesModel, models::HestonParameters,
                 RoughBergomiModel>;
using HedgeModel =
    std::variant<BlackScholesModel, models::HestonParameters,
                 RoughBergomiModel>;

struct MarketScenario {
    std::string name;
    MarketModel model;
};

struct HedgeStrategy {
    std::string name;
    HedgeModel model;
};

struct ExperimentConfig {
    core::AssetPrice initial_asset;
    double strike;
    core::Time maturity;
    pricing::OptionRight right;
    std::uint64_t path_count{1024U};
    std::uint64_t time_step_count{32U};
    std::vector<std::uint64_t> hedge_every_steps{1U, 4U, 8U};
    std::vector<std::uint64_t> seeds{731927U, 731928U, 731929U};
    double proportional_transaction_cost{0.0};
    double confidence_level{0.95};
    std::uint64_t rough_bergomi_pricing_paths{1024U};
    std::uint64_t rough_bergomi_pricing_steps{32U};
};

struct HedgeIntervalRecord {
    std::uint64_t step;
    core::Time time;
    core::AssetPrice spot;
    double delta_before;
    double delta_after;
    double shares_traded;
    double financing_cash_flow;
    double dividend_cash_flow;
    double transaction_cost;
    double cash_after;
    double shares_after;
    bool rebalanced;
    std::optional<bool> delta_pricer_converged;
};

struct PathAccounting {
    std::string market_scenario;
    std::string hedge_model;
    std::uint64_t hedge_every_steps;
    std::uint64_t seed;
    std::uint64_t path_index;
    core::OptionPrice initial_option_value;
    core::OptionPrice initial_option_value_numerical_error;
    std::string initial_option_value_error_kind;
    core::AssetPrice initial_asset;
    double initial_delta;
    bool initial_delta_pricer_converged;
    std::uint64_t nonconverged_delta_rebalances;
    double initial_cash_after_hedge;
    double initial_transaction_cost;
    double terminal_asset;
    double terminal_cash;
    double terminal_shares;
    double terminal_payoff;
    double hedging_pnl;
    double turnover;
    double transaction_costs;
    std::vector<HedgeIntervalRecord> interval_records;
};

struct HedgingSummary {
    std::string market_scenario;
    std::string hedge_model;
    std::uint64_t hedge_every_steps;
    std::uint64_t path_count;
    double mean_hedging_error;
    double standard_deviation;
    double root_mean_square_error;
    double expected_absolute_error;
    double fifth_percentile_pnl;
    double expected_shortfall_five_percent_loss;
    double mean_turnover;
    double mean_transaction_costs;
    double mean_standard_error;
    double confidence_interval_lower;
    double confidence_interval_upper;
};

struct HedgingExperimentResult {
    ExperimentConfig config;
    std::vector<HedgingSummary> summaries;
    std::vector<PathAccounting> path_accounting;
};

[[nodiscard]] HedgingExperimentResult run_hedging_experiment(
    const std::vector<MarketScenario>& scenarios,
    const std::vector<HedgeStrategy>& strategies,
    const ExperimentConfig& config);

void write_hedging_experiment_csv(const HedgingExperimentResult& result,
                                 const std::string& path);

[[nodiscard]] const char* market_model_name(const MarketModel& model) noexcept;
[[nodiscard]] const char* hedge_model_name(const HedgeModel& model) noexcept;

}  // namespace rvre::hedging
