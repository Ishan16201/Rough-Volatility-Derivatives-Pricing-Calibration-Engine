#include "rvre/hedging/experiment.hpp"
#include "rvre/pricing/black_scholes.hpp"
#include "rvre/risk/sensitivities.hpp"

#include <catch2/catch_approx.hpp>
#include <catch2/catch_test_macros.hpp>

#include <algorithm>
#include <cmath>
#include <cstdint>
#include <cstdio>
#include <fstream>
#include <string>
#include <vector>

namespace {

using namespace rvre;

TEST_CASE("Finite-difference Greeks refine against analytical Black-Scholes") {
    const risk::SensitivityState state{
        core::AssetPrice{100.0},
        100.0,
        core::Time{1.0},
        core::Volatility{0.25},
        core::InterestRate{0.03},
        core::InterestRate{0.01},
        models::HurstExponent{0.1},
        models::VolatilityOfVolatility{1.2},
        {core::Variance{0.04}, core::Variance{0.05}}};
    std::vector<std::uint64_t> observed_streams;
    const risk::OptionValueFunction price =
        [&observed_streams](const risk::SensitivityState& point,
                            std::uint64_t,
                            std::uint64_t stream) {
            observed_streams.push_back(stream);
            return core::OptionPrice{pricing::black_scholes_price(
                point.spot, point.strike, point.maturity,
                point.volatility.value * point.volatility.value *
                    point.maturity.value,
                core::DiscountFactor{std::exp(
                    -point.risk_free_rate.value * point.maturity.value)},
                core::DiscountFactor{std::exp(
                    -point.dividend_yield.value * point.maturity.value)},
                pricing::OptionRight::call)};
        };
    const auto estimates = risk::estimate_sensitivities(
        state, pricing::OptionRight::call,
        {risk::RiskFactor::delta, risk::RiskFactor::gamma,
         risk::RiskFactor::vega, risk::RiskFactor::theta,
         risk::RiskFactor::rho},
        risk::FiniteDifferenceConfig{.relative_step = 0.01,
                                     .absolute_step = 1.0e-3,
                                     .convergence_tolerance = 0.08,
                                     .refinement_count = 5U},
        price);
    const auto analytic = pricing::black_scholes_greeks(
        state.spot, state.strike, state.maturity, state.volatility,
        state.risk_free_rate, state.dividend_yield,
        pricing::OptionRight::call);
    REQUIRE(estimates.size() == 5U);
    REQUIRE(estimates[0].value == Catch::Approx(analytic.delta).margin(2.0e-4));
    REQUIRE(estimates[1].value == Catch::Approx(analytic.gamma).margin(2.0e-4));
    REQUIRE(estimates[2].value == Catch::Approx(analytic.vega).margin(2.0e-3));
    REQUIRE(estimates[3].value == Catch::Approx(analytic.theta).margin(2.0e-3));
    REQUIRE(estimates[4].value ==
            Catch::Approx(risk::black_scholes_rho(
                              state.spot, state.strike, state.maturity,
                              state.volatility, state.risk_free_rate,
                              state.dividend_yield,
                              pricing::OptionRight::call))
                .margin(2.0e-3));
    constexpr double rate_bump = 1.0e-5;
    const auto price_at_rate = [&state](double rate) {
        return pricing::black_scholes_price(
            state.spot, state.strike, state.maturity,
            state.volatility.value * state.volatility.value *
                state.maturity.value,
            core::DiscountFactor{
                std::exp(-rate * state.maturity.value)},
            core::DiscountFactor{std::exp(
                -state.dividend_yield.value * state.maturity.value)},
            pricing::OptionRight::call);
    };
    REQUIRE(risk::black_scholes_rho(
                state.spot, state.strike, state.maturity, state.volatility,
                state.risk_free_rate, state.dividend_yield,
                pricing::OptionRight::call) ==
            Catch::Approx((price_at_rate(state.risk_free_rate.value +
                                         rate_bump) -
                           price_at_rate(state.risk_free_rate.value -
                                         rate_bump)) /
                          (2.0 * rate_bump))
                .margin(2.0e-5));
    REQUIRE(observed_streams.size() > 5U);
    for (std::size_t index = 0U; index < 10U; index += 2U) {
        REQUIRE(observed_streams[index] == observed_streams[index + 1U]);
    }
    for (std::size_t index = 10U; index < 25U; index += 3U) {
        REQUIRE(observed_streams[index] == observed_streams[index + 1U]);
        REQUIRE(observed_streams[index] == observed_streams[index + 2U]);
    }
    for (std::size_t index = 25U; index + 1U < observed_streams.size();
         index += 2U) {
        REQUIRE(observed_streams[index] == observed_streams[index + 1U]);
    }
}

TEST_CASE("Discrete hedge reports full cash, turnover, and cost accounting") {
    const core::AssetPrice spot{100.0};
    const hedging::BlackScholesModel model{
        spot, core::Volatility{0.2}, core::InterestRate{0.03},
        core::InterestRate{0.01}};
    const hedging::ExperimentConfig config{
        spot,
        100.0,
        core::Time{1.0},
        pricing::OptionRight::call,
        16U,
        8U,
        {1U, 4U},
        {101U, 202U},
        0.001,
        0.95,
        16U,
        8U};
    const std::vector<hedging::MarketScenario> scenarios{
        {"BS test market", model}};
    const std::vector<hedging::HedgeStrategy> strategies{
        {"BS test hedge", model}};
    const auto result =
        hedging::run_hedging_experiment(scenarios, strategies, config);
    const auto replay =
        hedging::run_hedging_experiment(scenarios, strategies, config);
    REQUIRE(result.summaries.size() == 2U);
    REQUIRE(result.path_accounting.size() == 64U);
    REQUIRE(result.summaries[0].mean_hedging_error ==
            replay.summaries[0].mean_hedging_error);
    REQUIRE(result.path_accounting.front().transaction_costs ==
            Catch::Approx(0.001 *
                          result.path_accounting.front().turnover));
    REQUIRE(result.path_accounting.front().initial_delta_pricer_converged);
    REQUIRE(result.path_accounting.front().nonconverged_delta_rebalances == 0U);
    REQUIRE(result.path_accounting.front().hedging_pnl ==
            Catch::Approx(
                result.path_accounting.front().terminal_cash +
                result.path_accounting.front().terminal_shares *
                    result.path_accounting.front().terminal_asset -
                result.path_accounting.front().terminal_payoff));
    REQUIRE(result.path_accounting.front().interval_records.size() == 8U);
    REQUIRE(result.path_accounting.front().interval_records.back().rebalanced ==
            false);
    REQUIRE(result.path_accounting.front().interval_records[0].rebalanced);
    REQUIRE(result.summaries.front().confidence_interval_lower <=
            result.summaries.front().confidence_interval_upper);
    REQUIRE(result.summaries.front().root_mean_square_error >=
            std::abs(result.summaries.front().mean_hedging_error));
    const std::string report_path = "/tmp/rvre-hedging-accounting.csv";
    hedging::write_hedging_experiment_csv(result, report_path);
    std::ifstream report{report_path};
    std::string line;
    REQUIRE(static_cast<bool>(std::getline(report, line)));
    const std::size_t column_count =
        static_cast<std::size_t>(std::count(line.begin(), line.end(), ',')) + 1U;
    REQUIRE(column_count == 54U);
    while (std::getline(report, line)) {
        REQUIRE(static_cast<std::size_t>(
                    std::count(line.begin(), line.end(), ',')) + 1U ==
                column_count);
    }
    report.close();
    std::remove(report_path.c_str());
}

TEST_CASE("All hedge strategies use the same market option mark") {
    const core::AssetPrice spot{100.0};
    const hedging::BlackScholesModel market{
        spot, core::Volatility{0.25}, core::InterestRate{0.02},
        core::InterestRate{0.01}};
    const hedging::ExperimentConfig config{
        spot,
        100.0,
        core::Time{0.5},
        pricing::OptionRight::call,
        4U,
        2U,
        {1U},
        {91U},
        0.0,
        0.95,
        8U,
        4U};
    const std::vector<hedging::MarketScenario> scenarios{
        {"common mark market", market}};
    const std::vector<hedging::HedgeStrategy> strategies{
        {"market vol hedge", market},
        {"misspecified hedge",
         hedging::BlackScholesModel{
             spot, core::Volatility{0.4}, core::InterestRate{0.02},
             core::InterestRate{0.01}}}};

    const auto result =
        hedging::run_hedging_experiment(scenarios, strategies, config);
    REQUIRE(result.path_accounting.size() == 8U);
    for (std::size_t path = 0U; path < 4U; ++path) {
        REQUIRE(result.path_accounting[path].initial_option_value.value ==
                result.path_accounting[path + 4U].initial_option_value.value);
        REQUIRE(result.path_accounting[path].initial_option_value.value ==
                Catch::Approx(pricing::black_scholes_price(
                    spot, config.strike, config.maturity,
                    market.volatility.value * market.volatility.value *
                        config.maturity.value,
                    core::DiscountFactor{std::exp(
                        -market.risk_free_rate.value * config.maturity.value)},
                    core::DiscountFactor{std::exp(
                        -market.dividend_yield.value * config.maturity.value)},
                    config.right)));
    }
    REQUIRE(result.path_accounting[0].initial_delta !=
            result.path_accounting[4].initial_delta);
}

TEST_CASE("Rough-model parameter sensitivities use explicit curve knots") {
    const risk::SensitivityState state{
        core::AssetPrice{100.0},
        100.0,
        core::Time{1.0},
        core::Volatility{0.2},
        core::InterestRate{0.03},
        core::InterestRate{0.01},
        models::HurstExponent{0.12},
        models::VolatilityOfVolatility{1.1},
        {core::Variance{0.04}, core::Variance{0.06}}};
    const risk::OptionValueFunction synthetic_price =
        [](const risk::SensitivityState& point, std::uint64_t,
           std::uint64_t) {
            return core::OptionPrice{
                2.0 * point.hurst.value + 3.0 * point.eta.value +
                4.0 *
                    point.initial_forward_variance_curve[1].value};
        };
    const auto estimates = risk::estimate_sensitivities(
        state, pricing::OptionRight::call,
        {risk::RiskFactor::hurst, risk::RiskFactor::eta,
         risk::RiskFactor::initial_forward_variance},
        risk::FiniteDifferenceConfig{
            .relative_step = 0.01,
            .absolute_step = 1.0e-4,
            .convergence_tolerance = 1.0e-8,
            .refinement_count = 4U,
            .forward_variance_curve_index = 1U},
        synthetic_price);
    REQUIRE(estimates[0].value == Catch::Approx(2.0).margin(1.0e-8));
    REQUIRE(estimates[1].value == Catch::Approx(3.0).margin(1.0e-8));
    REQUIRE(estimates[2].value == Catch::Approx(4.0).margin(1.0e-8));
    REQUIRE(estimates[0].converged);
    REQUIRE(estimates[1].converged);
    REQUIRE(estimates[2].converged);
}

}  // namespace
