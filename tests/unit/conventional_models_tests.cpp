#include "rvre/models/heston.hpp"
#include "rvre/pricing/black_scholes.hpp"

#include <catch2/catch_approx.hpp>
#include <catch2/catch_test_macros.hpp>

#include <cmath>
#include <stdexcept>

namespace {

using namespace rvre;

[[nodiscard]] models::HestonParameters heston_parameters(
    double volatility_of_variance = 0.3) {
    return models::HestonParameters{
        .initial_asset = core::AssetPrice{100.0},
        .initial_variance = core::Variance{0.04},
        .mean_reversion_speed = 2.0,
        .long_run_variance = core::Variance{0.04},
        .volatility_of_variance = volatility_of_variance,
        .correlation = -0.7,
        .risk_free_rate = core::InterestRate{0.05},
        .dividend_yield = core::InterestRate{0.0}};
}

}  // namespace

TEST_CASE("Black-Scholes price, Greeks, and implied volatility are consistent") {
    const core::AssetPrice spot{100.0};
    constexpr double strike = 100.0;
    constexpr double maturity = 1.0;
    constexpr double volatility = 0.2;
    const core::InterestRate rate{0.05};
    const core::InterestRate dividend{0.02};
    const double discount_r = std::exp(-rate.value * maturity);
    const double discount_q = std::exp(-dividend.value * maturity);

    const double call = pricing::black_scholes_price(
        spot, strike, core::Time{maturity},
        volatility * volatility * maturity,
        core::DiscountFactor{discount_r}, core::DiscountFactor{discount_q},
        pricing::OptionRight::call);
    const double put = pricing::black_scholes_price(
        spot, strike, core::Time{maturity},
        volatility * volatility * maturity,
        core::DiscountFactor{discount_r}, core::DiscountFactor{discount_q},
        pricing::OptionRight::put);
    REQUIRE(call == Catch::Approx(9.227005508).margin(2.0e-8));
    REQUIRE(call - put ==
            Catch::Approx(spot.value * discount_q - strike * discount_r)
                .margin(1.0e-12));

    const auto greeks = pricing::black_scholes_greeks(
        spot, strike, core::Time{maturity}, core::Volatility{volatility}, rate,
        dividend, pricing::OptionRight::call);
    REQUIRE(greeks.delta == Catch::Approx(0.586851146).margin(2.0e-8));
    REQUIRE(greeks.gamma == Catch::Approx(0.018950579).margin(2.0e-8));
    REQUIRE(greeks.vega == Catch::Approx(37.901157).margin(2.0e-5));
    REQUIRE(greeks.theta == Catch::Approx(-5.089318).margin(2.0e-5));
    const auto put_greeks = pricing::black_scholes_greeks(
        spot, strike, core::Time{maturity}, core::Volatility{volatility}, rate,
        dividend, pricing::OptionRight::put);
    REQUIRE(greeks.delta - put_greeks.delta ==
            Catch::Approx(discount_q).margin(1.0e-12));
    REQUIRE(put_greeks.gamma == Catch::Approx(greeks.gamma).margin(1.0e-12));
    REQUIRE(put_greeks.vega == Catch::Approx(greeks.vega).margin(1.0e-12));

    const auto implied = pricing::black_scholes_implied_volatility(
        call, spot, strike, core::Time{maturity},
        core::DiscountFactor{discount_r}, core::DiscountFactor{discount_q},
        pricing::OptionRight::call);
    REQUIRE(implied.volatility.value == Catch::Approx(volatility).margin(1.0e-9));
    REQUIRE(std::abs(implied.price_residual) < 1.0e-10);
    const double put_iv = pricing::black_scholes_implied_volatility(
        put, spot, strike, core::Time{maturity},
        core::DiscountFactor{discount_r}, core::DiscountFactor{discount_q},
        pricing::OptionRight::put)
                              .volatility.value;
    REQUIRE(put_iv == Catch::Approx(volatility).margin(1.0e-9));
}

TEST_CASE("Black-Scholes handles no-arbitrage and zero-volatility boundaries") {
    REQUIRE(pricing::black_scholes_price(
                core::AssetPrice{100.0}, 120.0, core::Time{1.0}, 0.0,
                core::DiscountFactor{1.0}, core::DiscountFactor{1.0},
                pricing::OptionRight::call) == 0.0);
    REQUIRE_THROWS_AS(
        pricing::black_scholes_implied_volatility(
            101.0, core::AssetPrice{100.0}, 100.0, core::Time{1.0},
            core::DiscountFactor{1.0}, core::DiscountFactor{1.0},
            pricing::OptionRight::call),
        std::invalid_argument);
    REQUIRE_THROWS_AS(
        pricing::black_scholes_greeks(
            core::AssetPrice{100.0}, 100.0, core::Time{0.0},
            core::Volatility{0.2}, core::InterestRate{0.03},
            core::InterestRate{0.0}, pricing::OptionRight::call),
        std::invalid_argument);
}

TEST_CASE("Black-Scholes terminal Monte Carlo agrees with its analytical value") {
    const core::BlackScholesParameters parameters{
        .initial_asset = core::AssetPrice{100.0},
        .volatility = core::Volatility{0.25},
        .risk_free_rate = core::InterestRate{0.03},
        .dividend_yield = core::InterestRate{0.01},
        .maturity = core::Time{1.25}};
    const double analytic = pricing::black_scholes_price(
        parameters.initial_asset, 105.0, parameters.maturity,
        parameters.volatility.value * parameters.volatility.value *
            parameters.maturity.value,
        core::DiscountFactor{std::exp(-parameters.risk_free_rate.value *
                                      parameters.maturity.value)},
        core::DiscountFactor{std::exp(-parameters.dividend_yield.value *
                                      parameters.maturity.value)},
        pricing::OptionRight::call);
    const auto estimate = pricing::price_black_scholes_monte_carlo(
        parameters, 105.0, pricing::OptionRight::call, 50000U, 84217U);
    REQUIRE(std::abs(estimate.price.value - analytic) <
            5.0 * estimate.standard_error.value);
    REQUIRE(estimate.path_count == 50000U);
    REQUIRE(estimate.confidence_interval_lower.value <
            estimate.price.value);
    REQUIRE(estimate.confidence_interval_upper.value >
            estimate.price.value);
}

TEST_CASE("Heston validates boundaries and handles the zero-vol-of-variance limit") {
    auto parameters = heston_parameters(0.0);
    REQUIRE_NOTHROW(models::validate_heston_parameters(parameters));
    const auto characteristic = models::price_heston_characteristic_function(
        parameters, 100.0, core::Time{1.0}, pricing::OptionRight::call);
    const double analytic = pricing::black_scholes_price(
        parameters.initial_asset, 100.0, core::Time{1.0}, 0.04,
        core::DiscountFactor{std::exp(-0.05)},
        core::DiscountFactor{1.0}, pricing::OptionRight::call);
    REQUIRE(characteristic.price.value == Catch::Approx(analytic).margin(1.0e-12));
    REQUIRE(characteristic.diagnostics.converged);
    REQUIRE(characteristic.diagnostics.integrand_evaluation_count == 0U);

    parameters.correlation = 1.01;
    REQUIRE_THROWS_AS(models::validate_heston_parameters(parameters),
                      std::invalid_argument);
    parameters = heston_parameters();
    parameters.mean_reversion_speed = 0.0;
    REQUIRE_THROWS_AS(models::validate_heston_parameters(parameters),
                      std::invalid_argument);
    parameters = heston_parameters();
    parameters.initial_variance = core::Variance{-0.01};
    REQUIRE_THROWS_AS(models::validate_heston_parameters(parameters),
                      std::invalid_argument);
}

TEST_CASE("Heston characteristic function respects independent benchmark and parity") {
    const auto parameters = heston_parameters();
    const auto call = models::price_heston_characteristic_function(
        parameters, 100.0, core::Time{1.0}, pricing::OptionRight::call);
    const auto put = models::price_heston_characteristic_function(
        parameters, 100.0, core::Time{1.0}, pricing::OptionRight::put);
    REQUIRE(call.diagnostics.converged);
    REQUIRE(call.diagnostics.integrand_evaluation_count > 0U);
    REQUIRE(call.price.value == Catch::Approx(10.394).margin(0.01));
    REQUIRE(call.price.value - put.price.value ==
            Catch::Approx(100.0 - 100.0 * std::exp(-0.05))
                .margin(1.0e-7));
}

TEST_CASE("Heston full-truncation Monte Carlo is reproducible and agrees with characteristic pricing") {
    const auto parameters = heston_parameters();
    const auto config = models::HestonMonteCarloConfig{
        .path_count = 30000U,
        .time_step_count = 96U,
        .seed = 71223U,
        .confidence_level = 0.95};
    const auto first = models::price_heston_monte_carlo(
        parameters, 100.0, core::Time{1.0}, pricing::OptionRight::call, config);
    const auto replay = models::price_heston_monte_carlo(
        parameters, 100.0, core::Time{1.0}, pricing::OptionRight::call, config);
    const auto characteristic = models::price_heston_characteristic_function(
        parameters, 100.0, core::Time{1.0}, pricing::OptionRight::call);
    REQUIRE(first.price.value == replay.price.value);
    REQUIRE(first.standard_error.value == replay.standard_error.value);
    REQUIRE(first.diagnostics.minimum_variance < 0.0);
    REQUIRE(first.diagnostics.fraction_of_variance_updates_at_zero > 0.0);
    REQUIRE(std::abs(first.diagnostics.discounted_forward_mean.value - 100.0) <
            5.0 * first.diagnostics.discounted_forward_standard_error.value);
    REQUIRE(std::abs(first.price.value - characteristic.price.value) <
            5.0 * first.standard_error.value + 0.15);
}

TEST_CASE("Heston full truncation remains defined when the Feller condition fails") {
    auto parameters = heston_parameters(0.9);
    REQUIRE_NOTHROW(models::validate_heston_parameters(parameters));
    const auto result = models::price_heston_monte_carlo(
        parameters, 100.0, core::Time{1.0}, pricing::OptionRight::call,
        models::HestonMonteCarloConfig{5000U, 48U, 441U, 0.95});
    REQUIRE(std::isfinite(result.price.value));
    REQUIRE(result.diagnostics.minimum_variance < 0.0);
    REQUIRE(result.diagnostics.fraction_of_variance_updates_at_zero > 0.0);
}

TEST_CASE("Heston path-count and time-step changes expose sampling and discretization") {
    const auto parameters = heston_parameters();
    const auto coarse = models::price_heston_monte_carlo(
        parameters, 100.0, core::Time{1.0}, pricing::OptionRight::call,
        models::HestonMonteCarloConfig{6000U, 24U, 9131U, 0.95});
    const auto fine = models::price_heston_monte_carlo(
        parameters, 100.0, core::Time{1.0}, pricing::OptionRight::call,
        models::HestonMonteCarloConfig{6000U, 96U, 9131U, 0.95});
    const auto more_paths = models::price_heston_monte_carlo(
        parameters, 100.0, core::Time{1.0}, pricing::OptionRight::call,
        models::HestonMonteCarloConfig{24000U, 96U, 9131U, 0.95});
    REQUIRE(coarse.time_step_count == 24U);
    REQUIRE(fine.time_step_count == 96U);
    REQUIRE(more_paths.standard_error.value < fine.standard_error.value);
    REQUIRE(std::abs(fine.price.value - more_paths.price.value) <
            5.0 * std::hypot(fine.standard_error.value,
                             more_paths.standard_error.value) +
                0.1);
}
