#include "rvre/pricing/monte_carlo.hpp"

#include <catch2/catch_test_macros.hpp>

#include <cmath>
#include <memory>
#include <stdexcept>
#include <vector>

namespace {

using namespace rvre;

[[nodiscard]] models::RoughBergomiParameters parameters(
    double eta = 1.0e-8) {
    return models::RoughBergomiParameters{
        .initial_asset = core::AssetPrice{100.0},
        .hurst = models::HurstExponent{0.2},
        .eta = models::VolatilityOfVolatility{eta},
        .correlation = -0.5,
        .initial_forward_variance =
            std::make_shared<models::PiecewiseLinearForwardVarianceCurve>(
                std::vector<core::Time>{core::Time{0.0}, core::Time{2.0}},
                std::vector<core::Variance>{core::Variance{0.04},
                                            core::Variance{0.04}}),
        .risk_free_discount_curve =
            std::make_shared<models::PiecewiseLinearDiscountCurve>(
                std::vector<core::Time>{core::Time{0.0}, core::Time{2.0}},
                std::vector<core::DiscountFactor>{
                    core::DiscountFactor{1.0},
                    core::DiscountFactor{std::exp(-0.03 * 2.0)}}),
        .dividend_discount_curve =
            std::make_shared<models::PiecewiseLinearDiscountCurve>(
                std::vector<core::Time>{core::Time{0.0}, core::Time{2.0}},
                std::vector<core::DiscountFactor>{
                    core::DiscountFactor{1.0},
                    core::DiscountFactor{std::exp(-0.01 * 2.0)}})};
}

[[nodiscard]] pricing::MonteCarloPricingConfig config(
    std::uint64_t paths = 12000U) {
    return pricing::MonteCarloPricingConfig{
        .path_count = paths,
        .batch_size = 2000U,
        .time_step_count = 16U,
        .pilot_path_count = 6000U,
        .seed = 1234567U,
        .confidence_level = 0.95,
        .model = models::RoughBergomiConfig{
            .path_count = 2U,
            .seed = 0U,
            .volterra_scheme = simulation::VolterraScheme::cell_average_reference}};
}

}  // namespace

TEST_CASE("Black-Scholes pricing satisfies parity and boundary cases") {
    using pricing::OptionRight;
    const double call = pricing::black_scholes_price(
        core::AssetPrice{100.0}, 100.0, core::Time{1.0}, 0.04,
        core::DiscountFactor{std::exp(-0.03)},
        core::DiscountFactor{std::exp(-0.01)}, OptionRight::call);
    const double put = pricing::black_scholes_price(
        core::AssetPrice{100.0}, 100.0, core::Time{1.0}, 0.04,
        core::DiscountFactor{std::exp(-0.03)},
        core::DiscountFactor{std::exp(-0.01)}, OptionRight::put);

    REQUIRE(call > 0.0);
    REQUIRE(put > 0.0);
    REQUIRE(std::abs((call - put) -
                     (100.0 * std::exp(-0.01) -
                      100.0 * std::exp(-0.03))) < 1.0e-12);
    REQUIRE(pricing::black_scholes_price(
                core::AssetPrice{100.0}, 100.0, core::Time{1.0}, 0.0,
                core::DiscountFactor{1.0}, core::DiscountFactor{1.0},
                OptionRight::call) == 0.0);
    REQUIRE_THROWS_AS(
        pricing::black_scholes_price(
            core::AssetPrice{100.0}, 0.0, core::Time{1.0}, 0.04,
            core::DiscountFactor{1.0}, core::DiscountFactor{1.0},
            OptionRight::call),
        std::invalid_argument);
}

TEST_CASE("Rough Bergomi Monte Carlo converges to Black-Scholes for constant variance") {
    const auto model = parameters();
    const auto run = pricing::price_rough_bergomi(
        model, config(), {100.0}, {core::Time{1.0}},
        {pricing::VarianceReduction::standard});
    REQUIRE(run.estimates.size() == 2U);
    REQUIRE(run.seed == 1234567U);
    REQUIRE(run.path_count == 12000U);
    REQUIRE(run.estimates[0].method == pricing::VarianceReduction::standard);
    REQUIRE(run.estimates[0].independent_observation_count == run.path_count);

    const double analytic = pricing::black_scholes_price(
        model.initial_asset, 100.0, core::Time{1.0}, 0.04,
        model.risk_free_discount_curve->discount_factor(core::Time{1.0}),
        model.dividend_discount_curve->discount_factor(core::Time{1.0}),
        pricing::OptionRight::call);
    const auto& estimate = run.estimates[0];
    REQUIRE(std::abs(estimate.mean.value - analytic) <
            5.0 * estimate.standard_error.value + 0.03);
    REQUIRE(estimate.confidence_interval.lower <= estimate.mean.value);
    REQUIRE(estimate.confidence_interval.upper >= estimate.mean.value);
    REQUIRE(estimate.mean.value >= 0.0);
}

TEST_CASE("Constant-variance prices remain consistent under grid refinement") {
    const auto model = parameters();
    auto coarse_config = config(6000U);
    coarse_config.time_step_count = 8U;
    auto fine_config = coarse_config;
    fine_config.time_step_count = 32U;
    const auto coarse = pricing::price_rough_bergomi(
        model, coarse_config, {100.0}, {core::Time{1.0}},
        {pricing::VarianceReduction::standard});
    const auto fine = pricing::price_rough_bergomi(
        model, fine_config, {100.0}, {core::Time{1.0}},
        {pricing::VarianceReduction::standard});
    const double analytic = pricing::black_scholes_price(
        model.initial_asset, 100.0, core::Time{1.0}, 0.04,
        model.risk_free_discount_curve->discount_factor(core::Time{1.0}),
        model.dividend_discount_curve->discount_factor(core::Time{1.0}),
        pricing::OptionRight::call);
    for (const auto* estimate : {&coarse.estimates[0], &fine.estimates[0]}) {
        REQUIRE(std::abs(estimate->mean.value - analytic) <
                5.0 * estimate->standard_error.value + 0.03);
    }
}

TEST_CASE("Monte Carlo sampling error scales with the square root of path count") {
    auto smaller_config = config(2000U);
    auto larger_config = config(8000U);
    const auto smaller = pricing::price_rough_bergomi(
        parameters(), smaller_config, {100.0}, {core::Time{1.0}},
        {pricing::VarianceReduction::standard});
    const auto larger = pricing::price_rough_bergomi(
        parameters(), larger_config, {100.0}, {core::Time{1.0}},
        {pricing::VarianceReduction::standard});
    const double standard_error_ratio =
        larger.estimates[0].standard_error.value /
        smaller.estimates[0].standard_error.value;
    REQUIRE(standard_error_ratio > 0.35);
    REQUIRE(standard_error_ratio < 0.65);
    REQUIRE(std::abs(larger.estimates[0].mean.value -
                     smaller.estimates[0].mean.value) <
            5.0 * std::hypot(larger.estimates[0].standard_error.value,
                             smaller.estimates[0].standard_error.value));
}

TEST_CASE("Antithetic standard errors use independent pair averages") {
    const auto run = pricing::price_rough_bergomi(
        parameters(), config(), {90.0, 100.0, 110.0}, {core::Time{0.5}},
        {pricing::VarianceReduction::antithetic});
    REQUIRE(run.estimates.size() == 6U);
    for (const auto& estimate : run.estimates) {
        REQUIRE(estimate.path_count == run.path_count);
        REQUIRE(estimate.independent_observation_count == run.path_count / 2U);
    }
    REQUIRE(run.estimates[0].option.right == pricing::OptionRight::call);
    REQUIRE(run.estimates[2].option.strike == 100.0);
}

TEST_CASE("European payoff prices satisfy pathwise bounds and strike ordering") {
    const auto run = pricing::price_rough_bergomi(
        parameters(0.7), config(8000U), {90.0, 100.0, 110.0},
        {core::Time{1.0}},
        {pricing::VarianceReduction::standard});
    const double risk_free =
        std::exp(-0.03);
    const double dividend = std::exp(-0.01);
    for (std::size_t strike_index = 0; strike_index < 3U; ++strike_index) {
        const auto& call = run.estimates[strike_index * 2U];
        const auto& put = run.estimates[strike_index * 2U + 1U];
        REQUIRE(call.mean.value >=
                std::max(100.0 * dividend - call.option.strike * risk_free, 0.0) -
                    4.0 * call.standard_error.value);
        REQUIRE(put.mean.value >=
                std::max(call.option.strike * risk_free - 100.0 * dividend, 0.0) -
                    4.0 * put.standard_error.value);
        REQUIRE(call.mean.value >= 0.0);
        REQUIRE(put.mean.value >= 0.0);
    }
    REQUIRE(run.estimates[0].mean.value >= run.estimates[2].mean.value);
    REQUIRE(run.estimates[2].mean.value >= run.estimates[4].mean.value);
    REQUIRE(run.estimates[1].mean.value <= run.estimates[3].mean.value);
    REQUIRE(run.estimates[3].mean.value <= run.estimates[5].mean.value);

    const double observed_parity =
        run.estimates[2].mean.value - run.estimates[3].mean.value;
    const double expected_parity = 100.0 * dividend - 100.0 * risk_free;
    const double parity_tolerance =
        5.0 * (run.estimates[2].standard_error.value +
               run.estimates[3].standard_error.value);
    REQUIRE(std::abs(observed_parity - expected_parity) < parity_tolerance);
}

TEST_CASE("Control variate uses independent pilot paths and is replayable") {
    auto options = std::vector<pricing::VarianceReduction>{
        pricing::VarianceReduction::stock_control_variate};
    const auto first = pricing::price_rough_bergomi(
        parameters(0.7), config(4000U), {100.0}, {core::Time{0.75}}, options);
    const auto second = pricing::price_rough_bergomi(
        parameters(0.7), config(4000U), {100.0}, {core::Time{0.75}}, options);
    REQUIRE(first.pilot_path_count == 6000U);
    REQUIRE(first.estimates.size() == 2U);
    REQUIRE(std::isfinite(first.estimates[0].control_coefficient));
    REQUIRE(first.estimates[0].mean.value == second.estimates[0].mean.value);
    REQUIRE(first.estimates[0].standard_error.value ==
            second.estimates[0].standard_error.value);
}

TEST_CASE("Batching supports odd sample counts without singleton simulations") {
    auto run_config = config(3U);
    run_config.batch_size = 2U;
    const auto run = pricing::price_rough_bergomi(
        parameters(), run_config, {100.0}, {core::Time{0.25}},
        {pricing::VarianceReduction::standard});
    REQUIRE(run.estimates[0].path_count == 3U);
    REQUIRE(run.estimates[0].independent_observation_count == 3U);
}

TEST_CASE("Pricing rejects invalid configurations") {
    auto invalid = config();
    invalid.time_step_count = 0U;
    REQUIRE_THROWS_AS(
        pricing::price_rough_bergomi(
            parameters(), invalid, {100.0}, {core::Time{1.0}},
            {pricing::VarianceReduction::standard}),
        std::invalid_argument);

    invalid = config(1001U);
    REQUIRE_THROWS_AS(
        pricing::price_rough_bergomi(
            parameters(), invalid, {100.0}, {core::Time{1.0}},
            {pricing::VarianceReduction::antithetic}),
        std::invalid_argument);

    invalid = config(2U);
    invalid.batch_size = 2U;
    REQUIRE_THROWS_AS(
        pricing::price_rough_bergomi(
            parameters(), invalid, {100.0}, {core::Time{1.0}},
            {pricing::VarianceReduction::antithetic}),
        std::invalid_argument);

    REQUIRE_THROWS_AS(
        pricing::price_rough_bergomi(
            parameters(), config(), {-1.0}, {core::Time{1.0}},
            {pricing::VarianceReduction::standard}),
        std::invalid_argument);
}
