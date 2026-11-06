#include "rvre/models/rough_bergomi.hpp"

#include <catch2/catch_test_macros.hpp>

#include <cmath>
#include <memory>
#include <stdexcept>
#include <vector>

namespace {

using rvre::core::AssetPrice;
using rvre::core::DiscountFactor;
using rvre::core::Time;
using rvre::core::Variance;
using rvre::models::DiscountCurve;
using rvre::models::ForwardVarianceCurve;
using rvre::models::HurstExponent;
using rvre::models::PiecewiseLinearDiscountCurve;
using rvre::models::PiecewiseLinearForwardVarianceCurve;
using rvre::models::RoughBergomiConfig;
using rvre::models::RoughBergomiParameters;
using rvre::models::VolatilityOfVolatility;
using rvre::simulation::TimeGrid;
using rvre::simulation::VolterraScheme;

[[nodiscard]] std::shared_ptr<const ForwardVarianceCurve> constant_forward_variance(
    double variance = 0.04) {
    return std::make_shared<PiecewiseLinearForwardVarianceCurve>(
        std::vector<Time>{Time{0.0}, Time{2.0}},
        std::vector<Variance>{Variance{variance}, Variance{variance}});
}

[[nodiscard]] std::shared_ptr<const DiscountCurve> constant_discount_curve(
    double continuously_compounded_rate) {
    return std::make_shared<PiecewiseLinearDiscountCurve>(
        std::vector<Time>{Time{0.0}, Time{2.0}},
        std::vector<DiscountFactor>{
            DiscountFactor{1.0}, DiscountFactor{std::exp(-2.0 * continuously_compounded_rate)}});
}

[[nodiscard]] RoughBergomiParameters valid_parameters() {
    return RoughBergomiParameters{
        .initial_asset = AssetPrice{100.0},
        .hurst = HurstExponent{0.2},
        .eta = VolatilityOfVolatility{0.8},
        .correlation = -0.65,
        .initial_forward_variance = constant_forward_variance(),
        .risk_free_discount_curve = constant_discount_curve(0.03),
        .dividend_discount_curve = constant_discount_curve(0.01)};
}

[[nodiscard]] RoughBergomiConfig valid_config(std::size_t paths = 20000U,
                                               std::uint64_t seed = 198731U) {
    return RoughBergomiConfig{
        .path_count = paths,
        .seed = seed,
        .volterra_scheme = VolterraScheme::cell_average_reference};
}

class InvalidForwardVarianceCurve final : public ForwardVarianceCurve {
public:
    [[nodiscard]] Variance value(Time time) const override {
        return Variance{time.value == 0.0 ? 0.04 : -0.01};
    }
};

class InvalidDiscountCurve final : public DiscountCurve {
public:
    [[nodiscard]] DiscountFactor discount_factor(Time time) const override {
        return DiscountFactor{time.value == 0.0 ? 1.0 : 0.0};
    }
};

}  // namespace

TEST_CASE("Rough Bergomi validates model parameters and required curves") {
    auto parameters = valid_parameters();
    REQUIRE_NOTHROW(rvre::models::validate_rough_bergomi_parameters(parameters));

    parameters.hurst = HurstExponent{0.0};
    REQUIRE_THROWS_AS(rvre::models::validate_rough_bergomi_parameters(parameters),
                      std::invalid_argument);
    parameters = valid_parameters();
    parameters.hurst = HurstExponent{0.5};
    REQUIRE_THROWS_AS(rvre::models::validate_rough_bergomi_parameters(parameters),
                      std::invalid_argument);
    parameters = valid_parameters();
    parameters.eta = VolatilityOfVolatility{0.0};
    REQUIRE_THROWS_AS(rvre::models::validate_rough_bergomi_parameters(parameters),
                      std::invalid_argument);
    parameters = valid_parameters();
    parameters.correlation = -1.0;
    REQUIRE_THROWS_AS(rvre::models::validate_rough_bergomi_parameters(parameters),
                      std::invalid_argument);
    parameters = valid_parameters();
    parameters.correlation = 1.0;
    REQUIRE_THROWS_AS(rvre::models::validate_rough_bergomi_parameters(parameters),
                      std::invalid_argument);
    parameters = valid_parameters();
    parameters.initial_asset = AssetPrice{-1.0};
    REQUIRE_THROWS_AS(rvre::models::validate_rough_bergomi_parameters(parameters),
                      std::invalid_argument);
    parameters = valid_parameters();
    parameters.initial_forward_variance.reset();
    REQUIRE_THROWS_AS(rvre::models::validate_rough_bergomi_parameters(parameters),
                      std::invalid_argument);
    parameters = valid_parameters();
    parameters.risk_free_discount_curve.reset();
    REQUIRE_THROWS_AS(rvre::models::validate_rough_bergomi_parameters(parameters),
                      std::invalid_argument);
    parameters = valid_parameters();
    parameters.dividend_discount_curve.reset();
    REQUIRE_THROWS_AS(rvre::models::validate_rough_bergomi_parameters(parameters),
                      std::invalid_argument);
}

TEST_CASE("Piecewise curves enforce domains and use flat extrapolation") {
    using namespace rvre::models;

    const PiecewiseLinearForwardVarianceCurve forward_curve{
        {Time{0.0}, Time{1.0}},
        {Variance{0.04}, Variance{0.09}}};
    REQUIRE(forward_curve.value(Time{0.5}).value == 0.065);
    REQUIRE(forward_curve.value(Time{2.0}).value == 0.09);
    REQUIRE_THROWS_AS(
        (PiecewiseLinearForwardVarianceCurve{{Time{0.0}}, {Variance{0.0}}}),
        std::invalid_argument);

    REQUIRE_THROWS_AS(
        (PiecewiseLinearDiscountCurve{
            {Time{0.0}, Time{1.0}},
            {DiscountFactor{0.9}, DiscountFactor{0.8}}}),
        std::invalid_argument);
    REQUIRE_THROWS_AS(
        (PiecewiseLinearDiscountCurve{
            {Time{0.0}, Time{1.0}},
            {DiscountFactor{1.0}, DiscountFactor{0.0}}}),
        std::invalid_argument);
    REQUIRE_THROWS_AS(
        (PiecewiseLinearForwardVarianceCurve{
            {Time{0.0}, Time{0.0}},
            {Variance{0.04}, Variance{0.05}}}),
        std::invalid_argument);
}

TEST_CASE("Rough Bergomi paths are positive, reproducible, and preserve forward variance") {
    using namespace rvre;
    const auto grid = TimeGrid::uniform(Time{1.0}, 8U);
    const auto parameters = valid_parameters();
    const auto config = valid_config();
    const auto first = models::simulate_rough_bergomi(parameters, config, grid);
    const auto replay = models::simulate_rough_bergomi(parameters, config, grid);

    REQUIRE(first.asset_values() == replay.asset_values());
    REQUIRE(first.variance_values() == replay.variance_values());
    REQUIRE(first.volterra().values() == replay.volterra().values());
    REQUIRE(first.path_count() == config.path_count);
    REQUIRE(first.asset_price(0U, 0U).value == parameters.initial_asset.value);
    REQUIRE(first.variance(0U, 0U) == 0.04);

    for (std::size_t path = 0; path < first.path_count(); ++path) {
        for (std::size_t point = 0; point < grid.point_count(); ++point) {
            REQUIRE(first.asset_price(path, point).value > 0.0);
            REQUIRE(std::isfinite(first.asset_price(path, point).value));
            REQUIRE(first.variance(path, point) > 0.0);
            REQUIRE(std::isfinite(first.variance(path, point)));
        }
    }

    for (const auto& diagnostic : first.diagnostics().forward_variance_diagnostics) {
        const double difference =
            diagnostic.empirical_mean_variance.value -
            diagnostic.target_forward_variance.value;
        const double standard_error = diagnostic.mean_standard_error.value;
        REQUIRE(std::abs(difference) <= 6.0 * standard_error + 1.0e-14);
    }
}

TEST_CASE("Rough Bergomi Brownian drivers preserve the requested correlation") {
    const auto paths = rvre::models::simulate_rough_bergomi(
        valid_parameters(), valid_config(30000U, 821U),
        rvre::simulation::TimeGrid::uniform(rvre::core::Time{1.0}, 4U));
    const double rho = valid_parameters().correlation;
    const double observations =
        static_cast<double>(paths.path_count() * paths.grid().interval_count());
    REQUIRE(std::abs(paths.diagnostics().empirical_driver_correlation - rho) <
            6.0 / std::sqrt(observations));
}

TEST_CASE("Rough Bergomi discounted forward is a martingale within sampling error") {
    const auto parameters = valid_parameters();
    const auto paths = rvre::models::simulate_rough_bergomi(
        parameters,
        valid_config(50000U, 98122U),
        rvre::simulation::TimeGrid::uniform(rvre::core::Time{1.0}, 16U));
    const auto& diagnostics = paths.diagnostics();
    REQUIRE(std::abs(diagnostics.discounted_forward_mean_asset_price -
                     parameters.initial_asset.value) <
            6.0 * diagnostics.discounted_forward_standard_error);
    REQUIRE(std::abs(diagnostics.relative_discounted_forward_error) <
            6.0 * diagnostics.discounted_forward_standard_error /
                parameters.initial_asset.value);
}

TEST_CASE("Rough Bergomi supports the hybrid Volterra variance correction") {
    auto config = valid_config(20000U, 456U);
    config.volterra_scheme = VolterraScheme::hybrid_exact_diagonal;
    const auto paths = rvre::models::simulate_rough_bergomi(
        valid_parameters(), config,
        rvre::simulation::TimeGrid::uniform(rvre::core::Time{1.0}, 8U));
    REQUIRE(paths.volterra().scheme() == VolterraScheme::hybrid_exact_diagonal);
    for (const auto& diagnostic : paths.diagnostics().forward_variance_diagnostics) {
        REQUIRE(std::abs(diagnostic.empirical_mean_variance.value -
                         diagnostic.target_forward_variance.value) <
                6.0 * diagnostic.mean_standard_error.value + 1.0e-14);
    }
}

TEST_CASE("Rough Bergomi rejects invalid path counts and curves on its grid") {
    using namespace rvre;
    const auto grid = simulation::TimeGrid::uniform(core::Time{1.0}, 4U);
    REQUIRE_THROWS_AS(
        models::simulate_rough_bergomi(valid_parameters(), valid_config(1U), grid),
        std::invalid_argument);
    REQUIRE_THROWS_AS(
        models::simulate_rough_bergomi(valid_parameters(), valid_config(0U), grid),
        std::invalid_argument);

    auto parameters = valid_parameters();
    parameters.initial_forward_variance =
        std::make_shared<InvalidForwardVarianceCurve>();
    REQUIRE_THROWS_AS(
        models::simulate_rough_bergomi(parameters, valid_config(), grid),
        std::invalid_argument);

    parameters = valid_parameters();
    parameters.risk_free_discount_curve = std::make_shared<InvalidDiscountCurve>();
    REQUIRE_THROWS_AS(
        models::simulate_rough_bergomi(parameters, valid_config(), grid),
        std::invalid_argument);
}

TEST_CASE("Discrete variance correction converges under time-grid refinement") {
    using namespace rvre;
    for (const double hurst : {0.1, 0.2, 0.4}) {
        for (const auto scheme : {simulation::VolterraScheme::cell_average_reference,
                                  simulation::VolterraScheme::hybrid_exact_diagonal}) {
            const auto error = [hurst, scheme](std::size_t steps) {
                const auto grid = simulation::TimeGrid::uniform(core::Time{1.0}, steps);
                return std::abs(simulation::volterra_discrete_variance(
                                    grid, hurst, scheme, steps) -
                                1.0);
            };
            REQUIRE(error(64U) < error(16U));
            REQUIRE(error(256U) < error(64U));
        }
    }

}

TEST_CASE("Rough Bergomi discounted-forward consistency holds on refined grids") {
    using namespace rvre;
    const auto parameters = valid_parameters();
    const auto coarse = models::simulate_rough_bergomi(
        parameters,
        valid_config(30000U, 8219U),
        simulation::TimeGrid::uniform(core::Time{1.0}, 4U));
    const auto fine = models::simulate_rough_bergomi(
        parameters,
        valid_config(30000U, 18219U),
        simulation::TimeGrid::uniform(core::Time{1.0}, 32U));

    for (const auto* result : {&coarse, &fine}) {
        const auto& diagnostics = result->diagnostics();
        REQUIRE(std::abs(diagnostics.discounted_forward_mean_asset_price -
                         parameters.initial_asset.value) <
                6.0 * diagnostics.discounted_forward_standard_error);
    }
}
