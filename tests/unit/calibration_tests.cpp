#include "rvre/calibration/calibrator.hpp"

#include <catch2/catch_approx.hpp>
#include <catch2/catch_test_macros.hpp>

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <fstream>
#include <memory>
#include <optional>
#include <sstream>
#include <string>
#include <utility>
#include <vector>

namespace {

using namespace rvre;

[[nodiscard]] calibration::OptionObservation observation(
    double strike,
    double maturity,
    pricing::OptionRight right,
    double price,
    double reliability = 1.0,
    std::optional<double> bid = std::nullopt,
    std::optional<double> ask = std::nullopt) {
    return calibration::OptionObservation{
        .instrument_id = "T" + std::to_string(maturity) + "_K" +
                         std::to_string(strike),
        .right = right,
        .strike = strike,
        .maturity = core::Time{maturity},
        .market_price = core::OptionPrice{price},
        .bid = bid ? std::optional<core::OptionPrice>{core::OptionPrice{*bid}}
                   : std::nullopt,
        .ask = ask ? std::optional<core::OptionPrice>{core::OptionPrice{*ask}}
                   : std::nullopt,
        .implied_volatility = std::nullopt,
        .underlying_price = core::AssetPrice{100.0},
        .risk_free_discount =
            core::DiscountFactor{std::exp(-0.03 * maturity)},
        .dividend_discount =
            core::DiscountFactor{std::exp(-0.01 * maturity)},
        .reliability = reliability};
}

[[nodiscard]] calibration::CalibrationOptimizerConfig optimizer_config(
    std::uint32_t starts = 3U) {
    return calibration::CalibrationOptimizerConfig{
        .maximum_iterations = 300U,
        .starting_points = starts,
        .objective_tolerance = 1.0e-12,
        .parameter_tolerance = 1.0e-7,
        .initial_simplex_step = 0.1,
        .seed = 445566U};
}

[[nodiscard]] std::vector<calibration::OptionObservation>
make_black_scholes_surface(double volatility) {
    std::vector<calibration::OptionObservation> observations;
    for (const double maturity : {0.25, 0.75, 1.5}) {
        for (const double strike : {80.0, 100.0, 120.0}) {
            for (const auto right : {pricing::OptionRight::call,
                                     pricing::OptionRight::put}) {
                const double price = pricing::black_scholes_price(
                    core::AssetPrice{100.0}, strike, core::Time{maturity},
                    volatility * volatility * maturity,
                    core::DiscountFactor{std::exp(-0.03 * maturity)},
                    core::DiscountFactor{std::exp(-0.01 * maturity)}, right);
                observations.push_back(
                    observation(strike, maturity, right, price));
            }
        }
    }
    return observations;
}

[[nodiscard]] calibration::ParameterControl fixed(double value) {
    return calibration::ParameterControl{value, value, value, false};
}

}  // namespace

TEST_CASE("Synthetic Black-Scholes surface recovers volatility reproducibly") {
    constexpr double true_volatility = 0.31;
    calibration::CalibrationRequest request{
        .model = calibration::BlackScholesCalibrationModel{
            calibration::ParameterControl{0.16, 0.05, 0.8, true}},
        .observations = make_black_scholes_surface(true_volatility),
        .objective = calibration::ObjectiveConfig{
            .kind = calibration::ObjectiveKind::option_price,
            .weighting = calibration::WeightingScheme::equal},
        .optimizer = optimizer_config(),
        .numerics = {}};
    const auto report = calibration::calibrate(request);
    const auto replay = calibration::calibrate(request);
    REQUIRE(report.status == calibration::CalibrationStatus::converged);
    REQUIRE(report.parameters.size() == 1U);
    REQUIRE(report.parameters[0].name == "volatility");
    REQUIRE(report.parameters[0].value ==
            Catch::Approx(true_volatility).margin(2.0e-7));
    REQUIRE(report.objective_value < 1.0e-15);
    REQUIRE(report.validation_objective_value ==
            Catch::Approx(report.objective_value).margin(1.0e-15));
    REQUIRE(report.validation_random_seed != report.common_random_seed);
    REQUIRE(report.validation_random_stream_id !=
            report.common_random_stream_id);
    REQUIRE(report.residuals.size() == request.observations.size());
    REQUIRE(report.starting_points.size() == 3U);
    REQUIRE(report.sensitivities.size() == 1U);
    REQUIRE(replay.parameters[0].value == report.parameters[0].value);
    REQUIRE(replay.objective_value == report.objective_value);
    REQUIRE(report.objective_evaluations > 0U);
    REQUIRE(report.duration_seconds >= 0.0);
    const std::string report_path = "/tmp/rvre-calibration-report.csv";
    calibration::write_calibration_report_csv(report, report_path);
    std::ifstream report_file{report_path};
    std::string line;
    std::size_t column_count = 0U;
    REQUIRE(static_cast<bool>(std::getline(report_file, line)));
    column_count = static_cast<std::size_t>(
        std::count(line.begin(), line.end(), ',')) + 1U;
    while (std::getline(report_file, line)) {
        REQUIRE(static_cast<std::size_t>(
                    std::count(line.begin(), line.end(), ',')) + 1U ==
                column_count);
    }
    report_file.close();
    std::remove(report_path.c_str());
}

TEST_CASE("Implied-volatility objective uses observed volatilities and reports fits") {
    auto observations = make_black_scholes_surface(0.24);
    for (auto& row : observations) {
        row.implied_volatility = core::Volatility{0.24};
    }
    calibration::CalibrationRequest request{
        .model = calibration::BlackScholesCalibrationModel{
            calibration::ParameterControl{0.4, 0.05, 0.8, true}},
        .observations = std::move(observations),
        .objective = calibration::ObjectiveConfig{
            .kind = calibration::ObjectiveKind::implied_volatility,
            .weighting = calibration::WeightingScheme::equal},
        .optimizer = optimizer_config(2U),
        .numerics = {}};
    const auto report = calibration::calibrate(request);
    REQUIRE(report.parameters.front().value ==
            Catch::Approx(0.24).margin(1.0e-7));
    REQUIRE(report.residuals.front().market_implied_volatility.has_value());
    REQUIRE(*report.residuals.front().market_implied_volatility == 0.24);
    REQUIRE(report.residuals.front().model_implied_volatility.has_value());
    REQUIRE(*report.residuals.front().model_implied_volatility ==
            Catch::Approx(0.24).margin(1.0e-7));
    REQUIRE(report.objective_value < 1.0e-12);
}

TEST_CASE("Synthetic Heston surface recovers vol-of-variance with other parameters fixed") {
    constexpr double true_sigma = 0.62;
    std::vector<calibration::OptionObservation> observations;
    for (const double maturity : {0.5, 1.0, 2.0}) {
        for (const double strike : {80.0, 100.0, 120.0}) {
            for (const auto right : {pricing::OptionRight::call,
                                     pricing::OptionRight::put}) {
                const models::HestonParameters truth{
                    .initial_asset = core::AssetPrice{100.0},
                    .initial_variance = core::Variance{0.035},
                    .mean_reversion_speed = 1.7,
                    .long_run_variance = core::Variance{0.05},
                    .volatility_of_variance = true_sigma,
                    .correlation = -0.65,
                    .risk_free_rate = core::InterestRate{0.03},
                    .dividend_yield = core::InterestRate{0.01}};
                const auto price =
                    models::price_heston_characteristic_function(
                        truth, strike, core::Time{maturity}, right);
                observations.push_back(observation(
                    strike, maturity, right, price.price.value));
            }
        }
    }
    calibration::CalibrationRequest request{
        .model = calibration::HestonCalibrationModel{
            .initial_variance = fixed(0.035),
            .mean_reversion_speed = fixed(1.7),
            .long_run_variance = fixed(0.05),
            .volatility_of_variance =
                calibration::ParameterControl{0.25, 0.1, 1.2, true},
            .correlation = fixed(-0.65)},
        .observations = std::move(observations),
        .objective = calibration::ObjectiveConfig{},
        .optimizer = optimizer_config(4U),
        .numerics = {}};
    const auto report = calibration::calibrate(request);
    REQUIRE(report.status == calibration::CalibrationStatus::converged);
    const auto sigma = std::find_if(
        report.parameters.begin(), report.parameters.end(),
        [](const calibration::ParameterEstimate& parameter) {
            return parameter.name == "volatility_of_variance";
        });
    REQUIRE(sigma != report.parameters.end());
    REQUIRE(sigma->value == Catch::Approx(true_sigma).margin(2.0e-4));
    REQUIRE(report.objective_value < 1.0e-12);
    REQUIRE(report.maximum_model_price_standard_error == 0.0);
    REQUIRE(report.maximum_model_price_numerical_error > 0.0);
    REQUIRE(report.residuals.front().model_price_standard_error.value == 0.0);
    REQUIRE(report.residuals.front().model_price_numerical_error.value > 0.0);
    REQUIRE(report.residuals.front().heston_quadrature_converged.has_value());
}

TEST_CASE("Synthetic Rough Bergomi surface recovers eta with common random numbers") {
    constexpr double true_eta = 1.2;
    constexpr std::uint64_t path_count = 512U;
    constexpr std::uint64_t batch_size = 128U;
    constexpr std::uint64_t steps = 8U;
    constexpr std::uint64_t seed = 198202U;
    constexpr std::uint64_t stream_id = 718U;
    std::vector<core::Time> curve_times{core::Time{0.0}, core::Time{2.0}};
    std::vector<core::Variance> variances{core::Variance{0.04},
                                          core::Variance{0.04}};
    auto forward = std::make_shared<models::PiecewiseLinearForwardVarianceCurve>(
        curve_times, variances);
    auto risk = std::make_shared<models::PiecewiseLinearDiscountCurve>(
        std::vector<core::Time>{core::Time{0.0}, core::Time{2.0}},
        std::vector<core::DiscountFactor>{core::DiscountFactor{1.0},
                                          core::DiscountFactor{std::exp(-0.03 * 2.0)}});
    auto dividend = std::make_shared<models::PiecewiseLinearDiscountCurve>(
        std::vector<core::Time>{core::Time{0.0}, core::Time{2.0}},
        std::vector<core::DiscountFactor>{core::DiscountFactor{1.0},
                                          core::DiscountFactor{std::exp(-0.01 * 2.0)}});
    const models::RoughBergomiParameters truth{
        .initial_asset = core::AssetPrice{100.0},
        .hurst = models::HurstExponent{0.12},
        .eta = models::VolatilityOfVolatility{true_eta},
        .correlation = -0.65,
        .initial_forward_variance = forward,
        .risk_free_discount_curve = risk,
        .dividend_discount_curve = dividend};
    const pricing::MonteCarloPricingConfig price_config{
        .path_count = path_count,
        .batch_size = batch_size,
        .time_step_count = steps,
        .pilot_path_count = 2U,
        .seed = seed,
        .confidence_level = 0.95,
        .model = models::RoughBergomiConfig{
            .path_count = static_cast<std::size_t>(path_count),
            .seed = seed,
            .volterra_scheme =
                simulation::VolterraScheme::cell_average_reference,
            .stream_id = stream_id}};
    const std::vector<double> strikes{85.0, 100.0, 115.0};
    const std::vector<core::Time> maturities{core::Time{0.5},
                                              core::Time{1.0}};
    const auto surface = pricing::price_rough_bergomi(
        truth, price_config, strikes, maturities,
        {pricing::VarianceReduction::standard});
    std::vector<calibration::OptionObservation> observations;
    for (const auto& estimate : surface.estimates) {
        auto row = observation(
            estimate.option.strike, estimate.option.maturity.value,
            estimate.option.right, estimate.mean.value);
        row.risk_free_discount =
            risk->discount_factor(estimate.option.maturity);
        row.dividend_discount =
            dividend->discount_factor(estimate.option.maturity);
        observations.push_back(std::move(row));
    }
    calibration::CalibrationRequest request{
        .model = calibration::RoughBergomiCalibrationModel{
            .hurst = fixed(0.12),
            .eta = calibration::ParameterControl{0.7, 0.3, 2.2, true},
            .correlation = fixed(-0.65),
            .forward_variance_curve = {
                {core::Time{0.0}, fixed(0.04)},
                {core::Time{2.0}, fixed(0.04)}}},
        .observations = std::move(observations),
        .objective = calibration::ObjectiveConfig{},
        .optimizer = optimizer_config(2U),
        .numerics = calibration::CalibrationNumerics{
            .rough_path_count = path_count,
            .rough_batch_size = batch_size,
            .rough_time_step_count = steps,
            .common_random_seed = seed,
            .common_random_stream_id = stream_id,
            .rough_volterra_scheme =
                simulation::VolterraScheme::cell_average_reference}};
    request.optimizer.maximum_iterations = 120U;
    const auto report = calibration::calibrate(request);
    REQUIRE(report.status == calibration::CalibrationStatus::converged);
    const auto eta = std::find_if(
        report.parameters.begin(), report.parameters.end(),
        [](const calibration::ParameterEstimate& parameter) {
            return parameter.name == "eta";
        });
    REQUIRE(eta != report.parameters.end());
    REQUIRE(eta->value == Catch::Approx(true_eta).margin(0.02));
    REQUIRE(report.objective_value < 1.0e-12);
    REQUIRE(report.estimated_objective_sampling_error > 0.0);
    REQUIRE(report.maximum_model_price_standard_error > 0.0);
    REQUIRE(report.maximum_model_price_numerical_error == 0.0);
    REQUIRE(report.residuals.front().model_price_standard_error.value > 0.0);
}

TEST_CASE("Bid-ask objective filters unreliable wide quotes and reports residual rows") {
    auto observations = make_black_scholes_surface(0.25);
    for (auto& row : observations) {
        row.bid = core::OptionPrice{row.market_price.value - 0.1};
        row.ask = core::OptionPrice{row.market_price.value + 0.1};
    }
    observations.front().reliability = 0.05;
    observations.back().bid = core::OptionPrice{1.0};
    observations.back().ask = core::OptionPrice{90.0};
    calibration::CalibrationRequest request{
        .model = calibration::BlackScholesCalibrationModel{
            calibration::ParameterControl{0.3, 0.05, 0.8, true}},
        .observations = std::move(observations),
        .objective = calibration::ObjectiveConfig{
            .kind = calibration::ObjectiveKind::bid_ask_aware,
            .weighting = calibration::WeightingScheme::inverse_bid_ask_width,
            .minimum_reliability = 0.1,
            .maximum_relative_bid_ask_width = 0.5,
            .minimum_relative_bid_ask_width = 0.001,
            .maximum_inverse_width_weight = 10.0},
        .optimizer = optimizer_config(2U),
        .numerics = {}};
    const auto report = calibration::calibrate(request);
    REQUIRE(report.residuals.front().included == false);
    REQUIRE(report.residuals.back().included == false);
    REQUIRE(report.residuals[1].included);
    REQUIRE(report.residuals[1].applied_weight <=
            10.0 * static_cast<double>(report.residuals.size()));
}

TEST_CASE("Calibration CSV ingestion validates headers and quotes") {
    const std::string path = "/tmp/rvre-calibration-observations.csv";
    {
        std::ofstream output{path};
        output << "instrument_id,right,strike,maturity,market_price,bid,ask,"
                  "implied_volatility,underlying_price,risk_free_discount,"
                  "dividend_discount,reliability\n"
               << "\"quoted,id\",call,100,1,10,9.5,10.5,0.2,100,0.97,0.99,0.9\n"
               << ",put,100,1,10,,,,100,0.97,0.99,\n";
    }
    const auto rows = calibration::read_option_observations_csv(path);
    std::remove(path.c_str());
    REQUIRE(rows.size() == 2U);
    REQUIRE(rows.front().instrument_id == "quoted,id");
    REQUIRE(rows.front().bid->value == 9.5);
    REQUIRE(rows.front().reliability == 0.9);
    REQUIRE(rows.back().instrument_id.empty());
    REQUIRE_FALSE(rows.back().bid.has_value());
    REQUIRE_FALSE(rows.back().implied_volatility.has_value());
    REQUIRE(rows.back().reliability == 1.0);
}
