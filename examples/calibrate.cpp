#include "rvre/calibration/calibrator.hpp"

#include "rvre/models/rough_bergomi.hpp"

#include <algorithm>
#include <charconv>
#include <cmath>
#include <cstdint>
#include <iostream>
#include <memory>
#include <stdexcept>
#include <string>
#include <string_view>
#include <vector>

namespace {

using namespace rvre;

struct Arguments {
    std::string model{"black_scholes"};
    std::string input_path;
    std::string output_path{"calibration_report.csv"};
    std::string objective{"price"};
    std::uint64_t seed{20261106U};
    std::uint64_t paths{1024U};
    std::uint64_t steps{32U};
    std::uint32_t starts{4U};
};

[[nodiscard]] std::uint64_t parse_integer(std::string_view text,
                                          const char* name) {
    std::uint64_t value = 0U;
    const auto [end, error] =
        std::from_chars(text.data(), text.data() + text.size(), value);
    if (error != std::errc{} || end != text.data() + text.size() || value == 0U) {
        throw std::invalid_argument{std::string{name} +
                                    " must be a positive integer"};
    }
    return value;
}

void print_help() {
    std::cout
        << "Usage: rvre_calibrate [--model black_scholes|heston|rough_bergomi]\n"
           "                      [--input observations.csv] [--output report.csv]\n"
           "                      [--objective price|implied_volatility|bid_ask]\n"
           "                      [--seed N] [--starts N] [--paths N] [--steps N]\n"
           "Without --input, a documented synthetic surface is generated.\n";
}

[[nodiscard]] Arguments parse_arguments(int argc, char* argv[]) {
    Arguments arguments;
    for (int index = 1; index < argc; ++index) {
        const std::string option{argv[index]};
        if (option == "--help") {
            print_help();
            std::exit(0);
        }
        if (index + 1 >= argc) {
            throw std::invalid_argument{"Missing value for " + option};
        }
        const std::string value{argv[++index]};
        if (option == "--model") {
            arguments.model = value;
        } else if (option == "--input") {
            arguments.input_path = value;
        } else if (option == "--output") {
            arguments.output_path = value;
        } else if (option == "--objective") {
            arguments.objective = value;
        } else if (option == "--seed") {
            arguments.seed = parse_integer(value, "seed");
        } else if (option == "--starts") {
            const auto parsed = parse_integer(value, "starts");
            if (parsed > 128U) {
                throw std::invalid_argument{"starts cannot exceed 128"};
            }
            arguments.starts = static_cast<std::uint32_t>(parsed);
        } else if (option == "--paths") {
            arguments.paths = parse_integer(value, "paths");
        } else if (option == "--steps") {
            arguments.steps = parse_integer(value, "steps");
        } else {
            throw std::invalid_argument{"Unknown option: " + option};
        }
    }
    return arguments;
}

[[nodiscard]] calibration::OptionObservation make_observation(
    double strike,
    double maturity,
    pricing::OptionRight right,
    double price) {
    return calibration::OptionObservation{
        .instrument_id = "T" + std::to_string(maturity) + "_K" +
                         std::to_string(strike) + "_" +
                         pricing::to_string(right),
        .right = right,
        .strike = strike,
        .maturity = core::Time{maturity},
        .market_price = core::OptionPrice{price},
        .bid = std::nullopt,
        .ask = std::nullopt,
        .implied_volatility = std::nullopt,
        .underlying_price = core::AssetPrice{100.0},
        .risk_free_discount =
            core::DiscountFactor{std::exp(-0.03 * maturity)},
        .dividend_discount =
            core::DiscountFactor{std::exp(-0.01 * maturity)},
        .reliability = 1.0};
}

[[nodiscard]] std::vector<calibration::OptionObservation> synthetic_surface(
    const Arguments& arguments,
    calibration::CalibrationNumerics& numerics) {
    std::vector<calibration::OptionObservation> observations;
    const std::vector<double> strikes{80.0, 90.0, 100.0, 110.0, 120.0};
    const std::vector<double> maturities{0.25, 0.5, 1.0, 1.5};

    if (arguments.model == "black_scholes") {
        constexpr double volatility = 0.27;
        for (const double maturity : maturities) {
            for (const double strike : strikes) {
                for (const auto right : {pricing::OptionRight::call,
                                         pricing::OptionRight::put}) {
                    const double price = pricing::black_scholes_price(
                        core::AssetPrice{100.0}, strike, core::Time{maturity},
                        volatility * volatility * maturity,
                        core::DiscountFactor{std::exp(-0.03 * maturity)},
                        core::DiscountFactor{std::exp(-0.01 * maturity)}, right);
                    observations.push_back(
                        make_observation(strike, maturity, right, price));
                }
            }
        }
        return observations;
    }

    if (arguments.model == "heston") {
        numerics.heston_characteristic_function.integration_upper_bound =
            400.0;
        const models::HestonParameters truth{
            .initial_asset = core::AssetPrice{100.0},
            .initial_variance = core::Variance{0.04},
            .mean_reversion_speed = 2.0,
            .long_run_variance = core::Variance{0.05},
            .volatility_of_variance = 0.5,
            .correlation = -0.7,
            .risk_free_rate = core::InterestRate{0.03},
            .dividend_yield = core::InterestRate{0.01}};
        for (const double maturity : maturities) {
            for (const double strike : strikes) {
                for (const auto right : {pricing::OptionRight::call,
                                         pricing::OptionRight::put}) {
                    const auto price =
                        models::price_heston_characteristic_function(
                            truth, strike, core::Time{maturity}, right,
                            numerics.heston_characteristic_function);
                    if (!price.diagnostics.converged) {
                        throw std::runtime_error{
                            "Synthetic Heston surface quadrature did not converge"};
                    }
                    observations.push_back(make_observation(
                        strike, maturity, right, price.price.value));
                }
            }
        }
        return observations;
    }

    if (arguments.model == "rough_bergomi") {
        constexpr double maximum_maturity = 2.0;
        const auto forward =
            std::make_shared<models::PiecewiseLinearForwardVarianceCurve>(
                std::vector<core::Time>{core::Time{0.0},
                                        core::Time{maximum_maturity}},
                std::vector<core::Variance>{core::Variance{0.04},
                                            core::Variance{0.04}});
        const auto risk =
            std::make_shared<models::PiecewiseLinearDiscountCurve>(
                std::vector<core::Time>{core::Time{0.0},
                                        core::Time{maximum_maturity}},
                std::vector<core::DiscountFactor>{
                    core::DiscountFactor{1.0},
                    core::DiscountFactor{std::exp(-0.03 * maximum_maturity)}});
        const auto dividend =
            std::make_shared<models::PiecewiseLinearDiscountCurve>(
                std::vector<core::Time>{core::Time{0.0},
                                        core::Time{maximum_maturity}},
                std::vector<core::DiscountFactor>{
                    core::DiscountFactor{1.0},
                    core::DiscountFactor{std::exp(-0.01 * maximum_maturity)}});
        const models::RoughBergomiParameters truth{
            .initial_asset = core::AssetPrice{100.0},
            .hurst = models::HurstExponent{0.1},
            .eta = models::VolatilityOfVolatility{1.4},
            .correlation = -0.7,
            .initial_forward_variance = forward,
            .risk_free_discount_curve = risk,
            .dividend_discount_curve = dividend};
        numerics.rough_path_count = arguments.paths;
        numerics.rough_batch_size = std::min<std::uint64_t>(512U, arguments.paths);
        if (numerics.rough_batch_size < 2U) {
            throw std::invalid_argument{"rough calibration needs at least two paths"};
        }
        numerics.rough_time_step_count = arguments.steps;
        numerics.common_random_seed = arguments.seed;
        pricing::MonteCarloPricingConfig config{
            .path_count = numerics.rough_path_count,
            .batch_size = numerics.rough_batch_size,
            .time_step_count = numerics.rough_time_step_count,
            .pilot_path_count = 2U,
            .seed = numerics.common_random_seed,
            .confidence_level = 0.95,
            .model = models::RoughBergomiConfig{
                .path_count = static_cast<std::size_t>(numerics.rough_path_count),
                .seed = numerics.common_random_seed,
                .volterra_scheme = numerics.rough_volterra_scheme}};
        std::vector<core::Time> time_points;
        for (const double maturity : maturities) {
            time_points.emplace_back(maturity);
        }
        const auto surface = pricing::price_rough_bergomi(
            truth, config, strikes, time_points,
            {pricing::VarianceReduction::standard});
        for (const auto& estimate : surface.estimates) {
            auto row = make_observation(
                estimate.option.strike, estimate.option.maturity.value,
                estimate.option.right, estimate.mean.value);
            row.risk_free_discount = risk->discount_factor(estimate.option.maturity);
            row.dividend_discount =
                dividend->discount_factor(estimate.option.maturity);
            observations.push_back(std::move(row));
        }
        return observations;
    }
    throw std::invalid_argument{
        "Model must be black_scholes, heston, or rough_bergomi"};
}

[[nodiscard]] calibration::CalibrationModel model_specification(
    const std::string& model) {
    using calibration::ParameterControl;
    if (model == "black_scholes") {
        return calibration::BlackScholesCalibrationModel{
            ParameterControl{0.18, 0.05, 0.8, true}};
    }
    if (model == "heston") {
        return calibration::HestonCalibrationModel{
            ParameterControl{0.035, 0.005, 0.12, true},
            ParameterControl{1.2, 0.1, 8.0, true},
            ParameterControl{0.045, 0.005, 0.12, true},
            ParameterControl{0.3, 0.05, 1.2, true},
            ParameterControl{-0.5, -0.98, 0.98, true}};
    }
    if (model == "rough_bergomi") {
        return calibration::RoughBergomiCalibrationModel{
            ParameterControl{0.15, 0.03, 0.45, true},
            ParameterControl{1.0, 0.1, 3.0, true},
            ParameterControl{-0.5, -0.98, 0.98, true},
            {{core::Time{0.0}, ParameterControl{0.035, 0.01, 0.1, true}},
             {core::Time{2.0}, ParameterControl{0.045, 0.01, 0.1, true}}}};
    }
    throw std::invalid_argument{
        "Model must be black_scholes, heston, or rough_bergomi"};
}

[[nodiscard]] calibration::ObjectiveKind parse_objective(
    const std::string& name) {
    if (name == "price") {
        return calibration::ObjectiveKind::option_price;
    }
    if (name == "implied_volatility") {
        return calibration::ObjectiveKind::implied_volatility;
    }
    if (name == "bid_ask") {
        return calibration::ObjectiveKind::bid_ask_aware;
    }
    throw std::invalid_argument{
        "Objective must be price, implied_volatility, or bid_ask"};
}

}  // namespace

int main(int argc, char* argv[]) {
    try {
        const Arguments arguments = parse_arguments(argc, argv);
        calibration::CalibrationNumerics numerics;
        std::vector<calibration::OptionObservation> observations =
            arguments.input_path.empty()
                ? synthetic_surface(arguments, numerics)
                : calibration::read_option_observations_csv(
                      arguments.input_path);
        if (arguments.input_path.empty()) {
            numerics.common_random_seed = arguments.seed;
        }
        const calibration::CalibrationRequest request{
            .model = model_specification(arguments.model),
            .observations = std::move(observations),
            .objective = calibration::ObjectiveConfig{
                .kind = parse_objective(arguments.objective),
                .weighting = calibration::WeightingScheme::reliability},
            .optimizer = calibration::CalibrationOptimizerConfig{
                .maximum_iterations = 1000U,
                .starting_points = arguments.starts,
                .objective_tolerance = 1.0e-8,
                .parameter_tolerance = 1.0e-5,
                .initial_simplex_step = 0.08,
                .seed = arguments.seed},
            .numerics = numerics};
        const auto report = calibration::calibrate(request);
        const auto heston_quadrature_warnings = static_cast<std::size_t>(
            std::count_if(
                report.residuals.begin(), report.residuals.end(),
                [](const calibration::CalibrationResidual& residual) {
                    return residual.heston_quadrature_converged.has_value() &&
                           !*residual.heston_quadrature_converged;
                }));
        calibration::write_calibration_report_csv(report, arguments.output_path);
        std::cout << "status="
                  << (report.status == calibration::CalibrationStatus::converged
                          ? "converged"
                          : "maximum_iterations")
                  << " objective=" << report.objective_value
                  << " initial_objective=" << report.initial_objective_value
                  << " iterations=" << report.iterations
                  << " evaluations=" << report.objective_evaluations
                  << " seconds=" << report.duration_seconds
                  << " observations=" << report.residuals.size()
                  << " report=" << arguments.output_path << '\n';
        if (heston_quadrature_warnings > 0U) {
            std::cerr << "Warning: Heston quadrature reported non-convergence "
                      << "for " << heston_quadrature_warnings
                      << " fitted observations; inspect the report's "
                         "quadrature flag and numerical-error columns.\n";
        }
        for (const auto& parameter : report.parameters) {
            std::cout << parameter.name << '=' << parameter.value << '\n';
        }
        return report.status == calibration::CalibrationStatus::converged ? 0 : 2;
    } catch (const std::exception& error) {
        std::cerr << "Calibration error: " << error.what() << '\n';
        return 1;
    }
}
