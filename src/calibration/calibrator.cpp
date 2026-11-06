#include "rvre/calibration/calibrator.hpp"

#include "rvre/models/rough_bergomi.hpp"
#include "rvre/random/random_stream.hpp"

#include <algorithm>
#include <array>
#include <chrono>
#include <cmath>
#include <fstream>
#include <functional>
#include <iomanip>
#include <limits>
#include <map>
#include <sstream>
#include <numeric>
#include <stdexcept>
#include <string>
#include <type_traits>
#include <utility>

namespace rvre::calibration {
namespace {

constexpr std::uint64_t validation_seed_salt = 0x56414C4944415445ULL;
constexpr std::uint64_t validation_stream_salt = 0x5345434F4E445F53ULL;

struct ParameterDescriptor {
    std::string name;
    ParameterControl control;
};

struct ObjectiveRow {
    CalibrationResidual residual;
    double weighted_squared_error;
};

struct ObjectiveEvaluation {
    double value;
    double estimated_sampling_error;
    double maximum_model_price_standard_error;
    double maximum_model_price_numerical_error;
    std::vector<ObjectiveRow> rows;
};

struct ModelPriceSurface {
    std::vector<double> prices;
    std::vector<double> standard_errors;
    std::vector<double> numerical_errors;
    std::vector<std::optional<bool>> heston_quadrature_converged;
};

struct OptimizerResult {
    std::vector<double> point;
    double value;
    double initial_value;
    std::uint32_t iterations;
    CalibrationStatus status;
};

[[nodiscard]] std::string csv_number(double value) {
    std::ostringstream stream;
    stream << std::setprecision(17) << value;
    return stream.str();
}

void write_csv_record(std::ostream& output,
                      const std::vector<std::string>& fields) {
    for (std::size_t i = 0; i < fields.size(); ++i) {
        if (i != 0U) {
            output << ',';
        }
        const auto& field = fields[i];
        const bool quote =
            field.find_first_of(",\"\r\n") != std::string::npos;
        if (quote) {
            output << '"';
        }
        for (const char character : field) {
            if (character == '"') {
                output << "\"\"";
            } else {
                output << character;
            }
        }
        if (quote) {
            output << '"';
        }
    }
    output << '\n';
}

[[maybe_unused, nodiscard]] double normal_quantile(double probability) {
    constexpr double a1 = -3.969683028665376e+01;
    constexpr double a2 = 2.209460984245205e+02;
    constexpr double a3 = -2.759285104916674e+02;
    constexpr double a4 = 1.383577518672690e+02;
    constexpr double a5 = -3.066479806614716e+01;
    constexpr double a6 = 2.506628277459239e+00;
    constexpr double b1 = -5.447609879822406e+01;
    constexpr double b2 = 1.615858368580409e+02;
    constexpr double b3 = -1.556989798598866e+02;
    constexpr double b4 = 6.680131188771972e+01;
    constexpr double b5 = -1.328068155288572e+01;
    constexpr double c1 = -7.784894002430293e-03;
    constexpr double c2 = -3.223964580411365e-01;
    constexpr double c3 = -2.400758277161838e+00;
    constexpr double c4 = -2.549732539343734e+00;
    constexpr double c5 = 4.374664141464968e+00;
    constexpr double c6 = 2.938163982698783e+00;
    constexpr double d1 = 7.784695709041462e-03;
    constexpr double d2 = 3.224671290700398e-01;
    constexpr double d3 = 2.445134137142996e+00;
    constexpr double d4 = 3.754408661907416e+00;
    constexpr double tail = 0.02425;

    if (probability <= 0.0 || probability >= 1.0) {
        throw std::invalid_argument{"Confidence probability must be in (0, 1)"};
    }
    if (probability < tail) {
        const double q = std::sqrt(-2.0 * std::log(probability));
        return (((((c1 * q + c2) * q + c3) * q + c4) * q + c5) * q + c6) /
               ((((d1 * q + d2) * q + d3) * q + d4) * q + 1.0);
    }
    if (probability > 1.0 - tail) {
        const double q = std::sqrt(-2.0 * std::log1p(-probability));
        return -(((((c1 * q + c2) * q + c3) * q + c4) * q + c5) * q + c6) /
               ((((d1 * q + d2) * q + d3) * q + d4) * q + 1.0);
    }
    const double q = probability - 0.5;
    const double r = q * q;
    return (((((a1 * r + a2) * r + a3) * r + a4) * r + a5) * r + a6) * q /
           (((((b1 * r + b2) * r + b3) * r + b4) * r + b5) * r + 1.0);
}

[[nodiscard]] std::vector<ParameterDescriptor> parameter_descriptors(
    const CalibrationModel& model) {
    std::vector<ParameterDescriptor> parameters;
    std::visit(
        [&parameters](const auto& specification) {
            using ModelType = std::decay_t<decltype(specification)>;
            if constexpr (std::is_same_v<ModelType,
                                         BlackScholesCalibrationModel>) {
                parameters.push_back({"volatility", specification.volatility});
            } else if constexpr (std::is_same_v<ModelType,
                                                HestonCalibrationModel>) {
                parameters.push_back(
                    {"initial_variance", specification.initial_variance});
                parameters.push_back(
                    {"mean_reversion_speed", specification.mean_reversion_speed});
                parameters.push_back(
                    {"long_run_variance", specification.long_run_variance});
                parameters.push_back(
                    {"volatility_of_variance",
                     specification.volatility_of_variance});
                parameters.push_back({"correlation", specification.correlation});
            } else {
                parameters.push_back({"hurst", specification.hurst});
                parameters.push_back({"eta", specification.eta});
                parameters.push_back({"correlation", specification.correlation});
                for (std::size_t i = 0;
                     i < specification.forward_variance_curve.size(); ++i) {
                    parameters.push_back(
                        {"xi0_" + std::to_string(i),
                         specification.forward_variance_curve[i].variance});
                }
            }
        },
        model);
    return parameters;
}

void validate_parameter_control(const ParameterDescriptor& parameter) {
    const auto& control = parameter.control;
    if (!std::isfinite(control.initial) ||
        !std::isfinite(control.lower_bound) ||
        !std::isfinite(control.upper_bound) ||
        control.lower_bound > control.upper_bound ||
        control.initial < control.lower_bound ||
        control.initial > control.upper_bound ||
        (control.optimize && control.lower_bound == control.upper_bound)) {
        throw std::invalid_argument{"Invalid parameter bounds for " +
                                    parameter.name};
    }
}

void validate_request(const CalibrationRequest& request,
                      const std::vector<ParameterDescriptor>& descriptors) {
    validate_option_observations(request.observations);
    if (request.optimizer.maximum_iterations == 0U ||
        request.optimizer.starting_points == 0U ||
        request.optimizer.starting_points > 128U ||
        !std::isfinite(request.optimizer.objective_tolerance) ||
        request.optimizer.objective_tolerance <= 0.0 ||
        !std::isfinite(request.optimizer.parameter_tolerance) ||
        request.optimizer.parameter_tolerance <= 0.0 ||
        !std::isfinite(request.optimizer.initial_simplex_step) ||
        request.optimizer.initial_simplex_step <= 0.0 ||
        request.optimizer.initial_simplex_step > 1.0) {
        throw std::invalid_argument{"Invalid calibration optimizer settings"};
    }
    const auto& objective = request.objective;
    if (!std::isfinite(objective.minimum_reliability) ||
        objective.minimum_reliability < 0.0 ||
        objective.minimum_reliability > 1.0 ||
        !std::isfinite(objective.maximum_relative_bid_ask_width) ||
        objective.maximum_relative_bid_ask_width <= 0.0 ||
        !std::isfinite(objective.minimum_relative_bid_ask_width) ||
        objective.minimum_relative_bid_ask_width <= 0.0 ||
        objective.minimum_relative_bid_ask_width >
            objective.maximum_relative_bid_ask_width ||
        !std::isfinite(objective.maximum_inverse_width_weight) ||
        objective.maximum_inverse_width_weight <= 0.0 ||
        !std::isfinite(objective.rough_variance_curve_smoothness) ||
        objective.rough_variance_curve_smoothness < 0.0) {
        throw std::invalid_argument{"Invalid calibration objective settings"};
    }
    if (request.numerics.rough_path_count < 2U ||
        request.numerics.rough_batch_size < 2U ||
        request.numerics.rough_batch_size > request.numerics.rough_path_count ||
        request.numerics.rough_time_step_count == 0U ||
        !std::isfinite(request.numerics.heston_price_error_budget) ||
        request.numerics.heston_price_error_budget <= 0.0) {
        throw std::invalid_argument{"Invalid rough Bergomi calibration numerics"};
    }
    for (const auto& parameter : descriptors) {
        validate_parameter_control(parameter);
    }
    if (descriptors.empty()) {
        throw std::invalid_argument{"Calibration model has no parameters"};
    }

    if (const auto* model =
            std::get_if<BlackScholesCalibrationModel>(&request.model)) {
        if (model->volatility.lower_bound < 0.0) {
            throw std::invalid_argument{
                "Black-Scholes volatility bound must be nonnegative"};
        }
    } else if (const auto* heston_model =
                   std::get_if<HestonCalibrationModel>(&request.model)) {
        if (heston_model->initial_variance.lower_bound < 0.0 ||
            heston_model->mean_reversion_speed.lower_bound <= 0.0 ||
            heston_model->long_run_variance.lower_bound < 0.0 ||
            heston_model->volatility_of_variance.lower_bound < 0.0 ||
            heston_model->correlation.lower_bound < -1.0 ||
            heston_model->correlation.upper_bound > 1.0) {
            throw std::invalid_argument{
                "Heston calibration bounds violate parameter domains"};
        }
    } else {
        const auto& rough_model =
            std::get<RoughBergomiCalibrationModel>(request.model);
        if (rough_model.hurst.lower_bound <= 0.0 ||
            rough_model.hurst.upper_bound >= 0.5 ||
            rough_model.eta.lower_bound <= 0.0 ||
            rough_model.correlation.lower_bound < -1.0 ||
            rough_model.correlation.upper_bound > 1.0 ||
            rough_model.forward_variance_curve.size() < 2U ||
            rough_model.forward_variance_curve.size() > 5U) {
            throw std::invalid_argument{
                "Rough Bergomi calibration bounds or curve dimension are invalid"};
        }
        if (rough_model.forward_variance_curve.front().knot_time.value != 0.0) {
            throw std::invalid_argument{
                "Rough Bergomi initial forward variance curve must start at time zero"};
        }
        double previous_time = -1.0;
        for (const auto& knot : rough_model.forward_variance_curve) {
            if (!std::isfinite(knot.knot_time.value) ||
                knot.knot_time.value <= previous_time ||
                knot.variance.lower_bound <= 0.0) {
                throw std::invalid_argument{
                    "Rough Bergomi curve knots must be increasing with positive variance bounds"};
            }
            previous_time = knot.knot_time.value;
        }
    }

    if (request.objective.kind == ObjectiveKind::bid_ask_aware &&
        std::none_of(request.observations.begin(), request.observations.end(),
                     [](const OptionObservation& observation) {
                         return observation.bid.has_value();
                     })) {
        throw std::invalid_argument{
            "Bid-ask-aware objective requires at least one bid/ask quote"};
    }
}

[[nodiscard]] std::vector<std::size_t> free_parameter_indices(
    const std::vector<ParameterDescriptor>& descriptors) {
    std::vector<std::size_t> indices;
    for (std::size_t i = 0; i < descriptors.size(); ++i) {
        if (descriptors[i].control.optimize) {
            indices.push_back(i);
        }
    }
    return indices;
}

[[nodiscard]] std::vector<double> physical_values(
    const std::vector<ParameterDescriptor>& descriptors,
    const std::vector<std::size_t>& free_indices,
    const std::vector<double>& unit_point) {
    std::vector<double> values;
    values.reserve(descriptors.size());
    for (const auto& descriptor : descriptors) {
        values.push_back(descriptor.control.initial);
    }
    for (std::size_t i = 0; i < free_indices.size(); ++i) {
        const auto& control = descriptors[free_indices[i]].control;
        values[free_indices[i]] =
            control.lower_bound +
            unit_point[i] * (control.upper_bound - control.lower_bound);
    }
    return values;
}

[[nodiscard]] double value_for(const std::vector<ParameterDescriptor>& descriptors,
                               const std::vector<double>& values,
                               const std::string& name) {
    for (std::size_t i = 0; i < descriptors.size(); ++i) {
        if (descriptors[i].name == name) {
            return values[i];
        }
    }
    throw std::logic_error{"Missing calibration parameter: " + name};
}

[[nodiscard]] double intrinsic_value(const OptionObservation& observation) {
    const double discounted_spot =
        observation.underlying_price.value *
        observation.dividend_discount.value;
    const double discounted_strike =
        observation.strike * observation.risk_free_discount.value;
    return observation.right == pricing::OptionRight::call
               ? std::max(discounted_spot - discounted_strike, 0.0)
               : std::max(discounted_strike - discounted_spot, 0.0);
}

[[nodiscard]] ModelPriceSurface price_black_scholes_surface(
    const CalibrationRequest& request,
    const std::vector<ParameterDescriptor>& descriptors,
    const std::vector<double>& values) {
    const double volatility = value_for(descriptors, values, "volatility");
    ModelPriceSurface result;
    result.prices.reserve(request.observations.size());
    result.standard_errors.assign(request.observations.size(), 0.0);
    result.numerical_errors.assign(request.observations.size(), 0.0);
    result.heston_quadrature_converged.assign(
        request.observations.size(), std::nullopt);
    for (const auto& observation : request.observations) {
        result.prices.push_back(pricing::black_scholes_price(
            observation.underlying_price, observation.strike,
            observation.maturity,
            volatility * volatility * observation.maturity.value,
            observation.risk_free_discount, observation.dividend_discount,
            observation.right));
    }
    return result;
}

[[nodiscard]] ModelPriceSurface price_heston_surface(
    const CalibrationRequest& request,
    const std::vector<ParameterDescriptor>& descriptors,
    const std::vector<double>& values) {
    const double v0 = value_for(descriptors, values, "initial_variance");
    const double kappa = value_for(descriptors, values, "mean_reversion_speed");
    const double theta = value_for(descriptors, values, "long_run_variance");
    const double sigma = value_for(descriptors, values, "volatility_of_variance");
    const double rho = value_for(descriptors, values, "correlation");
    ModelPriceSurface surface;
    surface.prices.reserve(request.observations.size());
    surface.standard_errors.reserve(request.observations.size());
    surface.numerical_errors.reserve(request.observations.size());
    surface.heston_quadrature_converged.reserve(
        request.observations.size());
    for (const auto& observation : request.observations) {
        if (observation.maturity.value == 0.0) {
            surface.prices.push_back(intrinsic_value(observation));
            surface.standard_errors.push_back(0.0);
            surface.numerical_errors.push_back(0.0);
            surface.heston_quadrature_converged.push_back(std::nullopt);
            continue;
        }
        const double rate =
            -std::log(observation.risk_free_discount.value) /
            observation.maturity.value;
        const double dividend =
            -std::log(observation.dividend_discount.value) /
            observation.maturity.value;
        const models::HestonParameters parameters{
            .initial_asset = observation.underlying_price,
            .initial_variance = core::Variance{v0},
            .mean_reversion_speed = kappa,
            .long_run_variance = core::Variance{theta},
            .volatility_of_variance = sigma,
            .correlation = rho,
            .risk_free_rate = core::InterestRate{rate},
            .dividend_yield = core::InterestRate{dividend}};
        const auto result = models::price_heston_characteristic_function(
            parameters, observation.strike, observation.maturity,
            observation.right, request.numerics.heston_characteristic_function);
        if (result.diagnostics.estimated_quadrature_error >
                request.numerics.heston_price_error_budget) {
            throw std::runtime_error{
                "Heston characteristic-function error estimate exceeds "
                "calibration budget for " +
                observation.instrument_id};
        }
        surface.prices.push_back(result.price.value);
        surface.standard_errors.push_back(0.0);
        surface.numerical_errors.push_back(
            result.diagnostics.estimated_quadrature_error);
        surface.heston_quadrature_converged.push_back(
            result.diagnostics.converged);
    }
    return surface;
}

[[nodiscard]] std::vector<std::pair<double, double>> unique_discount_curve_points(
    const std::vector<OptionObservation>& observations,
    bool risk_free) {
    std::map<double, double> points;
    points.emplace(0.0, 1.0);
    for (const auto& observation : observations) {
        const double maturity = observation.maturity.value;
        if (maturity == 0.0) {
            continue;
        }
        const double discount =
            risk_free ? observation.risk_free_discount.value
                      : observation.dividend_discount.value;
        const auto [iterator, inserted] = points.emplace(maturity, discount);
        if (!inserted &&
            std::abs(iterator->second - discount) >
                1.0e-10 * std::max(1.0, std::abs(discount))) {
            throw std::invalid_argument{
                "Discount factors at the same maturity must be consistent"};
        }
    }
    return {points.begin(), points.end()};
}

[[nodiscard]] ModelPriceSurface price_rough_bergomi_surface(
    const CalibrationRequest& request,
    const std::vector<ParameterDescriptor>& descriptors,
    const std::vector<double>& values) {
    const auto& specification =
        std::get<RoughBergomiCalibrationModel>(request.model);
    const double first_spot = request.observations.front().underlying_price.value;
    for (const auto& observation : request.observations) {
        if (std::abs(observation.underlying_price.value - first_spot) >
            1.0e-10 * first_spot) {
            throw std::invalid_argument{
                "Rough Bergomi surface calibration requires one common underlying spot"};
        }
    }
    const double hurst = value_for(descriptors, values, "hurst");
    const double eta = value_for(descriptors, values, "eta");
    const double rho = value_for(descriptors, values, "correlation");
    std::vector<core::Time> curve_times;
    std::vector<core::Variance> curve_variances;
    curve_times.reserve(specification.forward_variance_curve.size());
    curve_variances.reserve(specification.forward_variance_curve.size());
    for (std::size_t i = 0; i < specification.forward_variance_curve.size(); ++i) {
        curve_times.push_back(specification.forward_variance_curve[i].knot_time);
        curve_variances.emplace_back(value_for(
            descriptors, values, "xi0_" + std::to_string(i)));
    }
    const auto forward_curve =
        std::make_shared<models::PiecewiseLinearForwardVarianceCurve>(
            std::move(curve_times), std::move(curve_variances));

    ModelPriceSurface surface{
        std::vector<double>(request.observations.size(), 0.0),
        std::vector<double>(request.observations.size(), 0.0),
        std::vector<double>(request.observations.size(), 0.0),
        std::vector<std::optional<bool>>(request.observations.size(),
                                         std::nullopt)};
    std::vector<core::Time> maturities;
    std::vector<double> strikes;
    for (const auto& observation : request.observations) {
        if (observation.maturity.value > 0.0) {
            maturities.push_back(observation.maturity);
            strikes.push_back(observation.strike);
        } else {
            const auto index = static_cast<std::size_t>(
                &observation - request.observations.data());
            surface.prices[index] = intrinsic_value(observation);
        }
    }
    if (maturities.empty()) {
        return surface;
    }
    const auto risk_points =
        unique_discount_curve_points(request.observations, true);
    const auto dividend_points =
        unique_discount_curve_points(request.observations, false);
    std::vector<core::Time> risk_times;
    std::vector<core::DiscountFactor> risk_factors;
    for (const auto& [time, discount] : risk_points) {
        risk_times.emplace_back(time);
        risk_factors.emplace_back(discount);
    }
    std::vector<core::Time> dividend_times;
    std::vector<core::DiscountFactor> dividend_factors;
    for (const auto& [time, discount] : dividend_points) {
        dividend_times.emplace_back(time);
        dividend_factors.emplace_back(discount);
    }
    const auto risk_curve =
        std::make_shared<models::PiecewiseLinearDiscountCurve>(
            std::move(risk_times), std::move(risk_factors));
    const auto dividend_curve =
        std::make_shared<models::PiecewiseLinearDiscountCurve>(
            std::move(dividend_times), std::move(dividend_factors));
    const models::RoughBergomiParameters parameters{
        .initial_asset = core::AssetPrice{first_spot},
        .hurst = models::HurstExponent{hurst},
        .eta = models::VolatilityOfVolatility{eta},
        .correlation = rho,
        .initial_forward_variance = forward_curve,
        .risk_free_discount_curve = risk_curve,
        .dividend_discount_curve = dividend_curve};

    std::sort(maturities.begin(), maturities.end(),
              [](core::Time lhs, core::Time rhs) {
                  return lhs.value < rhs.value;
              });
    maturities.erase(
        std::unique(maturities.begin(), maturities.end(),
                    [](core::Time lhs, core::Time rhs) {
                        return lhs.value == rhs.value;
                    }),
        maturities.end());
    std::sort(strikes.begin(), strikes.end());
    strikes.erase(std::unique(strikes.begin(), strikes.end()), strikes.end());
    const pricing::MonteCarloPricingConfig pricing_config{
        .path_count = request.numerics.rough_path_count,
        .batch_size = request.numerics.rough_batch_size,
        .time_step_count = request.numerics.rough_time_step_count,
        .pilot_path_count = 2U,
        .seed = request.numerics.common_random_seed,
        .confidence_level = 0.95,
        .model = models::RoughBergomiConfig{
            .path_count =
                static_cast<std::size_t>(request.numerics.rough_path_count),
            .seed = request.numerics.common_random_seed,
            .volterra_scheme = request.numerics.rough_volterra_scheme,
            .stream_id = request.numerics.common_random_stream_id}};
    const auto result = pricing::price_rough_bergomi(
        parameters, pricing_config, strikes, maturities,
        {pricing::VarianceReduction::standard});
    for (std::size_t i = 0; i < request.observations.size(); ++i) {
        const auto& observation = request.observations[i];
        if (observation.maturity.value == 0.0) {
            continue;
        }
        const auto maturity_iterator =
            std::lower_bound(maturities.begin(), maturities.end(),
                             observation.maturity.value,
                             [](core::Time lhs, double rhs) {
                                 return lhs.value < rhs;
                             });
        const auto strike_iterator =
            std::lower_bound(strikes.begin(), strikes.end(), observation.strike);
        const std::size_t maturity_index =
            static_cast<std::size_t>(std::distance(maturities.begin(),
                                                   maturity_iterator));
        const std::size_t strike_index =
            static_cast<std::size_t>(std::distance(strikes.begin(),
                                                   strike_iterator));
        const std::size_t rights_index =
            observation.right == pricing::OptionRight::call ? 0U : 1U;
        const std::size_t estimate_index =
            maturity_index * strikes.size() * 2U + strike_index * 2U +
            rights_index;
        if (estimate_index >= result.estimates.size()) {
            throw std::logic_error{"Rough Bergomi price surface index is inconsistent"};
        }
        surface.prices[i] = result.estimates[estimate_index].mean.value;
        surface.standard_errors[i] =
            result.estimates[estimate_index].standard_error.value;
    }
    return surface;
}

[[nodiscard]] ModelPriceSurface model_prices(
    const CalibrationRequest& request,
    const std::vector<ParameterDescriptor>& descriptors,
    const std::vector<double>& values) {
    return std::visit(
        [&request, &descriptors, &values](const auto& model) {
            using ModelType = std::decay_t<decltype(model)>;
            if constexpr (std::is_same_v<ModelType,
                                         BlackScholesCalibrationModel>) {
                return price_black_scholes_surface(request, descriptors, values);
            } else if constexpr (std::is_same_v<ModelType,
                                                HestonCalibrationModel>) {
                return price_heston_surface(request, descriptors, values);
            } else {
                return price_rough_bergomi_surface(request, descriptors, values);
            }
        },
        request.model);
}

[[nodiscard]] double market_implied_volatility(
    const OptionObservation& observation) {
    if (observation.implied_volatility) {
        return observation.implied_volatility->value;
    }
    return pricing::black_scholes_implied_volatility(
               observation.market_price.value, observation.underlying_price,
               observation.strike, observation.maturity,
               observation.risk_free_discount, observation.dividend_discount,
               observation.right)
        .volatility.value;
}

[[nodiscard]] double model_implied_volatility(
    const OptionObservation& observation,
    double model_price) {
    return pricing::black_scholes_implied_volatility(
               model_price, observation.underlying_price, observation.strike,
               observation.maturity, observation.risk_free_discount,
               observation.dividend_discount, observation.right)
        .volatility.value;
}

[[nodiscard]] double calculate_objective_residual(
    const CalibrationRequest& request,
    const OptionObservation& observation,
    double model_price,
    const std::optional<double>& market_iv,
    const std::optional<double>& model_iv) {
    if (request.objective.kind == ObjectiveKind::implied_volatility) {
        if (observation.maturity.value == 0.0) {
            return 0.0;
        }
        if (!market_iv || !model_iv) {
            throw std::logic_error{
                "Implied-volatility objective is missing a volatility value"};
        }
        return (*model_iv - *market_iv) / 0.20;
    }
    if (request.objective.kind == ObjectiveKind::bid_ask_aware &&
        observation.bid && observation.ask) {
        const double error =
            model_price < observation.bid->value
                ? model_price - observation.bid->value
                : (model_price > observation.ask->value
                       ? model_price - observation.ask->value
                       : 0.0);
        const double scale = std::max(
            0.5 * (observation.ask->value - observation.bid->value),
            request.objective.minimum_relative_bid_ask_width *
                observation.underlying_price.value);
        return error / scale;
    }
    return (model_price - observation.market_price.value) /
           observation.underlying_price.value;
}

[[nodiscard]] ObjectiveEvaluation evaluate_objective(
    const CalibrationRequest& request,
    const std::vector<ParameterDescriptor>& descriptors,
    const std::vector<double>& values) {
    const auto surface = model_prices(request, descriptors, values);
    if (surface.prices.size() != request.observations.size() ||
        surface.standard_errors.size() != request.observations.size() ||
        surface.numerical_errors.size() != request.observations.size() ||
        surface.heston_quadrature_converged.size() !=
            request.observations.size()) {
        throw std::logic_error{"Model price count differs from observation count"};
    }
    struct StagedRow {
        CalibrationResidual residual;
        double unnormalized_weight;
    };
    std::vector<StagedRow> staged;
    staged.reserve(surface.prices.size());
    double weight_sum = 0.0;
    std::size_t included_count = 0U;

    for (std::size_t i = 0; i < surface.prices.size(); ++i) {
        const auto& observation = request.observations[i];
        const double price = surface.prices[i];
        const double price_se = surface.standard_errors[i];
        const double price_error = surface.numerical_errors[i];
        if (!std::isfinite(price) || price < 0.0 ||
            !std::isfinite(price_se) || price_se < 0.0 ||
            !std::isfinite(price_error) || price_error < 0.0) {
            throw std::overflow_error{
                "Calibration model produced an invalid price for " +
                observation.instrument_id};
        }
        bool included =
            observation.reliability >= request.objective.minimum_reliability &&
            observation.reliability > 0.0;
        double relative_spread = 0.0;
        if (observation.bid && observation.ask) {
            relative_spread =
                (observation.ask->value - observation.bid->value) /
                observation.underlying_price.value;
            if (relative_spread >
                request.objective.maximum_relative_bid_ask_width) {
                included = false;
            }
        }
        if (request.objective.kind == ObjectiveKind::bid_ask_aware &&
            !observation.bid) {
            included = false;
        }
        if (request.objective.kind == ObjectiveKind::implied_volatility &&
            observation.maturity.value == 0.0) {
            included = false;
        }

        std::optional<double> market_iv;
        std::optional<double> fitted_iv;
        if (request.objective.kind == ObjectiveKind::implied_volatility &&
            observation.maturity.value > 0.0) {
            market_iv = market_implied_volatility(observation);
            fitted_iv = model_implied_volatility(observation, price);
        }
        double residual = calculate_objective_residual(
            request, observation, price, market_iv, fitted_iv);

        double weight = 1.0;
        if (request.objective.weighting == WeightingScheme::reliability) {
            weight = observation.reliability;
        } else if (request.objective.weighting ==
                   WeightingScheme::inverse_bid_ask_width) {
            if (observation.bid) {
                const double regularized_width =
                    std::max(relative_spread,
                             request.objective.minimum_relative_bid_ask_width);
                weight = std::min(1.0 / regularized_width,
                                  request.objective
                                      .maximum_inverse_width_weight);
            }
            weight *= observation.reliability;
        }
        if (!included) {
            weight = 0.0;
        } else {
            ++included_count;
        }
        weight_sum += weight;
        staged.push_back(StagedRow{
            CalibrationResidual{
                observation.instrument_id,
                observation.right,
                observation.strike,
                observation.maturity,
                observation.market_price,
                core::OptionPrice{price},
                core::OptionPrice{price_se},
                core::OptionPrice{price_error},
                surface.heston_quadrature_converged[i],
                price - observation.market_price.value,
                market_iv,
                fitted_iv,
                residual,
                weight,
                included},
            weight});
    }
    if (weight_sum <= 0.0 || included_count == 0U) {
        throw std::invalid_argument{
            "Calibration objective excluded every option observation"};
    }

    ObjectiveEvaluation result{0.0, 0.0, 0.0, 0.0, {}};
    result.rows.reserve(staged.size());
    double sampling_variance = 0.0;
    for (auto& staged_row : staged) {
        const double normalized_weight =
            staged_row.unnormalized_weight > 0.0
                ? staged_row.unnormalized_weight *
                      static_cast<double>(included_count) / weight_sum
                : 0.0;
        staged_row.residual.applied_weight = normalized_weight;
        const double contribution =
            normalized_weight * staged_row.residual.objective_residual *
            staged_row.residual.objective_residual;
        if (staged_row.residual.included) {
            result.value += contribution;
            result.maximum_model_price_standard_error =
                std::max(result.maximum_model_price_standard_error,
                         staged_row.residual.model_price_standard_error.value);
            result.maximum_model_price_numerical_error =
                std::max(result.maximum_model_price_numerical_error,
                         staged_row.residual.model_price_numerical_error.value);
            const auto& observation =
                request.observations[result.rows.size()];
            double metric_standard_error =
                staged_row.residual.model_price_standard_error.value /
                observation.underlying_price.value;
            if (request.objective.kind ==
                ObjectiveKind::implied_volatility) {
                if (staged_row.residual.model_implied_volatility &&
                    *staged_row.residual.model_implied_volatility > 0.0 &&
                    observation.maturity.value > 0.0) {
                    const double maturity = observation.maturity.value;
                    const double rate =
                        -std::log(observation.risk_free_discount.value) /
                        maturity;
                    const double dividend =
                        -std::log(observation.dividend_discount.value) /
                        maturity;
                    const auto greeks = pricing::black_scholes_greeks(
                        observation.underlying_price, observation.strike,
                        observation.maturity,
                        core::Volatility{
                            *staged_row.residual.model_implied_volatility},
                        core::InterestRate{rate},
                        core::InterestRate{dividend}, observation.right);
                    if (greeks.vega > 0.0) {
                        metric_standard_error =
                            staged_row.residual.model_price_standard_error.value /
                            (0.20 * greeks.vega);
                    } else {
                        metric_standard_error /= 0.20;
                    }
                } else {
                    metric_standard_error /= 0.20;
                }
            } else if (request.objective.kind ==
                           ObjectiveKind::bid_ask_aware &&
                       observation.bid && observation.ask) {
                const double spread_scale = std::max(
                    0.5 * (observation.ask->value - observation.bid->value),
                    request.objective.minimum_relative_bid_ask_width *
                        observation.underlying_price.value);
                metric_standard_error =
                    staged_row.residual.model_price_standard_error.value /
                    spread_scale;
            }
            const double objective_error =
                2.0 * std::abs(staged_row.residual.objective_residual) *
                    metric_standard_error +
                metric_standard_error * metric_standard_error;
            const double normalized_weight_for_mean =
                normalized_weight / static_cast<double>(included_count);
            sampling_variance +=
                normalized_weight_for_mean * normalized_weight_for_mean *
                objective_error * objective_error;
        }
        result.rows.push_back(
            ObjectiveRow{std::move(staged_row.residual), contribution});
    }
    result.value /= static_cast<double>(included_count);
    result.estimated_sampling_error = std::sqrt(sampling_variance);

    if (request.objective.rough_variance_curve_smoothness > 0.0) {
        if (std::holds_alternative<RoughBergomiCalibrationModel>(
                request.model)) {
            const auto& specification =
                std::get<RoughBergomiCalibrationModel>(request.model);
            if (specification.forward_variance_curve.size() > 2U) {
                std::vector<double> log_variances;
                std::vector<double> times;
                for (std::size_t i = 0;
                     i < specification.forward_variance_curve.size(); ++i) {
                    const double xi =
                        value_for(descriptors, values, "xi0_" +
                                                          std::to_string(i));
                    log_variances.push_back(std::log(xi));
                    times.push_back(
                        specification.forward_variance_curve[i].knot_time.value);
                }
                for (std::size_t i = 1; i + 1U < log_variances.size(); ++i) {
                    const double left_slope =
                        (log_variances[i] - log_variances[i - 1U]) /
                        (times[i] - times[i - 1U]);
                    const double right_slope =
                        (log_variances[i + 1U] - log_variances[i]) /
                        (times[i + 1U] - times[i]);
                    const double curvature = right_slope - left_slope;
                    result.value +=
                        request.objective.rough_variance_curve_smoothness *
                        curvature * curvature;
                }
            }
        }
    }
    if (!std::isfinite(result.value)) {
        throw std::overflow_error{"Calibration objective is not finite"};
    }
    return result;
}

[[nodiscard]] double clamp_unit(double value) noexcept {
    return std::clamp(value, 0.0, 1.0);
}

[[nodiscard]] std::vector<std::uint32_t> first_primes(std::size_t count) {
    std::vector<std::uint32_t> primes;
    for (std::uint32_t candidate = 2U; primes.size() < count; ++candidate) {
        bool prime = true;
        for (const std::uint32_t divisor : primes) {
            if (divisor > candidate / divisor) {
                break;
            }
            if (candidate % divisor == 0U) {
                prime = false;
                break;
            }
        }
        if (prime) {
            primes.push_back(candidate);
        }
        if (candidate == std::numeric_limits<std::uint32_t>::max()) {
            throw std::overflow_error{"Unable to generate optimizer prime bases"};
        }
    }
    return primes;
}

[[nodiscard]] double radical_inverse(std::uint32_t index,
                                     std::uint32_t base) noexcept {
    double value = 0.0;
    double factor = 1.0 / static_cast<double>(base);
    while (index > 0U) {
        value += static_cast<double>(index % base) * factor;
        index /= base;
        factor /= static_cast<double>(base);
    }
    return value;
}

[[nodiscard]] std::vector<std::vector<double>> make_starting_points(
    const std::vector<ParameterDescriptor>& descriptors,
    const std::vector<std::size_t>& free_indices,
    const CalibrationOptimizerConfig& config) {
    std::vector<std::vector<double>> points;
    if (free_indices.empty()) {
        points.emplace_back();
        return points;
    }
    std::vector<double> initial;
    for (const std::size_t index : free_indices) {
        const auto& control = descriptors[index].control;
        initial.push_back((control.initial - control.lower_bound) /
                          (control.upper_bound - control.lower_bound));
    }
    points.push_back(initial);
    const auto primes = first_primes(free_indices.size());
    random::RandomStream shifts{config.seed, 0x43414C4942524154ULL};
    std::vector<double> shift_values;
    shift_values.reserve(free_indices.size());
    for (std::size_t dimension = 0; dimension < free_indices.size();
         ++dimension) {
        shift_values.push_back(shifts.uniform_open());
    }
    for (std::uint32_t start = 1U; start < config.starting_points; ++start) {
        std::vector<double> point;
        point.reserve(free_indices.size());
        for (std::size_t dimension = 0; dimension < free_indices.size();
             ++dimension) {
            const double shifted =
                radical_inverse(start, primes[dimension]) +
                shift_values[dimension];
            point.push_back(shifted - std::floor(shifted));
        }
        points.push_back(std::move(point));
    }
    return points;
}

[[nodiscard]] OptimizerResult nelder_mead(
    const std::vector<double>& starting_point,
    const CalibrationOptimizerConfig& config,
    const std::function<double(const std::vector<double>&)>& objective) {
    const std::size_t dimensions = starting_point.size();
    if (dimensions == 0U) {
        const double value = objective({});
        return OptimizerResult{{}, value, value, 0U,
                               CalibrationStatus::converged};
    }
    struct Vertex {
        std::vector<double> point;
        double value;
    };
    std::vector<Vertex> simplex;
    simplex.reserve(dimensions + 1U);
    simplex.push_back(Vertex{starting_point, objective(starting_point)});
    for (std::size_t dimension = 0; dimension < dimensions; ++dimension) {
        auto point = starting_point;
        const double step = config.initial_simplex_step;
        point[dimension] =
            point[dimension] + step <= 1.0 ? point[dimension] + step
                                           : point[dimension] - step;
        point[dimension] = clamp_unit(point[dimension]);
        simplex.push_back(Vertex{point, objective(point)});
    }
    const double initial_value = simplex.front().value;
    std::uint32_t iteration = 0U;
    bool converged = false;
    const auto sort_simplex = [&simplex]() {
        std::sort(simplex.begin(), simplex.end(),
                  [](const Vertex& lhs, const Vertex& rhs) {
                      return lhs.value < rhs.value;
                  });
    };
    const auto clamp_point = [](std::vector<double>& point) {
        for (double& value : point) {
            value = clamp_unit(value);
        }
    };
    for (; iteration < config.maximum_iterations; ++iteration) {
        sort_simplex();
        double maximum_parameter_spread = 0.0;
        for (std::size_t vertex = 1U; vertex < simplex.size(); ++vertex) {
            for (std::size_t dimension = 0; dimension < dimensions; ++dimension) {
                maximum_parameter_spread =
                    std::max(maximum_parameter_spread,
                             std::abs(simplex[vertex].point[dimension] -
                                      simplex[0].point[dimension]));
            }
        }
        const double objective_spread =
            std::abs(simplex.back().value - simplex.front().value);
        if (maximum_parameter_spread <= config.parameter_tolerance &&
            objective_spread <= config.objective_tolerance) {
            converged = true;
            break;
        }

        std::vector<double> centroid(dimensions, 0.0);
        for (std::size_t vertex = 0; vertex < dimensions; ++vertex) {
            for (std::size_t dimension = 0; dimension < dimensions; ++dimension) {
                centroid[dimension] += simplex[vertex].point[dimension] /
                                       static_cast<double>(dimensions);
            }
        }
        const auto make_trial = [&centroid, &simplex, &clamp_point, dimensions](
                                    double factor) {
            std::vector<double> point(dimensions);
            for (std::size_t dimension = 0; dimension < dimensions; ++dimension) {
                point[dimension] =
                    centroid[dimension] +
                    factor * (centroid[dimension] -
                              simplex.back().point[dimension]);
            }
            clamp_point(point);
            return point;
        };
        auto reflected = make_trial(1.0);
        const double reflected_value = objective(reflected);
        if (reflected_value < simplex.front().value) {
            auto expanded = make_trial(2.0);
            const double expanded_value = objective(expanded);
            simplex.back() =
                expanded_value < reflected_value
                    ? Vertex{std::move(expanded), expanded_value}
                    : Vertex{std::move(reflected), reflected_value};
            continue;
        }
        if (reflected_value < simplex[dimensions - 1U].value) {
            simplex.back() = Vertex{std::move(reflected), reflected_value};
            continue;
        }

        const bool outside_contraction =
            reflected_value < simplex.back().value;
        std::vector<double> contracted(dimensions);
        for (std::size_t dimension = 0; dimension < dimensions; ++dimension) {
            contracted[dimension] =
                outside_contraction
                    ? centroid[dimension] +
                          0.5 * (reflected[dimension] - centroid[dimension])
                    : centroid[dimension] +
                          0.5 * (simplex.back().point[dimension] -
                                 centroid[dimension]);
        }
        clamp_point(contracted);
        const double contracted_value = objective(contracted);
        if (contracted_value <
            (outside_contraction ? reflected_value : simplex.back().value)) {
            simplex.back() = Vertex{std::move(contracted), contracted_value};
            continue;
        }
        for (std::size_t vertex = 1U; vertex < simplex.size(); ++vertex) {
            for (std::size_t dimension = 0; dimension < dimensions; ++dimension) {
                simplex[vertex].point[dimension] =
                    0.5 * (simplex[vertex].point[dimension] +
                           simplex.front().point[dimension]);
            }
            simplex[vertex].value = objective(simplex[vertex].point);
        }
    }
    sort_simplex();
    return OptimizerResult{
        simplex.front().point,
        simplex.front().value,
        initial_value,
        iteration,
        converged ? CalibrationStatus::converged
                  : CalibrationStatus::maximum_iterations};
}

[[nodiscard]] std::vector<ParameterEstimate> make_parameter_estimates(
    const std::vector<ParameterDescriptor>& descriptors,
    const std::vector<std::size_t>& free_indices,
    const std::vector<double>& point) {
    const auto values = physical_values(descriptors, free_indices, point);
    std::vector<ParameterEstimate> estimates;
    estimates.reserve(descriptors.size());
    for (std::size_t i = 0; i < descriptors.size(); ++i) {
        estimates.push_back(ParameterEstimate{
            descriptors[i].name,
            values[i],
            descriptors[i].control.lower_bound,
            descriptors[i].control.upper_bound,
            descriptors[i].control.optimize});
    }
    return estimates;
}

[[nodiscard]] std::vector<ParameterSensitivity> calculate_sensitivities(
    const std::vector<ParameterDescriptor>& descriptors,
    const std::vector<std::size_t>& free_indices,
    const std::vector<double>& point,
    const std::function<double(const std::vector<double>&)>& objective) {
    std::vector<ParameterSensitivity> sensitivities;
    for (std::size_t free_dimension = 0;
         free_dimension < free_indices.size(); ++free_dimension) {
        const auto& descriptor = descriptors[free_indices[free_dimension]];
        const double width =
            descriptor.control.upper_bound - descriptor.control.lower_bound;
        const double physical_value =
            descriptor.control.lower_bound +
            point[free_dimension] * width;
        const double perturbation =
            std::max(width * 1.0e-3, std::abs(physical_value) * 1.0e-4);
        auto up = point;
        auto down = point;
        up[free_dimension] =
            std::min(1.0, point[free_dimension] + perturbation / width);
        down[free_dimension] =
            std::max(0.0, point[free_dimension] - perturbation / width);
        const double up_value = objective(up);
        const double down_value = objective(down);
        const double physical_up =
            descriptor.control.lower_bound + up[free_dimension] * width;
        const double physical_down =
            descriptor.control.lower_bound + down[free_dimension] * width;
        const double derivative =
            physical_up > physical_down
                ? (up_value - down_value) / (physical_up - physical_down)
                : 0.0;
        sensitivities.push_back(ParameterSensitivity{
            descriptor.name,
            0.5 * (physical_up - physical_down),
            up_value,
            down_value,
            derivative});
    }
    return sensitivities;
}

}  // namespace

CalibrationReport calibrate(const CalibrationRequest& request) {
    const auto started = std::chrono::steady_clock::now();
    const auto descriptors = parameter_descriptors(request.model);
    validate_request(request, descriptors);
    const auto free_indices = free_parameter_indices(descriptors);
    std::uint64_t objective_evaluations = 0U;
    const auto objective = [&request, &descriptors, &free_indices,
                            &objective_evaluations](
                               const std::vector<double>& point) {
        ++objective_evaluations;
        const auto values = physical_values(descriptors, free_indices, point);
        return evaluate_objective(request, descriptors, values).value;
    };

    const auto start_points = make_starting_points(
        descriptors, free_indices, request.optimizer);
    std::vector<StartingPointResult> start_results;
    start_results.reserve(start_points.size());
    OptimizerResult best{
        {}, std::numeric_limits<double>::infinity(),
        std::numeric_limits<double>::infinity(), 0U,
        CalibrationStatus::maximum_iterations};
    for (std::size_t start_index = 0; start_index < start_points.size();
         ++start_index) {
        const auto run = nelder_mead(start_points[start_index],
                                     request.optimizer, objective);
        start_results.push_back(StartingPointResult{
            static_cast<std::uint32_t>(start_index),
            run.initial_value,
            run.value,
            run.iterations,
            run.status});
        if (run.value < best.value) {
            best = run;
        }
    }

    const auto best_values =
        physical_values(descriptors, free_indices, best.point);
    ++objective_evaluations;
    const auto final_evaluation =
        evaluate_objective(request, descriptors, best_values);
    CalibrationRequest validation_request = request;
    validation_request.numerics.common_random_seed ^= validation_seed_salt;
    validation_request.numerics.common_random_stream_id ^=
        validation_stream_salt;
    ++objective_evaluations;
    const auto validation_evaluation =
        evaluate_objective(validation_request, descriptors, best_values);
    const auto sensitivity = calculate_sensitivities(
        descriptors, free_indices, best.point, objective);
    std::vector<CalibrationResidual> residuals;
    residuals.reserve(final_evaluation.rows.size());
    for (const auto& row : final_evaluation.rows) {
        residuals.push_back(row.residual);
    }
    const double duration =
        std::chrono::duration<double>(std::chrono::steady_clock::now() -
                                      started)
            .count();

    return CalibrationReport{
        best.status,
        make_parameter_estimates(descriptors, free_indices, best.point),
        final_evaluation.value,
        best.initial_value,
        best.iterations,
        objective_evaluations,
        request.numerics.common_random_seed,
        request.numerics.common_random_stream_id,
        request.numerics.common_random_seed ^ validation_seed_salt,
        request.numerics.common_random_stream_id ^ validation_stream_salt,
        request.numerics.rough_path_count,
        request.numerics.rough_time_step_count,
        request.numerics.heston_price_error_budget,
        request.numerics.rough_volterra_scheme,
        validation_evaluation.value,
        validation_evaluation.estimated_sampling_error,
        std::max(final_evaluation.maximum_model_price_standard_error,
                 validation_evaluation.maximum_model_price_standard_error),
        std::max(final_evaluation.maximum_model_price_numerical_error,
                 validation_evaluation.maximum_model_price_numerical_error),
        duration,
        std::move(start_results),
        std::move(residuals),
        sensitivity};
}

void write_calibration_report_csv(const CalibrationReport& report,
                                  const std::string& path) {
    std::ofstream output{path};
    if (!output) {
        throw std::runtime_error{"Unable to open calibration report: " + path};
    }
    output << "record_type,status,objective,initial_objective,iterations,"
              "evaluations,duration_seconds,common_random_seed,rough_paths,"
              "rough_steps,name,value,lower_bound,upper_bound,optimized,"
              "instrument_id,right,strike,maturity,market_price,model_price,"
              "price_residual,market_iv,model_iv,objective_residual,weight,"
              "included,perturbation,objective_up,objective_down,"
              "objective_derivative,common_random_stream_id,"
              "heston_price_error_budget,rough_volterra_scheme,"
              "model_price_standard_error,model_price_numerical_error,"
              "validation_objective,estimated_objective_sampling_error,"
              "maximum_model_price_standard_error,"
              "maximum_model_price_numerical_error,"
              "validation_random_seed,validation_random_stream_id,"
              "heston_quadrature_converged\n";
    const std::string rough_scheme =
        report.rough_volterra_scheme ==
                simulation::VolterraScheme::cell_average_reference
            ? "cell_average_reference"
            : "hybrid_exact_diagonal";
    const char* status =
        report.status == CalibrationStatus::converged ? "converged"
                                                      : "maximum_iterations";
    for (const auto& parameter : report.parameters) {
        std::vector<std::string> row(43U);
        row[0] = "parameter";
        row[1] = status;
        row[2] = csv_number(report.objective_value);
        row[3] = csv_number(report.initial_objective_value);
        row[4] = std::to_string(report.iterations);
        row[5] = std::to_string(report.objective_evaluations);
        row[6] = csv_number(report.duration_seconds);
        row[7] = std::to_string(report.common_random_seed);
        row[8] = std::to_string(report.rough_path_count);
        row[9] = std::to_string(report.rough_time_step_count);
        row[10] = parameter.name;
        row[11] = csv_number(parameter.value);
        row[12] = csv_number(parameter.lower_bound);
        row[13] = csv_number(parameter.upper_bound);
        row[14] = parameter.optimized ? "1" : "0";
        row[31] = std::to_string(report.common_random_stream_id);
        row[32] = csv_number(report.heston_price_error_budget);
        row[33] = rough_scheme;
        row[36] = csv_number(report.validation_objective_value);
        row[37] = csv_number(report.estimated_objective_sampling_error);
        row[38] = csv_number(report.maximum_model_price_standard_error);
        row[39] = csv_number(report.maximum_model_price_numerical_error);
        row[40] = std::to_string(report.validation_random_seed);
        row[41] = std::to_string(report.validation_random_stream_id);
        write_csv_record(output, row);
    }
    for (const auto& residual : report.residuals) {
        std::vector<std::string> row(43U);
        row[0] = "residual";
        row[1] = status;
        row[2] = csv_number(report.objective_value);
        row[3] = csv_number(report.initial_objective_value);
        row[4] = std::to_string(report.iterations);
        row[5] = std::to_string(report.objective_evaluations);
        row[6] = csv_number(report.duration_seconds);
        row[7] = std::to_string(report.common_random_seed);
        row[8] = std::to_string(report.rough_path_count);
        row[9] = std::to_string(report.rough_time_step_count);
        row[15] = residual.instrument_id;
        row[16] =
            residual.right == pricing::OptionRight::call ? "call" : "put";
        row[17] = csv_number(residual.strike);
        row[18] = csv_number(residual.maturity.value);
        row[19] = csv_number(residual.market_price.value);
        row[20] = csv_number(residual.model_price.value);
        row[21] = csv_number(residual.price_residual);
        if (residual.market_implied_volatility) {
            row[22] = csv_number(*residual.market_implied_volatility);
        }
        if (residual.model_implied_volatility) {
            row[23] = csv_number(*residual.model_implied_volatility);
        }
        row[24] = csv_number(residual.objective_residual);
        row[25] = csv_number(residual.applied_weight);
        row[26] = residual.included ? "1" : "0";
        row[34] = csv_number(residual.model_price_standard_error.value);
        row[35] = csv_number(residual.model_price_numerical_error.value);
        if (residual.heston_quadrature_converged) {
            row[42] = *residual.heston_quadrature_converged ? "1" : "0";
        }
        row[31] = std::to_string(report.common_random_stream_id);
        row[32] = csv_number(report.heston_price_error_budget);
        row[33] = rough_scheme;
        row[36] = csv_number(report.validation_objective_value);
        row[37] = csv_number(report.estimated_objective_sampling_error);
        row[38] = csv_number(report.maximum_model_price_standard_error);
        row[39] = csv_number(report.maximum_model_price_numerical_error);
        row[40] = std::to_string(report.validation_random_seed);
        row[41] = std::to_string(report.validation_random_stream_id);
        write_csv_record(output, row);
    }
    for (const auto& sensitivity : report.sensitivities) {
        std::vector<std::string> row(43U);
        row[0] = "sensitivity";
        row[1] = status;
        row[2] = csv_number(report.objective_value);
        row[3] = csv_number(report.initial_objective_value);
        row[4] = std::to_string(report.iterations);
        row[5] = std::to_string(report.objective_evaluations);
        row[6] = csv_number(report.duration_seconds);
        row[7] = std::to_string(report.common_random_seed);
        row[8] = std::to_string(report.rough_path_count);
        row[9] = std::to_string(report.rough_time_step_count);
        row[10] = sensitivity.parameter_name;
        row[27] = csv_number(sensitivity.perturbation);
        row[28] = csv_number(sensitivity.objective_up);
        row[29] = csv_number(sensitivity.objective_down);
        row[30] = csv_number(sensitivity.central_objective_derivative);
        row[31] = std::to_string(report.common_random_stream_id);
        row[32] = csv_number(report.heston_price_error_budget);
        row[33] = rough_scheme;
        row[36] = csv_number(report.validation_objective_value);
        row[37] = csv_number(report.estimated_objective_sampling_error);
        row[38] = csv_number(report.maximum_model_price_standard_error);
        row[39] = csv_number(report.maximum_model_price_numerical_error);
        row[40] = std::to_string(report.validation_random_seed);
        row[41] = std::to_string(report.validation_random_stream_id);
        write_csv_record(output, row);
    }
    for (const auto& start : report.starting_points) {
        std::vector<std::string> row(43U);
        row[0] = "start";
        row[1] = start.status == CalibrationStatus::converged
                     ? "converged"
                     : "maximum_iterations";
        row[2] = csv_number(start.final_objective);
        row[3] = csv_number(start.initial_objective);
        row[4] = std::to_string(start.iterations);
        row[5] = std::to_string(report.objective_evaluations);
        row[6] = csv_number(report.duration_seconds);
        row[7] = std::to_string(report.common_random_seed);
        row[8] = std::to_string(report.rough_path_count);
        row[9] = std::to_string(report.rough_time_step_count);
        row[10] = "start_" + std::to_string(start.starting_point_index);
        row[31] = std::to_string(report.common_random_stream_id);
        row[32] = csv_number(report.heston_price_error_budget);
        row[33] = rough_scheme;
        row[36] = csv_number(report.validation_objective_value);
        row[37] = csv_number(report.estimated_objective_sampling_error);
        row[38] = csv_number(report.maximum_model_price_standard_error);
        row[39] = csv_number(report.maximum_model_price_numerical_error);
        row[40] = std::to_string(report.validation_random_seed);
        row[41] = std::to_string(report.validation_random_stream_id);
        write_csv_record(output, row);
    }
    output.close();
    if (!output) {
        throw std::runtime_error{"Failed writing calibration report: " + path};
    }
}

}  // namespace rvre::calibration
