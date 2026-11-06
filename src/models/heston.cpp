#include "rvre/models/heston.hpp"

#include "rvre/pricing/black_scholes.hpp"
#include "rvre/random/random_stream.hpp"

#include <algorithm>
#include <cmath>
#include <complex>
#include <cstddef>
#include <cstdint>
#include <functional>
#include <limits>
#include <numbers>
#include <stdexcept>
#include <string>

namespace rvre::models {
namespace {

using Complex = std::complex<double>;

struct RunningMoments {
    std::uint64_t count{0U};
    double mean{0.0};
    double sum_squares{0.0};

    void add(double value) {
        if (!std::isfinite(value)) {
            throw std::overflow_error{"Heston Monte Carlo observation is not finite"};
        }
        ++count;
        const double delta = value - mean;
        mean += delta / static_cast<double>(count);
        sum_squares += delta * (value - mean);
    }

    [[nodiscard]] double standard_error() const {
        if (count < 2U) {
            throw std::logic_error{"Heston Monte Carlo requires at least two observations"};
        }
        return std::sqrt(sum_squares /
                         (static_cast<double>(count - 1U) *
                          static_cast<double>(count)));
    }
};

[[nodiscard]] double normal_quantile(double probability) {
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

[[nodiscard]] Complex heston_characteristic_function(
    const HestonParameters& parameters,
    double maturity,
    Complex argument) {
    const double sigma = parameters.volatility_of_variance;
    const double kappa = parameters.mean_reversion_speed;
    const double theta = parameters.long_run_variance.value;
    const double v0 = parameters.initial_variance.value;
    const Complex iu{0.0, 1.0};
    const Complex b =
        kappa - parameters.correlation * sigma * iu * argument;
    const Complex d = std::sqrt(
        b * b + sigma * sigma * (argument * argument + iu * argument));
    const Complex denominator = b + d;
    if (std::abs(denominator) < std::numeric_limits<double>::min()) {
        throw std::runtime_error{
            "Heston characteristic function encountered a singular denominator"};
    }
    const Complex g = (b - d) / denominator;
    const Complex exp_minus_dt = std::exp(-d * maturity);
    const Complex one_minus_g_exp = 1.0 - g * exp_minus_dt;
    const Complex one_minus_g = 1.0 - g;
    const Complex log_ratio = std::log(one_minus_g_exp / one_minus_g);
    const double sigma_squared = sigma * sigma;
    const Complex c =
        (kappa * theta / sigma_squared) *
        ((b - d) * maturity - 2.0 * log_ratio);
    const Complex capital_d =
        ((b - d) / sigma_squared) *
        ((1.0 - exp_minus_dt) / one_minus_g_exp);
    const Complex log_asset{
        std::log(parameters.initial_asset.value) +
            (parameters.risk_free_rate.value -
             parameters.dividend_yield.value) *
                maturity,
        0.0};
    const Complex result = std::exp(iu * argument * log_asset + c + capital_d * v0);
    if (!std::isfinite(result.real()) || !std::isfinite(result.imag())) {
        throw std::overflow_error{
            "Heston characteristic function returned a non-finite value"};
    }
    return result;
}

struct IntegralResult {
    double value;
    double estimated_error;
    std::uint64_t evaluation_count;
    bool converged;
};

[[nodiscard]] IntegralResult adaptive_simpson(
    const std::function<double(double)>& function,
    double lower,
    double upper,
    double absolute_tolerance,
    std::uint32_t maximum_depth) {
    std::uint64_t evaluations = 0U;
    const auto evaluate = [&function, &evaluations](double x) {
        ++evaluations;
        const double value = function(x);
        if (!std::isfinite(value)) {
            throw std::overflow_error{"Heston quadrature integrand is not finite"};
        }
        return value;
    };
    const double midpoint = 0.5 * (lower + upper);
    const double f_lower = evaluate(lower);
    const double f_midpoint = evaluate(midpoint);
    const double f_upper = evaluate(upper);
    const double initial =
        (upper - lower) * (f_lower + 4.0 * f_midpoint + f_upper) / 6.0;
    bool converged = true;

    std::function<IntegralResult(double, double, double, double, double, double,
                                 double, std::uint32_t)>
        refine;
    refine = [&](double a, double b, double fa, double fm, double fb,
                 double whole, double tolerance, std::uint32_t depth) {
        const double middle = 0.5 * (a + b);
        const double left_middle = 0.5 * (a + middle);
        const double right_middle = 0.5 * (middle + b);
        const double f_left_middle = evaluate(left_middle);
        const double f_right_middle = evaluate(right_middle);
        const double left =
            (middle - a) * (fa + 4.0 * f_left_middle + fm) / 6.0;
        const double right =
            (b - middle) * (fm + 4.0 * f_right_middle + fb) / 6.0;
        const double difference = left + right - whole;
        const double error = std::abs(difference) / 15.0;
        if (error <= tolerance) {
            return IntegralResult{left + right + difference / 15.0, error, 0U,
                                  true};
        }
        if (depth == 0U) {
            converged = false;
            return IntegralResult{left + right + difference / 15.0, error, 0U,
                                  false};
        }
        const auto left_result =
            refine(a, middle, fa, f_left_middle, fm, left, tolerance * 0.5,
                   depth - 1U);
        const auto right_result =
            refine(middle, b, fm, f_right_middle, fb, right, tolerance * 0.5,
                   depth - 1U);
        return IntegralResult{
            left_result.value + right_result.value,
            left_result.estimated_error + right_result.estimated_error,
            0U,
            left_result.converged && right_result.converged};
    };
    auto result = refine(lower, upper, f_lower, f_midpoint, f_upper, initial,
                         absolute_tolerance, maximum_depth);
    result.evaluation_count = evaluations;
    result.converged = result.converged && converged;
    return result;
}

[[nodiscard]] double integrated_deterministic_variance(
    const HestonParameters& parameters,
    double maturity) {
    const double kappa = parameters.mean_reversion_speed;
    const double theta = parameters.long_run_variance.value;
    const double initial = parameters.initial_variance.value;
    return theta * maturity +
           (initial - theta) * (-std::expm1(-kappa * maturity)) / kappa;
}

[[nodiscard]] double payoff(double terminal_asset,
                            double strike,
                            pricing::OptionRight right) noexcept {
    return right == pricing::OptionRight::call
               ? std::max(terminal_asset - strike, 0.0)
               : std::max(strike - terminal_asset, 0.0);
}

}  // namespace

void validate_heston_parameters(const HestonParameters& parameters) {
    if (!std::isfinite(parameters.initial_asset.value) ||
        parameters.initial_asset.value <= 0.0 ||
        !std::isfinite(parameters.initial_variance.value) ||
        parameters.initial_variance.value < 0.0 ||
        !std::isfinite(parameters.mean_reversion_speed) ||
        parameters.mean_reversion_speed <= 0.0 ||
        !std::isfinite(parameters.long_run_variance.value) ||
        parameters.long_run_variance.value < 0.0 ||
        !std::isfinite(parameters.volatility_of_variance) ||
        parameters.volatility_of_variance < 0.0 ||
        !std::isfinite(parameters.correlation) ||
        parameters.correlation < -1.0 || parameters.correlation > 1.0 ||
        !std::isfinite(parameters.risk_free_rate.value) ||
        !std::isfinite(parameters.dividend_yield.value)) {
        throw std::invalid_argument{"Invalid Heston model parameters"};
    }
}

HestonCharacteristicFunctionResult price_heston_characteristic_function(
    const HestonParameters& parameters,
    double strike,
    core::Time maturity,
    pricing::OptionRight right,
    const HestonCharacteristicFunctionConfig& config) {
    validate_heston_parameters(parameters);
    if (!std::isfinite(strike) || strike <= 0.0 ||
        !std::isfinite(maturity.value) || maturity.value < 0.0 ||
        !std::isfinite(config.integration_upper_bound) ||
        config.integration_upper_bound <= 0.0 ||
        config.integration_upper_bound >
            std::numeric_limits<double>::max() * 0.5 ||
        !std::isfinite(config.absolute_tolerance) ||
        config.absolute_tolerance <= 0.0 ||
        config.maximum_recursion_depth == 0U ||
        config.maximum_recursion_depth > 30U ||
        (right != pricing::OptionRight::call &&
         right != pricing::OptionRight::put)) {
        throw std::invalid_argument{
            "Invalid Heston characteristic-function pricing input"};
    }

    const double discount_r =
        std::exp(-parameters.risk_free_rate.value * maturity.value);
    const double discount_q =
        std::exp(-parameters.dividend_yield.value * maturity.value);
    if (maturity.value == 0.0) {
        const double intrinsic =
            right == pricing::OptionRight::call
                ? std::max(parameters.initial_asset.value - strike, 0.0)
                : std::max(strike - parameters.initial_asset.value, 0.0);
        return HestonCharacteristicFunctionResult{
            core::OptionPrice{intrinsic},
            HestonCharacteristicFunctionDiagnostics{0.0, 0U, 0.0, true}};
    }

    if (parameters.volatility_of_variance == 0.0 ||
        (parameters.initial_variance.value == 0.0 &&
         parameters.long_run_variance.value == 0.0)) {
        const double total_variance =
            integrated_deterministic_variance(parameters, maturity.value);
        const double price = pricing::black_scholes_price(
            parameters.initial_asset, strike, maturity,
            std::max(total_variance, 0.0), core::DiscountFactor{discount_r},
            core::DiscountFactor{discount_q}, right);
        return HestonCharacteristicFunctionResult{
            core::OptionPrice{price},
            HestonCharacteristicFunctionDiagnostics{0.0, 0U, 0.0, true}};
    }

    const double log_strike = std::log(strike);
    const double epsilon = 1.0e-7;
    const auto probability_integrand = [&](double u, bool stock_measure) {
        const double safe_u = std::max(u, epsilon);
        const Complex argument{safe_u, stock_measure ? -1.0 : 0.0};
        Complex characteristic = heston_characteristic_function(
            parameters, maturity.value, argument);
        if (stock_measure) {
            const Complex normalization = heston_characteristic_function(
                parameters, maturity.value, Complex{0.0, -1.0});
            characteristic /= normalization;
        }
        const Complex rotated =
            std::exp(Complex{0.0, -safe_u * log_strike}) * characteristic;
        return rotated.imag() / safe_u;
    };
    const auto integrate_probability = [&](bool stock_measure, double upper_bound) {
        return adaptive_simpson(
            [&probability_integrand, stock_measure](double u) {
                return probability_integrand(u, stock_measure);
            },
            0.0, upper_bound, config.absolute_tolerance,
            config.maximum_recursion_depth);
    };
    const auto p1_integral = integrate_probability(
        true, config.integration_upper_bound);
    const auto p2_integral = integrate_probability(
        false, config.integration_upper_bound);
    const double doubled_upper_bound =
        2.0 * config.integration_upper_bound;
    const auto p1_doubled =
        integrate_probability(true, doubled_upper_bound);
    const auto p2_doubled =
        integrate_probability(false, doubled_upper_bound);
    const auto price_from_integrals = [&](const IntegralResult& p1_result,
                                          const IntegralResult& p2_result) {
        const double p1 =
            std::clamp(0.5 + p1_result.value / std::numbers::pi, 0.0, 1.0);
        const double p2 =
            std::clamp(0.5 + p2_result.value / std::numbers::pi, 0.0, 1.0);
        const double call =
            parameters.initial_asset.value * discount_q * p1 -
            strike * discount_r * p2;
        const double put = call - parameters.initial_asset.value * discount_q +
                           strike * discount_r;
        return std::max(right == pricing::OptionRight::call ? call : put, 0.0);
    };
    const double price = price_from_integrals(p1_integral, p2_integral);
    const double doubled_price = price_from_integrals(p1_doubled, p2_doubled);
    if (!std::isfinite(price)) {
        throw std::overflow_error{"Heston characteristic-function price is not finite"};
    }
    const double quadrature_error =
        (parameters.initial_asset.value * discount_q *
             p1_integral.estimated_error +
         strike * discount_r * p2_integral.estimated_error) /
        std::numbers::pi;
    const double cutoff_sensitivity = std::abs(doubled_price - price);
    const double low_frequency_sensitivity =
        epsilon *
        (parameters.initial_asset.value * discount_q *
             std::abs(probability_integrand(epsilon, true)) +
         strike * discount_r *
             std::abs(probability_integrand(epsilon, false))) /
        std::numbers::pi;
    const double error_estimate =
        quadrature_error + cutoff_sensitivity + low_frequency_sensitivity;
    const bool converged =
        p1_integral.converged && p2_integral.converged &&
        p1_doubled.converged && p2_doubled.converged &&
        cutoff_sensitivity <=
            std::max(config.absolute_tolerance, 1.0e-8 * std::abs(price));
    return HestonCharacteristicFunctionResult{
        core::OptionPrice{price},
        HestonCharacteristicFunctionDiagnostics{
            error_estimate,
            p1_integral.evaluation_count + p2_integral.evaluation_count +
                p1_doubled.evaluation_count + p2_doubled.evaluation_count,
            doubled_upper_bound,
            converged}};
}

HestonMonteCarloResult price_heston_monte_carlo(
    const HestonParameters& parameters,
    double strike,
    core::Time maturity,
    pricing::OptionRight right,
    const HestonMonteCarloConfig& config) {
    validate_heston_parameters(parameters);
    if (!std::isfinite(strike) || strike <= 0.0 ||
        !std::isfinite(maturity.value) || maturity.value <= 0.0 ||
        config.path_count < 2U ||
        config.path_count > std::numeric_limits<std::size_t>::max() ||
        config.time_step_count == 0U ||
        config.time_step_count > std::numeric_limits<std::size_t>::max() ||
        !std::isfinite(config.confidence_level) ||
        config.confidence_level <= 0.0 || config.confidence_level >= 1.0 ||
        config.variance_scheme != HestonVarianceScheme::full_truncation_euler ||
        (right != pricing::OptionRight::call &&
         right != pricing::OptionRight::put)) {
        throw std::invalid_argument{"Invalid Heston Monte Carlo configuration"};
    }

    const double step = maturity.value /
                       static_cast<double>(config.time_step_count);
    const double sqrt_step = std::sqrt(step);
    const double discount_r =
        std::exp(-parameters.risk_free_rate.value * maturity.value);
    const double discount_q =
        std::exp(-parameters.dividend_yield.value * maturity.value);
    const double orthogonal_scale =
        std::sqrt(std::max(0.0, 1.0 - parameters.correlation *
                                         parameters.correlation));
    random::RandomStream stream{config.seed};
    RunningMoments payoff_moments;
    RunningMoments forward_moments;
    std::uint64_t zero_updates = 0U;
    if (config.path_count >
        std::numeric_limits<std::uint64_t>::max() /
            config.time_step_count) {
        throw std::invalid_argument{"Heston total path-step count overflows uint64"};
    }
    const std::uint64_t total_updates =
        config.path_count * config.time_step_count;
    double minimum_variance = parameters.initial_variance.value;
    double maximum_variance = parameters.initial_variance.value;

    for (std::uint64_t path = 0U; path < config.path_count; ++path) {
        double log_asset = std::log(parameters.initial_asset.value);
        double variance = parameters.initial_variance.value;
        for (std::uint64_t step_index = 0U;
             step_index < config.time_step_count; ++step_index) {
            const double positive_variance = std::max(variance, 0.0);
            if (variance <= 0.0) {
                ++zero_updates;
            }
            const double variance_normal = stream.standard_normal();
            const double independent_normal = stream.standard_normal();
            const double asset_normal =
                parameters.correlation * variance_normal +
                orthogonal_scale * independent_normal;
            log_asset +=
                (parameters.risk_free_rate.value -
                 parameters.dividend_yield.value -
                 0.5 * positive_variance) *
                    step +
                std::sqrt(positive_variance) * sqrt_step * asset_normal;
            variance +=
                parameters.mean_reversion_speed *
                    (parameters.long_run_variance.value - positive_variance) *
                    step +
                parameters.volatility_of_variance *
                    std::sqrt(positive_variance) * sqrt_step *
                    variance_normal;
            if (!std::isfinite(variance) || !std::isfinite(log_asset)) {
                throw std::overflow_error{
                    "Heston full-truncation path became non-finite at path " +
                    std::to_string(path) + ", step " +
                    std::to_string(step_index)};
            }
            minimum_variance = std::min(minimum_variance, variance);
            maximum_variance = std::max(maximum_variance, variance);
        }
        const double terminal_asset = std::exp(log_asset);
        if (!std::isfinite(terminal_asset) || terminal_asset <= 0.0) {
            throw std::overflow_error{"Heston terminal asset is non-finite or non-positive"};
        }
        payoff_moments.add(discount_r * payoff(terminal_asset, strike, right));
        forward_moments.add(terminal_asset * discount_r / discount_q);
    }

    const double payoff_se = payoff_moments.standard_error();
    const double forward_se = forward_moments.standard_error();
    const double quantile = normal_quantile(0.5 + 0.5 * config.confidence_level);
    return HestonMonteCarloResult{
        core::OptionPrice{payoff_moments.mean},
        core::OptionPrice{payoff_se},
        core::OptionPrice{payoff_moments.mean - quantile * payoff_se},
        core::OptionPrice{payoff_moments.mean + quantile * payoff_se},
        config.path_count,
        config.seed,
        config.time_step_count,
        HestonMonteCarloDiagnostics{
            core::AssetPrice{forward_moments.mean},
            core::OptionPrice{forward_se},
            static_cast<double>(zero_updates) /
                static_cast<double>(total_updates),
            minimum_variance,
            maximum_variance}};
}

}  // namespace rvre::models
