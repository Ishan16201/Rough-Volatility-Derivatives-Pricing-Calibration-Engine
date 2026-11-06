#include "rvre/risk/sensitivities.hpp"

#include "rvre/pricing/black_scholes.hpp"

#include <algorithm>
#include <cmath>
#include <limits>
#include <numbers>
#include <stdexcept>

namespace rvre::risk {
namespace {

[[nodiscard]] double normal_cdf(double value) {
    return 0.5 * std::erfc(-value / std::sqrt(2.0));
}

[[nodiscard]] std::uint64_t shifted_stream(std::uint64_t stream,
                                            std::size_t factor_index,
                                            unsigned int side) noexcept {
    std::uint64_t value = stream ^
                          (0x9E3779B97F4A7C15ULL *
                           (static_cast<std::uint64_t>(factor_index) + 1U));
    value ^= 0xD1B54A32D192ED03ULL * (static_cast<std::uint64_t>(side) + 1U);
    value = (value ^ (value >> 30U)) * 0xBF58476D1CE4E5B9ULL;
    value = (value ^ (value >> 27U)) * 0x94D049BB133111EBULL;
    return value ^ (value >> 31U);
}

[[nodiscard]] double& selected_value(SensitivityState& state,
                                     RiskFactor factor,
                                     std::size_t curve_index) {
    switch (factor) {
        case RiskFactor::delta:
        case RiskFactor::gamma:
            return state.spot.value;
        case RiskFactor::vega:
            return state.volatility.value;
        case RiskFactor::theta:
            return state.maturity.value;
        case RiskFactor::rho:
            return state.risk_free_rate.value;
        case RiskFactor::hurst:
            return state.hurst.value;
        case RiskFactor::eta:
            return state.eta.value;
        case RiskFactor::initial_forward_variance:
            if (curve_index >= state.initial_forward_variance_curve.size()) {
                throw std::invalid_argument{
                    "Forward-variance sensitivity index is out of range"};
            }
            return state.initial_forward_variance_curve[curve_index].value;
    }
    throw std::invalid_argument{"Unknown risk factor"};
}

[[nodiscard]] double initial_step(const SensitivityState& state,
                                  RiskFactor factor,
                                  std::size_t curve_index,
                                  const FiniteDifferenceConfig& config) {
    auto mutable_state = state;
    const double current =
        selected_value(mutable_state, factor, curve_index);
    return std::max(std::abs(current) * config.relative_step,
                    config.absolute_step);
}

void validate_config(const SensitivityState& state,
                     const std::vector<RiskFactor>& factors,
                     const FiniteDifferenceConfig& config,
                     const OptionValueFunction& price) {
    if (!price || factors.empty() || !std::isfinite(state.spot.value) ||
        state.spot.value <= 0.0 || !std::isfinite(state.strike) ||
        state.strike <= 0.0 || !std::isfinite(state.maturity.value) ||
        state.maturity.value <= 0.0 ||
        !std::isfinite(config.relative_step) || config.relative_step <= 0.0 ||
        !std::isfinite(config.absolute_step) || config.absolute_step <= 0.0 ||
        !std::isfinite(config.convergence_tolerance) ||
        config.convergence_tolerance <= 0.0 ||
        config.refinement_count < 2U || config.refinement_count > 12U ||
        !std::isfinite(state.volatility.value) ||
        !std::isfinite(state.risk_free_rate.value) ||
        !std::isfinite(state.dividend_yield.value) ||
        !std::isfinite(state.hurst.value) || !std::isfinite(state.eta.value)) {
        throw std::invalid_argument{"Invalid finite-difference sensitivity input"};
    }
}

}  // namespace

std::vector<SensitivityEstimate> estimate_sensitivities(
    const SensitivityState& state,
    pricing::OptionRight right,
    const std::vector<RiskFactor>& factors,
    const FiniteDifferenceConfig& config,
    const OptionValueFunction& price) {
    validate_config(state, factors, config, price);
    if (right != pricing::OptionRight::call &&
        right != pricing::OptionRight::put) {
        throw std::invalid_argument{"Invalid option right for sensitivity"};
    }
    std::vector<SensitivityEstimate> estimates;
    estimates.reserve(factors.size());
    for (std::size_t factor_index = 0; factor_index < factors.size();
         ++factor_index) {
        const RiskFactor factor = factors[factor_index];
        const std::size_t curve_index =
            factor == RiskFactor::initial_forward_variance
                ? config.forward_variance_curve_index
                : 0U;
        auto base_state = state;
        const double base =
            selected_value(base_state, factor, curve_index);
        const double step0 =
            initial_step(state, factor, curve_index, config);
        std::vector<double> approximations;
        approximations.reserve(config.refinement_count);
        double used_step = step0;
        for (std::uint32_t refinement = 0U;
             refinement < config.refinement_count; ++refinement) {
            const double step = std::ldexp(step0, -static_cast<int>(refinement));
            SensitivityState up = state;
            SensitivityState down = state;
            double& up_value =
                selected_value(up, factor, curve_index);
            double& down_value =
                selected_value(down, factor, curve_index);
            up_value = base + step;
            down_value = base - step;
            if (factor == RiskFactor::theta) {
                down_value = std::max(base - step, 0.0);
            } else if (factor == RiskFactor::hurst) {
                up_value = std::min(base + step, 0.499999999);
                down_value = std::max(base - step, 1.0e-8);
            } else if (factor == RiskFactor::eta ||
                       factor == RiskFactor::vega ||
                       factor == RiskFactor::initial_forward_variance) {
                down_value = std::max(base - step, 0.0);
            }
            const double actual_up = up_value - base;
            const double actual_down = base - down_value;
            if (actual_up <= 0.0 || actual_down <= 0.0) {
                throw std::invalid_argument{
                    "Finite-difference step cannot form a two-sided perturbation "
                    "inside the selected parameter domain"};
            }
            const std::uint64_t pair_stream =
                shifted_stream(config.common_random_stream_id, factor_index, 0U);
            const auto up_price =
                price(up, config.common_random_seed, pair_stream).value;
            const auto down_price =
                price(down, config.common_random_seed, pair_stream).value;
            const bool symmetric_step =
                std::abs(actual_up - actual_down) <=
                16.0 * std::numeric_limits<double>::epsilon() *
                    std::max(actual_up, actual_down);
            double estimate = 0.0;
            if (factor == RiskFactor::gamma || !symmetric_step) {
                const double base_price =
                    price(state, config.common_random_seed, pair_stream).value;
                if (factor == RiskFactor::gamma) {
                    estimate =
                        2.0 *
                        (actual_down * up_price -
                         (actual_up + actual_down) * base_price +
                         actual_up * down_price) /
                        (actual_up * actual_down *
                         (actual_up + actual_down));
                } else {
                    estimate =
                        (actual_down * actual_down * up_price +
                         (actual_up * actual_up -
                          actual_down * actual_down) *
                             base_price -
                         actual_up * actual_up * down_price) /
                        (actual_up * actual_down *
                         (actual_up + actual_down));
                }
            } else {
                estimate =
                    (up_price - down_price) /
                    (actual_up + actual_down);
            }
            if (factor == RiskFactor::theta) {
                estimate = -estimate;
            }
            if (!std::isfinite(estimate)) {
                throw std::overflow_error{
                    "Finite-difference sensitivity is not finite"};
            }
            approximations.push_back(estimate);
            used_step = std::min(actual_up, actual_down);
        }
        const double coarse = approximations[approximations.size() - 2U];
        const double fine = approximations.back();
        const double disagreement =
            std::abs(fine - coarse) /
            std::max({std::abs(fine), std::abs(coarse), 1.0e-12});
        estimates.push_back(SensitivityEstimate{
            factor, fine, used_step, disagreement, config.refinement_count,
            disagreement <= config.convergence_tolerance});
    }
    return estimates;
}

double black_scholes_rho(core::AssetPrice spot,
                          double strike,
                          core::Time maturity,
                          core::Volatility volatility,
                          core::InterestRate risk_free_rate,
                          core::InterestRate dividend_yield,
                          pricing::OptionRight right) {
    const auto greeks = pricing::black_scholes_greeks(
        spot, strike, maturity, volatility, risk_free_rate, dividend_yield,
        right);
    const double root_time = std::sqrt(maturity.value);
    const double d1 =
        (std::log(spot.value / strike) +
         (risk_free_rate.value - dividend_yield.value +
          0.5 * volatility.value * volatility.value) *
             maturity.value) /
        (volatility.value * root_time);
    const double d2 = d1 - volatility.value * root_time;
    const double discounted_strike =
        strike * std::exp(-risk_free_rate.value * maturity.value);
    const double rho =
        right == pricing::OptionRight::call
            ? maturity.value * discounted_strike * normal_cdf(d2)
            : -maturity.value * discounted_strike * normal_cdf(-d2);
    (void)greeks;
    return rho;
}

const char* to_string(RiskFactor factor) noexcept {
    switch (factor) {
        case RiskFactor::delta: return "delta";
        case RiskFactor::gamma: return "gamma";
        case RiskFactor::vega: return "vega";
        case RiskFactor::theta: return "theta";
        case RiskFactor::rho: return "rho";
        case RiskFactor::hurst: return "hurst";
        case RiskFactor::eta: return "eta";
        case RiskFactor::initial_forward_variance:
            return "initial_forward_variance";
    }
    return "unknown";
}

}  // namespace rvre::risk
