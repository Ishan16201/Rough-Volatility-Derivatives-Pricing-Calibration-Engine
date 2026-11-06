#pragma once

#include "rvre/core/types.hpp"
#include "rvre/models/rough_bergomi.hpp"
#include "rvre/pricing/option.hpp"

#include <cstddef>
#include <cstdint>
#include <functional>
#include <vector>

namespace rvre::risk {

struct SensitivityState {
    core::AssetPrice spot;
    double strike;
    core::Time maturity;
    core::Volatility volatility;
    core::InterestRate risk_free_rate;
    core::InterestRate dividend_yield;
    models::HurstExponent hurst;
    models::VolatilityOfVolatility eta;
    std::vector<core::Variance> initial_forward_variance_curve;
};

enum class RiskFactor {
    delta,
    gamma,
    vega,
    theta,
    rho,
    hurst,
    eta,
    initial_forward_variance
};

struct FiniteDifferenceConfig {
    double relative_step{1.0e-2};
    double absolute_step{1.0e-4};
    double convergence_tolerance{5.0e-2};
    std::uint32_t refinement_count{4U};
    std::size_t forward_variance_curve_index{0U};
    std::uint64_t common_random_seed{731927U};
    std::uint64_t common_random_stream_id{0U};
};

struct SensitivityEstimate {
    RiskFactor factor;
    double value;
    double finest_step;
    double refinement_disagreement;
    std::uint32_t refinements_used;
    bool converged;
};

using OptionValueFunction =
    std::function<core::OptionPrice(const SensitivityState&,
                                    std::uint64_t,
                                    std::uint64_t)>;

[[nodiscard]] std::vector<SensitivityEstimate> estimate_sensitivities(
    const SensitivityState& state,
    pricing::OptionRight right,
    const std::vector<RiskFactor>& factors,
    const FiniteDifferenceConfig& config,
    const OptionValueFunction& price);

[[nodiscard]] double black_scholes_rho(
    core::AssetPrice spot,
    double strike,
    core::Time maturity,
    core::Volatility volatility,
    core::InterestRate risk_free_rate,
    core::InterestRate dividend_yield,
    pricing::OptionRight right);

[[nodiscard]] const char* to_string(RiskFactor factor) noexcept;

}  // namespace rvre::risk
