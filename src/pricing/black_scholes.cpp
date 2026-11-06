#include "rvre/pricing/black_scholes.hpp"

#include "rvre/random/random_stream.hpp"

#include <algorithm>
#include <cmath>
#include <limits>
#include <numbers>
#include <stdexcept>

namespace rvre::pricing {
namespace {

[[nodiscard]] double normal_cdf(double value) {
    return 0.5 * std::erfc(-value / std::sqrt(2.0));
}

[[nodiscard]] double normal_pdf(double value) {
    return std::exp(-0.5 * value * value) /
           std::sqrt(2.0 * std::numbers::pi);
}

void validate_option_inputs(core::AssetPrice initial_asset,
                            double strike,
                            core::Time maturity,
                            OptionRight right) {
    if (!std::isfinite(initial_asset.value) || initial_asset.value <= 0.0 ||
        !std::isfinite(strike) || strike <= 0.0 ||
        !std::isfinite(maturity.value) || maturity.value < 0.0 ||
        (right != OptionRight::call && right != OptionRight::put)) {
        throw std::invalid_argument{"Invalid Black-Scholes option input"};
    }
}

struct RunningMoments {
    std::uint64_t count{0U};
    double mean{0.0};
    double sum_squares{0.0};

    void add(double value) {
        ++count;
        const double delta = value - mean;
        mean += delta / static_cast<double>(count);
        sum_squares += delta * (value - mean);
    }
};

[[nodiscard]] double normal_quantile(double probability) {
    constexpr double a1 = -3.969683028665376e+01;
    constexpr double a2 = 2.209460984245205e+02;
    constexpr double a3 = -2.759285104469687e+02;
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

}  // namespace

double black_scholes_price(core::AssetPrice initial_asset,
                           double strike,
                           core::Time maturity,
                           double total_variance,
                           core::DiscountFactor risk_free_discount,
                           core::DiscountFactor dividend_discount,
                           OptionRight right) {
    validate_option_inputs(initial_asset, strike, maturity, right);
    if (!std::isfinite(total_variance) || total_variance < 0.0 ||
        !std::isfinite(risk_free_discount.value) ||
        risk_free_discount.value <= 0.0 ||
        !std::isfinite(dividend_discount.value) ||
        dividend_discount.value <= 0.0) {
        throw std::invalid_argument{"Invalid Black-Scholes variance or discount factor"};
    }
    if (maturity.value == 0.0 &&
        (std::abs(risk_free_discount.value - 1.0) > 1.0e-12 ||
         std::abs(dividend_discount.value - 1.0) > 1.0e-12)) {
        throw std::invalid_argument{
            "Zero-maturity Black-Scholes discount factors must equal one"};
    }

    const double discounted_asset = initial_asset.value * dividend_discount.value;
    const double discounted_strike = strike * risk_free_discount.value;
    if (maturity.value == 0.0 || total_variance == 0.0) {
        return right == OptionRight::call
                   ? std::max(discounted_asset - discounted_strike, 0.0)
                   : std::max(discounted_strike - discounted_asset, 0.0);
    }
    const double root_variance = std::sqrt(total_variance);
    const double d1 = (std::log(discounted_asset / discounted_strike) +
                       0.5 * total_variance) /
                      root_variance;
    const double d2 = d1 - root_variance;
    const double value =
        right == OptionRight::call
            ? discounted_asset * normal_cdf(d1) -
                  discounted_strike * normal_cdf(d2)
            : discounted_strike * normal_cdf(-d2) -
                  discounted_asset * normal_cdf(-d1);
    if (!std::isfinite(value)) {
        throw std::overflow_error{"Black-Scholes price is not finite"};
    }
    return std::max(value, 0.0);
}

BlackScholesGreeks black_scholes_greeks(
    core::AssetPrice initial_asset,
    double strike,
    core::Time maturity,
    core::Volatility volatility,
    core::InterestRate risk_free_rate,
    core::InterestRate dividend_yield,
    OptionRight right) {
    validate_option_inputs(initial_asset, strike, maturity, right);
    if (!std::isfinite(volatility.value) || volatility.value <= 0.0 ||
        !std::isfinite(risk_free_rate.value) ||
        !std::isfinite(dividend_yield.value) || maturity.value <= 0.0) {
        throw std::invalid_argument{
            "Black-Scholes Greeks require positive volatility and maturity and finite rates"};
    }
    const double root_time = std::sqrt(maturity.value);
    const double sigma_root_time = volatility.value * root_time;
    const double d1 =
        (std::log(initial_asset.value / strike) +
         (risk_free_rate.value - dividend_yield.value +
          0.5 * volatility.value * volatility.value) *
             maturity.value) /
        sigma_root_time;
    const double d2 = d1 - sigma_root_time;
    const double dividend_discount =
        std::exp(-dividend_yield.value * maturity.value);
    const double risk_free_discount =
        std::exp(-risk_free_rate.value * maturity.value);
    const double density = normal_pdf(d1);
    const double delta =
        right == OptionRight::call
            ? dividend_discount * normal_cdf(d1)
            : dividend_discount * (normal_cdf(d1) - 1.0);
    const double gamma =
        dividend_discount * density /
        (initial_asset.value * sigma_root_time);
    const double vega = initial_asset.value * dividend_discount *
                        density * root_time;
    const double common_theta =
        -initial_asset.value * dividend_discount * density *
        volatility.value / (2.0 * root_time);
    const double theta =
        right == OptionRight::call
            ? common_theta -
                  risk_free_rate.value * strike * risk_free_discount *
                      normal_cdf(d2) +
                  dividend_yield.value * initial_asset.value *
                      dividend_discount * normal_cdf(d1)
            : common_theta +
                  risk_free_rate.value * strike * risk_free_discount *
                      normal_cdf(-d2) -
                  dividend_yield.value * initial_asset.value *
                      dividend_discount * normal_cdf(-d1);
    if (!std::isfinite(delta) || !std::isfinite(gamma) ||
        !std::isfinite(vega) || !std::isfinite(theta)) {
        throw std::overflow_error{"Black-Scholes Greek is not finite"};
    }
    return BlackScholesGreeks{delta, gamma, vega, theta};
}

ImpliedVolatilityResult black_scholes_implied_volatility(
    double option_price,
    core::AssetPrice initial_asset,
    double strike,
    core::Time maturity,
    core::DiscountFactor risk_free_discount,
    core::DiscountFactor dividend_discount,
    OptionRight right,
    double price_tolerance,
    std::uint32_t maximum_iterations) {
    validate_option_inputs(initial_asset, strike, maturity, right);
    if (!std::isfinite(option_price) || !std::isfinite(price_tolerance) ||
        price_tolerance <= 0.0 || maximum_iterations == 0U ||
        !std::isfinite(risk_free_discount.value) ||
        risk_free_discount.value <= 0.0 ||
        !std::isfinite(dividend_discount.value) ||
        dividend_discount.value <= 0.0) {
        throw std::invalid_argument{"Invalid implied-volatility input"};
    }
    if (maturity.value == 0.0 &&
        (std::abs(risk_free_discount.value - 1.0) > 1.0e-12 ||
         std::abs(dividend_discount.value - 1.0) > 1.0e-12)) {
        throw std::invalid_argument{
            "Zero-maturity implied-volatility discount factors must equal one"};
    }
    const double discounted_asset = initial_asset.value * dividend_discount.value;
    const double discounted_strike = strike * risk_free_discount.value;
    const double lower_bound =
        right == OptionRight::call
            ? std::max(discounted_asset - discounted_strike, 0.0)
            : std::max(discounted_strike - discounted_asset, 0.0);
    const double upper_bound =
        right == OptionRight::call ? discounted_asset : discounted_strike;
    if (option_price < lower_bound - price_tolerance ||
        option_price > upper_bound + price_tolerance) {
        throw std::invalid_argument{
            "Option price violates Black-Scholes no-arbitrage bounds"};
    }
    if (maturity.value == 0.0) {
        if (std::abs(option_price - lower_bound) <= price_tolerance) {
            return ImpliedVolatilityResult{core::Volatility{0.0}, 0U,
                                           lower_bound - option_price};
        }
        throw std::invalid_argument{
            "Implied volatility is undefined at zero maturity"};
    }
    if (option_price <= lower_bound + price_tolerance) {
        return ImpliedVolatilityResult{
            core::Volatility{0.0}, 0U, lower_bound - option_price};
    }
    if (option_price >= upper_bound - price_tolerance) {
        throw std::invalid_argument{
            "The upper no-arbitrage price bound has no finite implied volatility"};
    }

    const auto price_residual = [&](double volatility_value) {
        return black_scholes_price(
                   initial_asset, strike, maturity,
                   volatility_value * volatility_value * maturity.value,
                   risk_free_discount, dividend_discount, right) -
               option_price;
    };
    double a = 0.0;
    double b = 1.0;
    double fa = price_residual(a);
    double fb = price_residual(b);
    while (fb < 0.0 && b < 16.0) {
        b *= 2.0;
        fb = price_residual(b);
    }
    if (fb < 0.0) {
        throw std::runtime_error{
            "Unable to bracket implied volatility below volatility 16"};
    }

    double c = a;
    double fc = fa;
    double d = b - a;
    double e = d;
    double root = b;
    std::uint32_t iteration = 0U;
    bool converged = false;
    for (; iteration < maximum_iterations; ++iteration) {
        if ((fb > 0.0 && fc > 0.0) || (fb < 0.0 && fc < 0.0)) {
            c = a;
            fc = fa;
            d = b - a;
            e = d;
        }
        if (std::abs(fc) < std::abs(fb)) {
            const double previous_b = b;
            const double previous_fb = fb;
            a = b;
            fa = fb;
            b = c;
            fb = fc;
            c = previous_b;
            fc = previous_fb;
        }
        const double x_tolerance =
            2.0 * std::numeric_limits<double>::epsilon() * std::abs(b) +
            1.0e-12;
        const double midpoint = 0.5 * (c - b);
        if (std::abs(fb) <= price_tolerance ||
            std::abs(midpoint) <= x_tolerance) {
            root = b;
            converged = true;
            break;
        }

        if (std::abs(e) >= x_tolerance && std::abs(fa) > std::abs(fb)) {
            const double s = fb / fa;
            double p = 0.0;
            double q = 0.0;
            if (a == c) {
                p = 2.0 * midpoint * s;
                q = 1.0 - s;
            } else {
                const double r = fb / fc;
                const double ratio = fa / fc;
                p = s * (2.0 * midpoint * ratio *
                             (ratio - r) -
                         (b - a) * (r - 1.0));
                q = (ratio - 1.0) * (r - 1.0) * (s - 1.0);
            }
            if (p > 0.0) {
                q = -q;
            } else {
                p = -p;
            }
            const double previous_e = e;
            e = d;
            if (q != 0.0 &&
                2.0 * p <
                    std::min(3.0 * midpoint * q -
                                 std::abs(x_tolerance * q),
                             std::abs(previous_e * q))) {
                d = p / q;
            } else {
                d = midpoint;
                e = midpoint;
            }
        } else {
            d = midpoint;
            e = midpoint;
        }
        a = b;
        fa = fb;
        b += std::abs(d) > x_tolerance
                 ? d
                 : std::copysign(x_tolerance, midpoint);
        fb = price_residual(b);
        root = b;
    }
    if (!converged) {
        const double final_residual = price_residual(root);
        if (std::abs(final_residual) > price_tolerance) {
            throw std::runtime_error{"Brent implied-volatility solver did not converge"};
        }
    }
    const double model_price =
        option_price + price_residual(root);
    return ImpliedVolatilityResult{
        core::Volatility{root},
        std::min(iteration + 1U, maximum_iterations),
        model_price - option_price};
}

BlackScholesMonteCarloResult price_black_scholes_monte_carlo(
    const core::BlackScholesParameters& parameters,
    double strike,
    OptionRight right,
    std::uint64_t path_count,
    std::uint64_t seed,
    double confidence_level) {
    validate_option_inputs(parameters.initial_asset, strike,
                           parameters.maturity, right);
    if (path_count < 2U || path_count > std::numeric_limits<std::size_t>::max() ||
        !std::isfinite(parameters.volatility.value) ||
        parameters.volatility.value < 0.0 ||
        !std::isfinite(parameters.risk_free_rate.value) ||
        !std::isfinite(parameters.dividend_yield.value) ||
        !std::isfinite(confidence_level) || confidence_level <= 0.0 ||
        confidence_level >= 1.0) {
        throw std::invalid_argument{"Invalid Black-Scholes Monte Carlo configuration"};
    }
    const double maturity = parameters.maturity.value;
    const double risk_free_discount =
        std::exp(-parameters.risk_free_rate.value * maturity);
    const double drift =
        (parameters.risk_free_rate.value - parameters.dividend_yield.value -
         0.5 * parameters.volatility.value * parameters.volatility.value) *
        maturity;
    const double diffusion =
        parameters.volatility.value * std::sqrt(maturity);
    if (!std::isfinite(risk_free_discount) || !std::isfinite(drift) ||
        !std::isfinite(diffusion)) {
        throw std::overflow_error{
            "Black-Scholes Monte Carlo inputs overflowed during discounting"};
    }
    random::RandomStream stream{seed};
    RunningMoments moments;
    for (std::uint64_t path = 0U; path < path_count; ++path) {
        const double terminal =
            parameters.initial_asset.value *
            std::exp(drift + diffusion * stream.standard_normal());
        const double payoff =
            right == OptionRight::call
                ? std::max(terminal - strike, 0.0)
                : std::max(strike - terminal, 0.0);
        const double discounted_payoff = risk_free_discount * payoff;
        if (!std::isfinite(discounted_payoff)) {
            throw std::overflow_error{"Black-Scholes Monte Carlo payoff is not finite"};
        }
        moments.add(discounted_payoff);
    }
    const double standard_error =
        std::sqrt(moments.sum_squares /
                  (static_cast<double>(path_count - 1U) *
                   static_cast<double>(path_count)));
    if (!std::isfinite(standard_error)) {
        throw std::overflow_error{
            "Black-Scholes Monte Carlo standard error is not finite"};
    }
    const double quantile = normal_quantile(0.5 + 0.5 * confidence_level);
    return BlackScholesMonteCarloResult{
        core::OptionPrice{moments.mean},
        core::OptionPrice{standard_error},
        core::OptionPrice{moments.mean - quantile * standard_error},
        core::OptionPrice{moments.mean + quantile * standard_error},
        path_count,
        seed};
}

}  // namespace rvre::pricing
