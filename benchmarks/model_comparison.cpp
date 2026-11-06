#include "rvre/models/heston.hpp"
#include "rvre/pricing/black_scholes.hpp"
#include "rvre/pricing/monte_carlo.hpp"

#include <cmath>
#include <cstdint>
#include <iomanip>
#include <iostream>
#include <memory>
#include <vector>

int main() {
    using namespace rvre;
    constexpr double spot = 100.0;
    constexpr double strike = 100.0;
    constexpr double maturity = 1.0;
    constexpr double rate = 0.03;
    constexpr double dividend = 0.01;
    constexpr double initial_variance = 0.04;
    constexpr std::uint64_t path_count = 6000U;
    constexpr std::uint64_t seed = 261106U;
    constexpr std::uint64_t steps = 32U;

    const double discount_r = std::exp(-rate * maturity);
    const double discount_q = std::exp(-dividend * maturity);
    const double black_scholes = pricing::black_scholes_price(
        core::AssetPrice{spot}, strike, core::Time{maturity},
        initial_variance * maturity, core::DiscountFactor{discount_r},
        core::DiscountFactor{discount_q}, pricing::OptionRight::call);
    const core::BlackScholesParameters black_scholes_parameters{
        .initial_asset = core::AssetPrice{spot},
        .volatility = core::Volatility{std::sqrt(initial_variance)},
        .risk_free_rate = core::InterestRate{rate},
        .dividend_yield = core::InterestRate{dividend},
        .maturity = core::Time{maturity}};
    const auto black_scholes_mc = pricing::price_black_scholes_monte_carlo(
        black_scholes_parameters, strike, pricing::OptionRight::call,
        path_count, seed);

    const models::HestonParameters heston{
        .initial_asset = core::AssetPrice{spot},
        .initial_variance = core::Variance{initial_variance},
        .mean_reversion_speed = 2.0,
        .long_run_variance = core::Variance{initial_variance},
        .volatility_of_variance = 0.5,
        .correlation = -0.7,
        .risk_free_rate = core::InterestRate{rate},
        .dividend_yield = core::InterestRate{dividend}};
    const auto heston_cf = models::price_heston_characteristic_function(
        heston, strike, core::Time{maturity}, pricing::OptionRight::call);
    const auto heston_mc = models::price_heston_monte_carlo(
        heston, strike, core::Time{maturity}, pricing::OptionRight::call,
        models::HestonMonteCarloConfig{path_count, steps, seed, 0.95});

    const auto forward_curve =
        std::make_shared<models::PiecewiseLinearForwardVarianceCurve>(
            std::vector<core::Time>{core::Time{0.0}, core::Time{maturity}},
            std::vector<core::Variance>{core::Variance{initial_variance},
                                        core::Variance{initial_variance}});
    const auto risk_free_curve =
        std::make_shared<models::PiecewiseLinearDiscountCurve>(
            std::vector<core::Time>{core::Time{0.0}, core::Time{maturity}},
            std::vector<core::DiscountFactor>{core::DiscountFactor{1.0},
                                              core::DiscountFactor{discount_r}});
    const auto dividend_curve =
        std::make_shared<models::PiecewiseLinearDiscountCurve>(
            std::vector<core::Time>{core::Time{0.0}, core::Time{maturity}},
            std::vector<core::DiscountFactor>{core::DiscountFactor{1.0},
                                              core::DiscountFactor{discount_q}});
    const models::RoughBergomiParameters rough_bergomi{
        .initial_asset = core::AssetPrice{spot},
        .hurst = models::HurstExponent{0.1},
        .eta = models::VolatilityOfVolatility{1.5},
        .correlation = -0.7,
        .initial_forward_variance = forward_curve,
        .risk_free_discount_curve = risk_free_curve,
        .dividend_discount_curve = dividend_curve};
    const pricing::MonteCarloPricingConfig rough_config{
        .path_count = path_count,
        .batch_size = 1000U,
        .time_step_count = steps,
        .pilot_path_count = 2000U,
        .seed = seed,
        .confidence_level = 0.95,
        .model = models::RoughBergomiConfig{
            .path_count = 2U,
            .seed = seed,
            .volterra_scheme =
                simulation::VolterraScheme::cell_average_reference}};
    const auto rough_mc = pricing::price_rough_bergomi(
        rough_bergomi, rough_config, {strike}, {core::Time{maturity}},
        {pricing::VarianceReduction::standard});

    std::cout << std::setprecision(10)
              << "# spot=" << spot << ",strike=" << strike
              << ",maturity_years=" << maturity << ",rate=" << rate
              << ",dividend=" << dividend
              << ",initial_variance_all_models=" << initial_variance
              << ",paths=" << path_count << ",steps_stochastic_models=" << steps
              << ",seed=" << seed << '\n'
              << "model,method,price,standard_error,ci_lower,ci_upper,notes\n"
              << "Black-Scholes,analytic," << black_scholes
              << ",0," << black_scholes << ',' << black_scholes
              << ",constant volatility sqrt(initial variance)\n"
              << "Black-Scholes,exact-terminal Monte Carlo,"
              << black_scholes_mc.price.value << ','
              << black_scholes_mc.standard_error.value << ','
              << black_scholes_mc.confidence_interval_lower.value << ','
              << black_scholes_mc.confidence_interval_upper.value
              << ",terminal sampling has no time-grid bias\n"
              << "Heston,characteristic function,"
              << heston_cf.price.value << ','
              << heston_cf.diagnostics.estimated_quadrature_error << ','
              << heston_cf.price.value -
                     heston_cf.diagnostics.estimated_quadrature_error
              << ','
              << heston_cf.price.value +
                     heston_cf.diagnostics.estimated_quadrature_error
              << ",adaptive quadrature plus doubled-cutoff sensitivity estimate\n"
              << "Heston,full truncation Euler Monte Carlo,"
              << heston_mc.price.value << ','
              << heston_mc.standard_error.value << ','
              << heston_mc.confidence_interval_lower.value << ','
              << heston_mc.confidence_interval_upper.value
              << ",time-step bias and Monte Carlo error\n"
              << "Rough Bergomi,cell-average reference Monte Carlo,"
              << rough_mc.estimates.front().mean.value << ','
              << rough_mc.estimates.front().standard_error.value << ','
              << rough_mc.estimates.front().confidence_interval.lower << ','
              << rough_mc.estimates.front().confidence_interval.upper
              << ",time-step bias and Monte Carlo error\n"
              << "# Heston parameters: kappa=2,theta=0.04,sigma=0.5,rho=-0.7\n"
              << "# Rough Bergomi parameters: H=0.1,eta=1.5,rho=-0.7\n";
}
