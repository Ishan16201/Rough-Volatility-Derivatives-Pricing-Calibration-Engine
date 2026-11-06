#include "rvre/pricing/monte_carlo.hpp"

#include <chrono>
#include <cmath>
#include <cstdint>
#include <iomanip>
#include <iostream>
#include <memory>
#include <vector>

int main() {
    using namespace rvre;
    constexpr double maturity = 1.0;
    const models::RoughBergomiParameters model{
        .initial_asset = core::AssetPrice{100.0},
        .hurst = models::HurstExponent{0.1},
        .eta = models::VolatilityOfVolatility{1.5},
        .correlation = -0.7,
        .initial_forward_variance =
            std::make_shared<models::PiecewiseLinearForwardVarianceCurve>(
                std::vector<core::Time>{core::Time{0.0}, core::Time{maturity}},
                std::vector<core::Variance>{core::Variance{0.04},
                                            core::Variance{0.04}}),
        .risk_free_discount_curve =
            std::make_shared<models::PiecewiseLinearDiscountCurve>(
                std::vector<core::Time>{core::Time{0.0}, core::Time{maturity}},
                std::vector<core::DiscountFactor>{
                    core::DiscountFactor{1.0},
                    core::DiscountFactor{std::exp(-0.03 * maturity)}}),
        .dividend_discount_curve =
            std::make_shared<models::PiecewiseLinearDiscountCurve>(
                std::vector<core::Time>{core::Time{0.0}, core::Time{maturity}},
                std::vector<core::DiscountFactor>{
                    core::DiscountFactor{1.0},
                    core::DiscountFactor{std::exp(-0.01 * maturity)}})};
    const pricing::MonteCarloPricingConfig config{
        .path_count = 5000U,
        .batch_size = 1000U,
        .time_step_count = 32U,
        .pilot_path_count = 2000U,
        .seed = 20261106U,
        .confidence_level = 0.95,
        .model = models::RoughBergomiConfig{
            .path_count = 2U,
            .seed = 0U,
            .volterra_scheme = simulation::VolterraScheme::cell_average_reference}};

    std::cout << "method,price,se,ci_lower,ci_upper,seconds\n"
              << std::setprecision(10);
    for (const auto method :
         {pricing::VarianceReduction::standard,
          pricing::VarianceReduction::antithetic,
          pricing::VarianceReduction::stock_control_variate}) {
        const auto start = std::chrono::steady_clock::now();
        const auto result = pricing::price_rough_bergomi(
            model, config, {100.0}, {core::Time{maturity}}, {method});
        const double seconds = std::chrono::duration<double>(
                                   std::chrono::steady_clock::now() - start)
                                   .count();
        const auto& estimate = result.estimates.front();
        std::cout << pricing::to_string(method) << ',' << estimate.mean.value << ','
                  << estimate.standard_error.value << ','
                  << estimate.confidence_interval.lower << ','
                  << estimate.confidence_interval.upper << ',' << seconds << '\n';
    }
}
