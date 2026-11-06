#include "rvre/models/rough_bergomi.hpp"

#include <charconv>
#include <cmath>
#include <cstddef>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <memory>
#include <stdexcept>
#include <string>
#include <vector>

namespace {

[[nodiscard]] std::size_t parse_positive_size(const char* text, const char* name) {
    const std::string input{text};
    std::size_t value = 0U;
    const auto [end, error] =
        std::from_chars(input.data(), input.data() + input.size(), value);
    if (error != std::errc{} || end != input.data() + input.size() || value == 0U) {
        throw std::invalid_argument{std::string{name} + " must be a positive integer"};
    }
    return value;
}

void write_paths(std::ostream& output,
                 const rvre::models::RoughBergomiPaths& paths) {
    output << "path,time_years,asset_price,variance,volterra_driver\n";
    output << std::setprecision(17);
    for (std::size_t path = 0; path < paths.path_count(); ++path) {
        for (std::size_t point = 0; point < paths.grid().point_count(); ++point) {
            output << path << ',' << paths.grid().time_years(point) << ','
                   << paths.asset_price(path, point).value << ','
                   << paths.variance(path, point) << ','
                   << paths.volterra().value(path, point) << '\n';
        }
    }
}

}  // namespace

int main(int argc, char* argv[]) {
    if (argc > 4) {
        std::cerr << "Usage: rvre_rough_bergomi_example [output.csv] [paths] [time_steps]\n";
        return 2;
    }

    try {
        const std::size_t path_count =
            argc > 2 ? parse_positive_size(argv[2], "paths") : 8U;
        const std::size_t time_steps =
            argc > 3 ? parse_positive_size(argv[3], "time_steps") : 32U;
        auto forward_variance =
            std::make_shared<rvre::models::PiecewiseLinearForwardVarianceCurve>(
                std::vector<rvre::core::Time>{rvre::core::Time{0.0},
                                              rvre::core::Time{2.0}},
                std::vector<rvre::core::Variance>{rvre::core::Variance{0.04},
                                                  rvre::core::Variance{0.04}});
        auto risk_free_discount =
            std::make_shared<rvre::models::PiecewiseLinearDiscountCurve>(
                std::vector<rvre::core::Time>{rvre::core::Time{0.0},
                                              rvre::core::Time{2.0}},
                std::vector<rvre::core::DiscountFactor>{
                    rvre::core::DiscountFactor{1.0},
                    rvre::core::DiscountFactor{std::exp(-0.03 * 2.0)}});
        auto dividend_discount =
            std::make_shared<rvre::models::PiecewiseLinearDiscountCurve>(
                std::vector<rvre::core::Time>{rvre::core::Time{0.0},
                                              rvre::core::Time{2.0}},
                std::vector<rvre::core::DiscountFactor>{
                    rvre::core::DiscountFactor{1.0},
                    rvre::core::DiscountFactor{std::exp(-0.01 * 2.0)}});

        const rvre::models::RoughBergomiParameters parameters{
            .initial_asset = rvre::core::AssetPrice{100.0},
            .hurst = rvre::models::HurstExponent{0.1},
            .eta = rvre::models::VolatilityOfVolatility{1.5},
            .correlation = -0.7,
            .initial_forward_variance = std::move(forward_variance),
            .risk_free_discount_curve = std::move(risk_free_discount),
            .dividend_discount_curve = std::move(dividend_discount)};
        const rvre::models::RoughBergomiConfig config{
            .path_count = path_count,
            .seed = 20261106U,
            .volterra_scheme = rvre::simulation::VolterraScheme::cell_average_reference};
        const auto grid = rvre::simulation::TimeGrid::uniform(
            rvre::core::Time{1.0}, time_steps);
        const auto paths =
            rvre::models::simulate_rough_bergomi(parameters, config, grid);

        if (argc > 1) {
            std::ofstream output{argv[1]};
            if (!output) {
                throw std::runtime_error{"Unable to open output CSV: " +
                                         std::string{argv[1]}};
            }
            write_paths(output, paths);
            output.close();
            if (!output) {
                throw std::runtime_error{"Failed while writing output CSV: " +
                                         std::string{argv[1]}};
            }
            std::cout << "Wrote " << path_count << " asset/variance paths to "
                      << argv[1] << '\n';
        } else {
            write_paths(std::cout, paths);
        }

        const auto& diagnostics = paths.diagnostics();
        std::cerr << std::setprecision(8)
                  << "seed=" << config.seed << " paths=" << path_count
                  << " time_steps=" << time_steps << " H=" << parameters.hurst.value
                  << " eta=" << parameters.eta.value
                  << " rho=" << parameters.correlation << '\n'
                  << "empirical_driver_correlation="
                  << diagnostics.empirical_driver_correlation << '\n'
                  << "discounted_forward_mean="
                  << diagnostics.discounted_forward_mean_asset_price << '\n'
                  << "discounted_forward_standard_error="
                  << diagnostics.discounted_forward_standard_error << '\n'
                  << "relative_discounted_forward_error="
                  << diagnostics.relative_discounted_forward_error << '\n';
        for (const auto& variance_diagnostic :
             diagnostics.forward_variance_diagnostics) {
            std::cerr << "forward_variance time="
                      << variance_diagnostic.time.value
                      << " target="
                      << variance_diagnostic.target_forward_variance.value
                      << " empirical_mean="
                      << variance_diagnostic.empirical_mean_variance.value
                      << " standard_error="
                      << variance_diagnostic.mean_standard_error.value << '\n';
        }
        return 0;
    } catch (const std::exception& error) {
        std::cerr << "Rough Bergomi example error: " << error.what() << '\n';
        return 1;
    }
}
