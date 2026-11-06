#include "rvre/pricing/monte_carlo.hpp"

#include <algorithm>
#include <charconv>
#include <cmath>
#include <cstdint>
#include <iomanip>
#include <iostream>
#include <memory>
#include <stdexcept>
#include <string>
#include <string_view>
#include <vector>

namespace {

using namespace rvre;

[[nodiscard]] double parse_double(std::string_view text, const char* name) {
    std::string input{text};
    std::size_t parsed = 0U;
    double value = 0.0;
    try {
        value = std::stod(input, &parsed);
    } catch (const std::exception&) {
        throw std::invalid_argument{std::string{name} + " must be a number"};
    }
    if (parsed != input.size() || !std::isfinite(value)) {
        throw std::invalid_argument{std::string{name} + " must be finite"};
    }
    return value;
}

[[nodiscard]] std::uint64_t parse_integer(std::string_view text,
                                          const char* name) {
    std::uint64_t value = 0U;
    const auto [end, error] =
        std::from_chars(text.data(), text.data() + text.size(), value);
    if (error != std::errc{} || end != text.data() + text.size()) {
        throw std::invalid_argument{std::string{name} + " must be an unsigned integer"};
    }
    return value;
}

[[nodiscard]] std::vector<double> parse_list(std::string_view text,
                                             const char* name) {
    std::vector<double> values;
    std::size_t start = 0U;
    while (start < text.size()) {
        const std::size_t comma = text.find(',', start);
        const std::size_t end =
            comma == std::string_view::npos ? text.size() : comma;
        if (end == start) {
            throw std::invalid_argument{std::string{name} + " contains an empty item"};
        }
        values.push_back(parse_double(text.substr(start, end - start), name));
        if (comma == std::string_view::npos) {
            start = text.size();
        } else {
            if (comma + 1U == text.size()) {
                throw std::invalid_argument{std::string{name} + " has a trailing comma"};
            }
            start = comma + 1U;
        }
    }
    if (values.empty()) {
        throw std::invalid_argument{std::string{name} + " must not be empty"};
    }
    return values;
}

[[nodiscard]] std::vector<pricing::VarianceReduction> parse_methods(
    const std::string& value) {
    std::vector<pricing::VarianceReduction> methods;
    std::size_t start = 0U;
    while (start < value.size()) {
        const std::size_t comma = value.find(',', start);
        const std::string item = value.substr(
            start, comma == std::string::npos ? std::string::npos : comma - start);
        if (item.empty()) {
            throw std::invalid_argument{"methods must not contain empty items"};
        }
        if (item == "standard") {
            methods.push_back(pricing::VarianceReduction::standard);
        } else if (item == "antithetic") {
            methods.push_back(pricing::VarianceReduction::antithetic);
        } else if (item == "control") {
            methods.push_back(pricing::VarianceReduction::stock_control_variate);
        } else {
            throw std::invalid_argument{
                "methods must be a comma-separated subset of standard,antithetic,control"};
        }
        if (comma == std::string::npos) {
            start = value.size();
        } else {
            if (comma + 1U == value.size()) {
                throw std::invalid_argument{"methods has a trailing comma"};
            }
            start = comma + 1U;
        }
    }
    return methods;
}

void print_help() {
    std::cout
        << "Usage: rvre_price [--spot N] [--hurst N] [--eta N] [--rho N] "
           "[--variance N] [--rate N] [--dividend N]\n"
           "                  [--strikes K1,K2] [--maturities T1,T2] "
           "[--paths N] [--batch-size N]\n"
           "                  [--pilot-paths N] [--steps N] [--seed N] "
           "[--confidence N]\n"
           "                  [--methods standard,antithetic,control] "
           "[--scheme reference|hybrid]\n";
}

[[nodiscard]] std::vector<pricing::VarianceReduction> default_methods() {
    return {pricing::VarianceReduction::standard,
            pricing::VarianceReduction::antithetic,
            pricing::VarianceReduction::stock_control_variate};
}

}  // namespace

int main(int argc, char* argv[]) {
    try {
        double spot = 100.0;
        double hurst = 0.1;
        double eta = 1.5;
        double rho = -0.7;
        double variance = 0.04;
        double rate = 0.03;
        double dividend = 0.01;
        double confidence = 0.95;
        std::vector<double> strikes{90.0, 100.0, 110.0};
        std::vector<double> maturity_values{0.5, 1.0};
        std::uint64_t paths = 20000U;
        std::uint64_t batch_size = 2000U;
        std::uint64_t pilot_paths = 10000U;
        std::uint64_t steps = 64U;
        std::uint64_t seed = 20261106U;
        auto methods = default_methods();
        simulation::VolterraScheme scheme =
            simulation::VolterraScheme::cell_average_reference;

        for (int argument = 1; argument < argc; ++argument) {
            const std::string option{argv[argument]};
            if (option == "--help") {
                print_help();
                return 0;
            }
            if (argument + 1 >= argc) {
                throw std::invalid_argument{"Missing value for " + option};
            }
            const std::string value{argv[++argument]};
            if (option == "--spot") {
                spot = parse_double(value, "spot");
            } else if (option == "--hurst") {
                hurst = parse_double(value, "hurst");
            } else if (option == "--eta") {
                eta = parse_double(value, "eta");
            } else if (option == "--rho") {
                rho = parse_double(value, "rho");
            } else if (option == "--variance") {
                variance = parse_double(value, "variance");
            } else if (option == "--rate") {
                rate = parse_double(value, "rate");
            } else if (option == "--dividend") {
                dividend = parse_double(value, "dividend");
            } else if (option == "--confidence") {
                confidence = parse_double(value, "confidence");
            } else if (option == "--strikes") {
                strikes = parse_list(value, "strikes");
            } else if (option == "--maturities") {
                maturity_values = parse_list(value, "maturities");
            } else if (option == "--paths") {
                paths = parse_integer(value, "paths");
            } else if (option == "--batch-size") {
                batch_size = parse_integer(value, "batch-size");
            } else if (option == "--pilot-paths") {
                pilot_paths = parse_integer(value, "pilot-paths");
            } else if (option == "--steps") {
                steps = parse_integer(value, "steps");
            } else if (option == "--seed") {
                seed = parse_integer(value, "seed");
            } else if (option == "--methods") {
                methods = parse_methods(value);
            } else if (option == "--scheme" && value == "reference") {
                scheme = simulation::VolterraScheme::cell_average_reference;
            } else if (option == "--scheme" && value == "hybrid") {
                scheme = simulation::VolterraScheme::hybrid_exact_diagonal;
            } else {
                throw std::invalid_argument{"Unknown option or value: " + option};
            }
        }

        std::vector<core::Time> maturities;
        maturities.reserve(maturity_values.size());
        double maximum_maturity = 0.0;
        for (const double maturity : maturity_values) {
            maturities.push_back(core::Time{maturity});
            maximum_maturity = std::max(maximum_maturity, maturity);
        }
        if (maximum_maturity <= 0.0) {
            throw std::invalid_argument{"Maturities must be positive"};
        }
        const auto forward_curve =
            std::make_shared<models::PiecewiseLinearForwardVarianceCurve>(
                std::vector<core::Time>{core::Time{0.0},
                                        core::Time{maximum_maturity}},
                std::vector<core::Variance>{core::Variance{variance},
                                            core::Variance{variance}});
        const auto risk_free_curve =
            std::make_shared<models::PiecewiseLinearDiscountCurve>(
                std::vector<core::Time>{core::Time{0.0},
                                        core::Time{maximum_maturity}},
                std::vector<core::DiscountFactor>{
                    core::DiscountFactor{1.0},
                    core::DiscountFactor{std::exp(-rate * maximum_maturity)}});
        const auto dividend_curve =
            std::make_shared<models::PiecewiseLinearDiscountCurve>(
                std::vector<core::Time>{core::Time{0.0},
                                        core::Time{maximum_maturity}},
                std::vector<core::DiscountFactor>{
                    core::DiscountFactor{1.0},
                    core::DiscountFactor{std::exp(-dividend * maximum_maturity)}});
        const models::RoughBergomiParameters model{
            .initial_asset = core::AssetPrice{spot},
            .hurst = models::HurstExponent{hurst},
            .eta = models::VolatilityOfVolatility{eta},
            .correlation = rho,
            .initial_forward_variance = forward_curve,
            .risk_free_discount_curve = risk_free_curve,
            .dividend_discount_curve = dividend_curve};
        const pricing::MonteCarloPricingConfig config{
            .path_count = paths,
            .batch_size = batch_size,
            .time_step_count = steps,
            .pilot_path_count = pilot_paths,
            .seed = seed,
            .confidence_level = confidence,
            .model = models::RoughBergomiConfig{
                .path_count = 2U,
                .seed = seed,
                .volterra_scheme = scheme}};
        const auto result =
            pricing::price_rough_bergomi(model, config, strikes, maturities, methods);

        std::cout << std::setprecision(12)
                  << "# seed=" << result.seed << ",paths=" << result.path_count
                  << ",batch_size=" << result.batch_size
                  << ",pilot_paths=" << result.pilot_path_count
                  << ",time_steps=" << result.time_step_count
                  << ",confidence=" << result.confidence_level
                  << ",spot=" << spot << ",forward_variance=" << variance
                  << ",rate=" << rate << ",dividend=" << dividend
                  << ",H=" << hurst << ",eta=" << eta << ",rho=" << rho
                  << ",scheme="
                  << (scheme == simulation::VolterraScheme::cell_average_reference
                          ? "reference"
                          : "hybrid")
                  << '\n'
                  << "right,strike,maturity,method,price,standard_error,ci_lower,"
                     "ci_upper,independent_observations,control_beta\n";
        for (const auto& estimate : result.estimates) {
            std::cout << pricing::to_string(estimate.option.right) << ','
                      << estimate.option.strike << ','
                      << estimate.option.maturity.value << ','
                      << pricing::to_string(estimate.method) << ','
                      << estimate.mean.value << ','
                      << estimate.standard_error.value << ','
                      << estimate.confidence_interval.lower << ','
                      << estimate.confidence_interval.upper << ','
                      << estimate.independent_observation_count << ','
                      << estimate.control_coefficient << '\n';
        }
        return 0;
    } catch (const std::exception& error) {
        std::cerr << "Pricing error: " << error.what() << '\n';
        return 1;
    }
}
