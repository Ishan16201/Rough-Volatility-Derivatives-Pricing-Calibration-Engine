#include "rvre/hedging/experiment.hpp"

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
    std::uint64_t paths{256U};
    std::uint64_t steps{32U};
    std::uint64_t seed{731927U};
    std::uint64_t seed_count{3U};
    double transaction_cost{0.001};
    std::string output{"hedging_experiment.csv"};
};

[[nodiscard]] std::uint64_t parse_integer(std::string_view text,
                                          const char* name) {
    std::uint64_t value = 0U;
    const auto [end, error] =
        std::from_chars(text.data(), text.data() + text.size(), value);
    if (error != std::errc{} || end != text.data() + text.size() ||
        value == 0U) {
        throw std::invalid_argument{std::string{name} +
                                    " must be a positive integer"};
    }
    return value;
}

[[nodiscard]] double parse_double(std::string_view text, const char* name) {
    double value = 0.0;
    const auto [end, error] =
        std::from_chars(text.data(), text.data() + text.size(), value);
    if (error != std::errc{} || end != text.data() + text.size() ||
        !std::isfinite(value) || value < 0.0) {
        throw std::invalid_argument{std::string{name} +
                                    " must be finite and non-negative"};
    }
    return value;
}

[[nodiscard]] Arguments parse_arguments(int argc, char* argv[]) {
    Arguments arguments;
    for (int index = 1; index < argc; ++index) {
        const std::string option{argv[index]};
        if (option == "--help") {
            std::cout
                << "Usage: rvre_hedge [--paths N] [--steps N] [--seed N]\n"
                   "                   [--seed-count N] [--cost RATE] [--output FILE]\n"
                   "Runs paired Black-Scholes, Heston, and Rough Bergomi "
                   "market/hedger comparisons.\n";
            std::exit(0);
        }
        if (index + 1 >= argc) {
            throw std::invalid_argument{"Missing value for " + option};
        }
        const std::string value{argv[++index]};
        if (option == "--paths") {
            arguments.paths = parse_integer(value, "paths");
        } else if (option == "--steps") {
            arguments.steps = parse_integer(value, "steps");
        } else if (option == "--seed") {
            arguments.seed = parse_integer(value, "seed");
        } else if (option == "--seed-count") {
            arguments.seed_count = parse_integer(value, "seed-count");
        } else if (option == "--cost") {
            arguments.transaction_cost = parse_double(value, "cost");
        } else if (option == "--output") {
            arguments.output = value;
        } else {
            throw std::invalid_argument{"Unknown option: " + option};
        }
    }
    if (arguments.paths < 2U || arguments.steps < 4U ||
        arguments.seed_count > 1000U) {
        throw std::invalid_argument{
            "paths must be >= 2, steps >= 4, and seed-count <= 1000"};
    }
    return arguments;
}

}  // namespace

int main(int argc, char* argv[]) {
    try {
        const Arguments arguments = parse_arguments(argc, argv);
        const core::AssetPrice spot{100.0};
        const core::Time maturity{1.0};
        const hedging::BlackScholesModel black_scholes{
            spot, core::Volatility{0.2}, core::InterestRate{0.03},
            core::InterestRate{0.01}};
        const models::HestonParameters heston{
            spot, core::Variance{0.04}, 2.0, core::Variance{0.05}, 0.5,
            -0.7, core::InterestRate{0.03}, core::InterestRate{0.01}};
        auto forward = std::make_shared<
            models::PiecewiseLinearForwardVarianceCurve>(
            std::vector<core::Time>{core::Time{0.0}, maturity},
            std::vector<core::Variance>{core::Variance{0.04},
                                        core::Variance{0.04}});
        auto risk_free = std::make_shared<
            models::PiecewiseLinearDiscountCurve>(
            std::vector<core::Time>{core::Time{0.0}, maturity},
            std::vector<core::DiscountFactor>{
                core::DiscountFactor{1.0},
                core::DiscountFactor{std::exp(-0.03)}});
        auto dividend = std::make_shared<
            models::PiecewiseLinearDiscountCurve>(
            std::vector<core::Time>{core::Time{0.0}, maturity},
            std::vector<core::DiscountFactor>{
                core::DiscountFactor{1.0},
                core::DiscountFactor{std::exp(-0.01)}});
        const hedging::RoughBergomiModel rough{
            models::RoughBergomiParameters{
                spot, models::HurstExponent{0.1},
                models::VolatilityOfVolatility{1.4}, -0.7, forward,
                risk_free, dividend}};
        std::vector<hedging::MarketScenario> scenarios{
            {"Black-Scholes market", black_scholes},
            {"Heston market", heston},
            {"Rough Bergomi market", rough}};
        std::vector<hedging::HedgeStrategy> strategies{
            {"Black-Scholes hedge", black_scholes},
            {"Heston hedge", heston},
            {"Rough Bergomi frozen-IV hedge", rough}};
        const std::uint64_t quarter = std::max<std::uint64_t>(
            1U, arguments.steps / 4U);
        const std::uint64_t half = std::max<std::uint64_t>(
            quarter + 1U, arguments.steps / 2U);
        if (half >= arguments.steps) {
            throw std::invalid_argument{
                "steps must be large enough for three distinct hedge frequencies"};
        }
        std::vector<std::uint64_t> seeds;
        seeds.reserve(static_cast<std::size_t>(arguments.seed_count));
        for (std::uint64_t i = 0U; i < arguments.seed_count; ++i) {
            seeds.push_back(arguments.seed + i);
        }
        const hedging::ExperimentConfig config{
            spot,
            100.0,
            maturity,
            pricing::OptionRight::call,
            arguments.paths,
            arguments.steps,
            {1U, quarter, half},
            std::move(seeds),
            arguments.transaction_cost,
            0.95,
            std::max<std::uint64_t>(arguments.paths, 256U),
            std::max<std::uint64_t>(arguments.steps, 16U)};
        const auto result =
            hedging::run_hedging_experiment(scenarios, strategies, config);
        hedging::write_hedging_experiment_csv(result, arguments.output);
        std::uint64_t nonconverged_delta_evaluations = 0U;
        for (const auto& path : result.path_accounting) {
            nonconverged_delta_evaluations +=
                path.nonconverged_delta_rebalances;
            if (!path.initial_delta_pricer_converged) {
                ++nonconverged_delta_evaluations;
            }
        }
        std::cout << "Paired synthetic out-of-sample hedging experiment"
                  << " paths_per_seed=" << config.path_count
                  << " seeds=" << config.seeds.size()
                  << " steps=" << config.time_step_count
                  << " transaction_cost=" << config.proportional_transaction_cost
                  << " results=" << arguments.output << '\n';
        if (nonconverged_delta_evaluations > 0U) {
            std::cerr << "Warning: " << nonconverged_delta_evaluations
                      << " Heston characteristic-function delta evaluations "
                         "were flagged non-converged; their prices remain in "
                         "the CSV with convergence flags.\n";
        }
        for (const auto& summary : result.summaries) {
            std::cout << summary.market_scenario << " | " << summary.hedge_model
                      << " | every " << summary.hedge_every_steps
                      << " steps: mean=" << summary.mean_hedging_error
                      << " RMSE=" << summary.root_mean_square_error
                      << " ES5% loss="
                      << summary.expected_shortfall_five_percent_loss
                      << " turnover=" << summary.mean_turnover << '\n';
        }
        return 0;
    } catch (const std::exception& error) {
        std::cerr << "Hedging experiment error: " << error.what() << '\n';
        return 1;
    }
}
