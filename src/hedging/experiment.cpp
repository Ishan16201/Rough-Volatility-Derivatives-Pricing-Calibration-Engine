#include "rvre/hedging/experiment.hpp"

#include "rvre/models/rough_bergomi.hpp"
#include "rvre/pricing/black_scholes.hpp"
#include "rvre/pricing/monte_carlo.hpp"
#include "rvre/random/random_stream.hpp"
#include "rvre/simulation/time_grid.hpp"

#include <algorithm>
#include <cmath>
#include <fstream>
#include <functional>
#include <iomanip>
#include <limits>
#include <map>
#include <numbers>
#include <numeric>
#include <sstream>
#include <stdexcept>
#include <type_traits>

namespace rvre::hedging {
namespace {

struct MarketPaths {
    simulation::TimeGrid grid;
    std::size_t path_count;
    std::vector<double> assets;

    [[nodiscard]] double at(std::size_t path, std::size_t point) const {
        return assets.at(path * grid.point_count() + point);
    }
};

struct StrategyPricing {
    double initial_delta;
    bool initial_delta_pricer_converged;
    std::function<std::pair<double, bool>(double, double)> delta;
};

struct MarketOptionMark {
    double value;
    double numerical_error;
    std::string error_kind;
};

struct PathStatistics {
    std::vector<double> pnl;
    std::vector<double> turnover;
    std::vector<double> costs;
    double squared_error_sum{0.0};
    double absolute_error_sum{0.0};
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
    constexpr double d4 = 3.754408661698783e+00;
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

[[nodiscard]] double payoff(double spot,
                            double strike,
                            pricing::OptionRight right) noexcept {
    return right == pricing::OptionRight::call
               ? std::max(spot - strike, 0.0)
               : std::max(strike - spot, 0.0);
}

[[nodiscard]] double scenario_spot(const MarketModel& model) {
    return std::visit(
        [](const auto& specification) {
            using T = std::decay_t<decltype(specification)>;
            if constexpr (std::is_same_v<T, BlackScholesModel>) {
                return specification.initial_asset.value;
            } else if constexpr (std::is_same_v<T, models::HestonParameters>) {
                return specification.initial_asset.value;
            } else {
                return specification.parameters.initial_asset.value;
            }
        },
        model);
}

[[nodiscard]] MarketPaths simulate_market_paths(const MarketModel& model,
                                                const simulation::TimeGrid& grid,
                                                std::size_t path_count,
                                                std::uint64_t seed) {
    const std::size_t point_count = grid.point_count();
    if (path_count > std::numeric_limits<std::size_t>::max() / point_count) {
        throw std::length_error{"Hedging market path storage overflows size_t"};
    }
    MarketPaths paths{grid, path_count,
                      std::vector<double>(path_count * point_count, 0.0)};
    std::visit(
        [&](const auto& specification) {
            using T = std::decay_t<decltype(specification)>;
            if constexpr (std::is_same_v<T, BlackScholesModel>) {
                if (!std::isfinite(specification.volatility.value) ||
                    specification.volatility.value < 0.0 ||
                    !std::isfinite(specification.risk_free_rate.value) ||
                    !std::isfinite(specification.dividend_yield.value)) {
                    throw std::invalid_argument{
                        "Invalid Black-Scholes market scenario parameters"};
                }
                for (std::size_t path = 0; path < path_count; ++path) {
                    random::RandomStream stream{
                        seed, static_cast<std::uint64_t>(path)};
                    double spot = scenario_spot(model);
                    paths.assets[path * point_count] = spot;
                    for (std::size_t step = 0; step < grid.interval_count();
                         ++step) {
                        const double dt = grid.interval_years(step);
                        const double vol = specification.volatility.value;
                        spot *= std::exp(
                            (specification.risk_free_rate.value -
                             specification.dividend_yield.value -
                             0.5 * vol * vol) *
                                    dt +
                                vol * std::sqrt(dt) *
                                    stream.standard_normal());
                        paths.assets[path * point_count + step + 1U] = spot;
                    }
                }
            } else if constexpr (std::is_same_v<T, models::HestonParameters>) {
                models::validate_heston_parameters(specification);
                const double orthogonal =
                    std::sqrt(std::max(0.0, 1.0 -
                                               specification.correlation *
                                                   specification.correlation));
                for (std::size_t path = 0; path < path_count; ++path) {
                    random::RandomStream stream{
                        seed, static_cast<std::uint64_t>(path)};
                    double spot = specification.initial_asset.value;
                    double variance = specification.initial_variance.value;
                    paths.assets[path * point_count] = spot;
                    for (std::size_t step = 0; step < grid.interval_count();
                         ++step) {
                        const double dt = grid.interval_years(step);
                        const double root_dt = std::sqrt(dt);
                        const double zv = stream.standard_normal();
                        const double zs = specification.correlation * zv +
                                          orthogonal *
                                              stream.standard_normal();
                        const double positive_variance =
                            std::max(variance, 0.0);
                        spot *= std::exp(
                            (specification.risk_free_rate.value -
                             specification.dividend_yield.value -
                             0.5 * positive_variance) *
                                    dt +
                                std::sqrt(positive_variance) * root_dt * zs);
                        variance +=
                            specification.mean_reversion_speed *
                                (specification.long_run_variance.value -
                                 positive_variance) *
                                dt +
                            specification.volatility_of_variance *
                                std::sqrt(positive_variance) * root_dt * zv;
                        paths.assets[path * point_count + step + 1U] = spot;
                    }
                }
            } else {
                models::RoughBergomiParameters parameters =
                    specification.parameters;
                models::RoughBergomiConfig rough_config{
                    path_count, seed,
                    simulation::VolterraScheme::cell_average_reference};
                const auto simulated =
                    models::simulate_rough_bergomi(parameters, rough_config, grid);
                paths.assets = simulated.asset_values();
            }
        },
        model);
    return paths;
}

[[nodiscard]] double discount_ratio(const models::DiscountCurve& curve,
                                    double from,
                                    double to) {
    const double discount_from = curve.discount_factor(core::Time{from}).value;
    const double discount_to = curve.discount_factor(core::Time{to}).value;
    if (!std::isfinite(discount_from) || discount_from <= 0.0 ||
        !std::isfinite(discount_to) || discount_to <= 0.0) {
        throw std::runtime_error{"Discount curve returned an invalid factor"};
    }
    return discount_from / discount_to;
}

[[nodiscard]] double scenario_growth_factor(const MarketModel& model,
                                           double from,
                                           double to,
                                           bool risk_free) {
    return std::visit(
        [from, to, risk_free](const auto& specification) {
            using T = std::decay_t<decltype(specification)>;
            if constexpr (std::is_same_v<T, BlackScholesModel>) {
                const double rate = risk_free
                                        ? specification.risk_free_rate.value
                                        : specification.dividend_yield.value;
                return std::exp(rate * (to - from));
            } else if constexpr (std::is_same_v<T,
                                                models::HestonParameters>) {
                const double rate = risk_free
                                        ? specification.risk_free_rate.value
                                        : specification.dividend_yield.value;
                return std::exp(rate * (to - from));
            } else {
                const auto& curves = specification.parameters;
                const auto& curve = risk_free
                                        ? *curves.risk_free_discount_curve
                                        : *curves.dividend_discount_curve;
                return discount_ratio(curve, from, to);
            }
        },
        model);
}

[[nodiscard]] MarketOptionMark market_option_mark(
    const MarketModel& market,
    const ExperimentConfig& config,
    std::uint64_t seed) {
    return std::visit(
        [&](const auto& specification) -> MarketOptionMark {
            using T = std::decay_t<decltype(specification)>;
            if constexpr (std::is_same_v<T, BlackScholesModel>) {
                const double value = pricing::black_scholes_price(
                    specification.initial_asset, config.strike,
                    config.maturity,
                    specification.volatility.value *
                        specification.volatility.value *
                        config.maturity.value,
                    core::DiscountFactor{std::exp(
                        -specification.risk_free_rate.value *
                        config.maturity.value)},
                    core::DiscountFactor{std::exp(
                        -specification.dividend_yield.value *
                        config.maturity.value)},
                    config.right);
                return MarketOptionMark{value, 0.0, "analytic"};
            } else if constexpr (std::is_same_v<
                                     T, models::HestonParameters>) {
                const models::HestonCharacteristicFunctionConfig cf_config{
                    400.0, 1.0e-9, 26U};
                const auto result =
                    models::price_heston_characteristic_function(
                        specification, config.strike, config.maturity,
                        config.right, cf_config);
                return MarketOptionMark{
                    result.price.value,
                    result.diagnostics.estimated_quadrature_error,
                    result.diagnostics.converged
                        ? "characteristic_function_error_estimate"
                        : "characteristic_function_nonconverged"};
            } else {
                const auto& parameters = specification.parameters;
                const std::uint64_t mark_seed =
                    seed ^ 0x4D41524B45545F50ULL;
                pricing::MonteCarloPricingConfig mc{
                    config.rough_bergomi_pricing_paths,
                    std::min<std::uint64_t>(
                        512U, config.rough_bergomi_pricing_paths),
                    config.rough_bergomi_pricing_steps,
                    2U,
                    mark_seed,
                    config.confidence_level,
                    models::RoughBergomiConfig{
                        static_cast<std::size_t>(
                            config.rough_bergomi_pricing_paths),
                        mark_seed,
                        simulation::VolterraScheme::cell_average_reference}};
                const auto estimates = pricing::price_rough_bergomi(
                    parameters, mc, {config.strike}, {config.maturity},
                    {pricing::VarianceReduction::standard})
                                           .estimates;
                const auto estimate = std::find_if(
                    estimates.begin(), estimates.end(),
                    [&config](const pricing::PriceEstimate& candidate) {
                        return candidate.option.right == config.right;
                    });
                if (estimate == estimates.end()) {
                    throw std::logic_error{
                        "Rough Bergomi market pricing omitted option right"};
                }
                return MarketOptionMark{
                    estimate->mean.value,
                    estimate->standard_error.value,
                    "monte_carlo_standard_error"};
            }
        },
        market);
}

[[nodiscard]] StrategyPricing make_pricing(
    const HedgeModel& model,
    const ExperimentConfig& config,
    std::uint64_t seed) {
    return std::visit(
        [&](const auto& specification) -> StrategyPricing {
            using T = std::decay_t<decltype(specification)>;
            if constexpr (std::is_same_v<T, BlackScholesModel>) {
                const auto initial_greeks = pricing::black_scholes_greeks(
                    config.initial_asset, config.strike, config.maturity,
                    specification.volatility, specification.risk_free_rate,
                    specification.dividend_yield, config.right);
                return StrategyPricing{
                    initial_greeks.delta, true,
                    [specification, &config](double spot, double time) {
                        const double remaining =
                            config.maturity.value - time;
                        if (remaining <= 0.0) {
                            return std::pair{0.0, true};
                        }
                        return std::pair{
                            pricing::black_scholes_greeks(
                                   core::AssetPrice{spot}, config.strike,
                                   core::Time{remaining},
                                   specification.volatility,
                                   specification.risk_free_rate,
                                   specification.dividend_yield, config.right)
                                .delta,
                            true};
                    }};
            } else if constexpr (std::is_same_v<T,
                                                models::HestonParameters>) {
                models::HestonParameters parameters = specification;
                parameters.initial_asset = config.initial_asset;
                const models::HestonCharacteristicFunctionConfig cf_config{
                    400.0, 1.0e-9, 26U};
                const auto heston_delta =
                    [parameters, &config, cf_config](double spot, double time) {
                        const double remaining =
                            config.maturity.value - time;
                        if (remaining <= 0.0) {
                            return std::pair{0.0, true};
                        }
                        const double bump =
                            std::max(1.0e-4, spot * 1.0e-4);
                        auto up = parameters;
                        auto down = parameters;
                        up.initial_asset = core::AssetPrice{spot + bump};
                        down.initial_asset = core::AssetPrice{spot - bump};
                        const auto up_result =
                            models::price_heston_characteristic_function(
                                up, config.strike, core::Time{remaining},
                                config.right, cf_config);
                        const auto down_result =
                            models::price_heston_characteristic_function(
                                down, config.strike, core::Time{remaining},
                                config.right, cf_config);
                        const double up_price = up_result.price.value;
                        const double down_price = down_result.price.value;
                        return std::pair{
                            (up_price - down_price) / (2.0 * bump),
                            up_result.diagnostics.converged &&
                                down_result.diagnostics.converged};
                    };
                const auto initial_delta = heston_delta(
                    config.initial_asset.value, 0.0);
                return StrategyPricing{initial_delta.first,
                                       initial_delta.second, heston_delta};
            } else {
                const auto& parameters = specification.parameters;
                models::validate_rough_bergomi_parameters(parameters);
                pricing::MonteCarloPricingConfig mc{
                    config.rough_bergomi_pricing_paths,
                    std::min<std::uint64_t>(
                        512U, config.rough_bergomi_pricing_paths),
                    config.rough_bergomi_pricing_steps,
                    2U,
                    seed ^ 0x52425F4845444745ULL,
                    config.confidence_level,
                    models::RoughBergomiConfig{
                        static_cast<std::size_t>(
                            config.rough_bergomi_pricing_paths),
                        seed ^ 0x52425F4845444745ULL,
                        simulation::VolterraScheme::cell_average_reference}};
                const auto estimate = pricing::price_rough_bergomi(
                    parameters, mc, {config.strike}, {config.maturity},
                    {pricing::VarianceReduction::standard})
                                      .estimates;
                const auto selected_estimate =
                    std::find_if(
                    estimate.begin(), estimate.end(),
                    [&config](const pricing::PriceEstimate& candidate) {
                        return candidate.option.right == config.right;
                    });
                if (selected_estimate == estimate.end()) {
                    throw std::logic_error{
                    "Rough Bergomi pricing omitted the requested option right"};
                }
                const auto implied = pricing::black_scholes_implied_volatility(
                    selected_estimate->mean.value, config.initial_asset,
                    config.strike, config.maturity,
                    parameters.risk_free_discount_curve->discount_factor(
                        config.maturity),
                    parameters.dividend_discount_curve->discount_factor(
                        config.maturity),
                    config.right);
                const double frozen_volatility = implied.volatility.value;
                const auto delta =
                    [parameters, &config, frozen_volatility](double spot,
                                                              double time) {
                        const double remaining =
                            config.maturity.value - time;
                        if (remaining <= 0.0) {
                            return std::pair{0.0, true};
                        }
                        const double risk_discount =
                            discount_ratio(
                                *parameters.risk_free_discount_curve, time,
                                config.maturity.value);
                        const double dividend_discount =
                            discount_ratio(
                                *parameters.dividend_discount_curve, time,
                                config.maturity.value);
                        const double rate =
                            -std::log(risk_discount) / remaining;
                        const double dividend =
                            -std::log(dividend_discount) / remaining;
                        return std::pair{
                            pricing::black_scholes_greeks(
                                   core::AssetPrice{spot}, config.strike,
                                   core::Time{remaining},
                                   core::Volatility{frozen_volatility},
                                   core::InterestRate{rate},
                                   core::InterestRate{dividend}, config.right)
                                .delta,
                            true};
                    };
                const auto initial_delta = delta(
                    config.initial_asset.value, 0.0);
                return StrategyPricing{initial_delta.first,
                                       initial_delta.second, delta};
            }
        },
        model);
}

struct LowerTailStatistics {
    double quantile;
    double expected_shortfall_loss;
};

[[nodiscard]] LowerTailStatistics lower_tail_five_percent(
    std::vector<double> sorted_values) {
    std::sort(sorted_values.begin(), sorted_values.end());
    const std::size_t tail_count = std::max<std::size_t>(
        1U, (sorted_values.size() + 19U) / 20U);
    double sum = 0.0;
    for (std::size_t i = 0; i < tail_count; ++i) {
        sum += sorted_values[i];
    }
    const std::size_t quantile_index =
        static_cast<std::size_t>(
            std::floor(0.05 * static_cast<double>(sorted_values.size() - 1U)));
    return LowerTailStatistics{
        sorted_values[quantile_index],
        -sum / static_cast<double>(tail_count)};
}

[[nodiscard]] std::string csv_escape(const std::string& value) {
    if (value.find_first_of(",\"\r\n") == std::string::npos) {
        return value;
    }
    std::string escaped{"\""};
    for (const char character : value) {
        if (character == '"') {
            escaped += "\"\"";
        } else {
            escaped += character;
        }
    }
    escaped += '"';
    return escaped;
}

}  // namespace

HedgingExperimentResult run_hedging_experiment(
    const std::vector<MarketScenario>& scenarios,
    const std::vector<HedgeStrategy>& strategies,
    const ExperimentConfig& config) {
    if (scenarios.empty() || strategies.empty() || config.path_count < 2U ||
        config.path_count > std::numeric_limits<std::size_t>::max() ||
        config.time_step_count == 0U ||
        config.time_step_count > std::numeric_limits<std::size_t>::max() ||
        config.seeds.empty() || config.hedge_every_steps.empty() ||
        !std::isfinite(config.initial_asset.value) ||
        config.initial_asset.value <= 0.0 || !std::isfinite(config.strike) ||
        config.strike <= 0.0 || !std::isfinite(config.maturity.value) ||
        config.maturity.value <= 0.0 ||
        !std::isfinite(config.proportional_transaction_cost) ||
        config.proportional_transaction_cost < 0.0 ||
        !std::isfinite(config.confidence_level) ||
        config.confidence_level <= 0.0 || config.confidence_level >= 1.0 ||
        config.rough_bergomi_pricing_paths < 2U ||
        config.rough_bergomi_pricing_steps == 0U) {
        throw std::invalid_argument{"Invalid hedging experiment configuration"};
    }
    for (const std::uint64_t frequency : config.hedge_every_steps) {
        if (frequency == 0U || frequency > config.time_step_count) {
            throw std::invalid_argument{
                "Hedge frequency must be positive and no greater than the time-step count"};
        }
    }
    for (std::size_t i = 0U; i < config.hedge_every_steps.size(); ++i) {
        if (std::find(config.hedge_every_steps.begin(),
                      config.hedge_every_steps.begin() +
                          static_cast<std::ptrdiff_t>(i),
                      config.hedge_every_steps[i]) !=
            config.hedge_every_steps.begin() +
                static_cast<std::ptrdiff_t>(i)) {
            throw std::invalid_argument{
                "Hedge frequencies must be unique"};
        }
    }
    for (std::size_t i = 0U; i < config.seeds.size(); ++i) {
        if (std::find(config.seeds.begin(),
                      config.seeds.begin() + static_cast<std::ptrdiff_t>(i),
                      config.seeds[i]) !=
            config.seeds.begin() + static_cast<std::ptrdiff_t>(i)) {
            throw std::invalid_argument{
                "Hedging experiment seeds must be unique"};
        }
    }
    for (const auto& scenario : scenarios) {
        if (scenario.name.empty()) {
            throw std::invalid_argument{
                "Market scenario names must not be empty"};
        }
        const double model_spot = scenario_spot(scenario.model);
        if (std::abs(model_spot - config.initial_asset.value) >
            1.0e-10 * config.initial_asset.value) {
            throw std::invalid_argument{
                "Scenario initial spot must equal experiment initial spot"};
        }
    }
    for (const auto& strategy : strategies) {
        if (strategy.name.empty()) {
            throw std::invalid_argument{
                "Hedge strategy names must not be empty"};
        }
        const double strategy_spot = std::visit(
            [](const auto& specification) {
                using T = std::decay_t<decltype(specification)>;
                if constexpr (std::is_same_v<T, BlackScholesModel>) {
                    return specification.initial_asset.value;
                } else if constexpr (std::is_same_v<
                                         T, models::HestonParameters>) {
                    return specification.initial_asset.value;
                } else {
                    return specification.parameters.initial_asset.value;
                }
            },
            strategy.model);
        if (std::abs(strategy_spot - config.initial_asset.value) >
            1.0e-10 * config.initial_asset.value) {
            throw std::invalid_argument{
                "Hedge model initial spot must equal experiment spot"};
        }
    }

    const auto grid = simulation::TimeGrid::uniform(
        config.maturity,
        static_cast<std::size_t>(config.time_step_count));
    HedgingExperimentResult result{config, {}, {}};
    for (const auto& scenario : scenarios) {
        std::map<std::string, std::map<std::uint64_t, PathStatistics>> statistics;
        for (const auto seed : config.seeds) {
            const auto paths = simulate_market_paths(
                scenario.model, grid,
                static_cast<std::size_t>(config.path_count), seed);
            const auto market_mark =
                market_option_mark(scenario.model, config, seed);
            for (const auto& strategy : strategies) {
                const auto pricing =
                    make_pricing(strategy.model, config, seed);
                for (const std::uint64_t frequency : config.hedge_every_steps) {
                    auto& stats = statistics[strategy.name][frequency];
                    for (std::size_t path = 0; path < paths.path_count; ++path) {
                        const double initial_spot = paths.at(path, 0U);
                        double shares = pricing.initial_delta;
                        std::uint64_t nonconverged_delta_rebalances = 0U;
                        double cash = market_mark.value - shares * initial_spot;
                        double turnover = std::abs(shares) * initial_spot;
                        double costs =
                            config.proportional_transaction_cost * turnover;
                        cash -= costs;
                        std::vector<HedgeIntervalRecord> interval_records;
                        for (std::size_t step = 1U;
                             step <= grid.interval_count(); ++step) {
                            const double previous_time =
                                grid.time_years(step - 1U);
                            const double time = grid.time_years(step);
                            const double previous_spot =
                                paths.at(path, step - 1U);
                            const double spot = paths.at(path, step);
                            const double financing =
                                cash * (scenario_growth_factor(
                                            scenario.model, previous_time,
                                            time, true) -
                                        1.0);
                            cash += financing;
                            const double dividend_flow =
                                shares * previous_spot *
                                (scenario_growth_factor(
                                     scenario.model, previous_time, time,
                                     false) -
                                 1.0);
                            cash += dividend_flow;
                            const double old_shares = shares;
                            double traded = 0.0;
                            double transaction_cost = 0.0;
                            std::optional<bool> delta_pricer_converged;
                            const bool rebalanced =
                                step < grid.interval_count() &&
                                step % frequency == 0U;
                            if (rebalanced) {
                                const auto delta = pricing.delta(spot, time);
                                shares = delta.first;
                                delta_pricer_converged = delta.second;
                                if (!delta.second) {
                                    ++nonconverged_delta_rebalances;
                                }
                                traded = shares - old_shares;
                                transaction_cost =
                                    config.proportional_transaction_cost *
                                    std::abs(traded) * spot;
                                cash -= traded * spot + transaction_cost;
                                turnover += std::abs(traded) * spot;
                                costs += transaction_cost;
                            }
                            interval_records.push_back(HedgeIntervalRecord{
                                static_cast<std::uint64_t>(step),
                                core::Time{time},
                                core::AssetPrice{spot},
                                old_shares,
                                shares,
                                traded,
                                financing,
                                dividend_flow,
                                transaction_cost,
                                cash,
                                shares,
                                rebalanced,
                                delta_pricer_converged});
                        }
                        const double final_spot =
                            paths.at(path, grid.point_count() - 1U);
                        const double terminal_payoff =
                            payoff(final_spot, config.strike, config.right);
                        const double pnl =
                            cash + shares * final_spot - terminal_payoff;
                        PathAccounting accounting{
                            scenario.name,
                            strategy.name,
                            frequency,
                            seed,
                            static_cast<std::uint64_t>(path),
                            core::OptionPrice{market_mark.value},
                            core::OptionPrice{
                                market_mark.numerical_error},
                            market_mark.error_kind,
                            core::AssetPrice{initial_spot},
                            pricing.initial_delta,
                            pricing.initial_delta_pricer_converged,
                            nonconverged_delta_rebalances,
                            market_mark.value -
                                pricing.initial_delta * initial_spot -
                                config.proportional_transaction_cost *
                                    std::abs(pricing.initial_delta) *
                                    initial_spot,
                            config.proportional_transaction_cost *
                                std::abs(pricing.initial_delta) * initial_spot,
                            final_spot,
                            cash,
                            shares,
                            terminal_payoff,
                            pnl,
                            turnover,
                            costs,
                            std::move(interval_records)};
                        stats.pnl.push_back(pnl);
                        stats.turnover.push_back(turnover);
                        stats.costs.push_back(costs);
                        stats.squared_error_sum += pnl * pnl;
                        stats.absolute_error_sum += std::abs(pnl);
                        result.path_accounting.push_back(
                            std::move(accounting));
                    }
                }
            }
        }
        for (const auto& strategy : strategies) {
            for (const auto frequency : config.hedge_every_steps) {
                const auto& stats = statistics.at(strategy.name).at(frequency);
                const std::size_t count = stats.pnl.size();
                const double n = static_cast<double>(count);
                const double mean =
                    std::accumulate(stats.pnl.begin(), stats.pnl.end(), 0.0) /
                    n;
                double squared_deviation = 0.0;
                for (const double value : stats.pnl) {
                    squared_deviation += (value - mean) * (value - mean);
                }
                const double standard_deviation =
                    std::sqrt(squared_deviation / (n - 1.0));
                const double standard_error =
                    standard_deviation / std::sqrt(n);
                const double quantile =
                    normal_quantile(0.5 + 0.5 * config.confidence_level);
                const auto average = [n](const std::vector<double>& values) {
                    return std::accumulate(values.begin(), values.end(), 0.0) /
                           n;
                };
                const auto tail = lower_tail_five_percent(stats.pnl);
                result.summaries.push_back(HedgingSummary{
                    scenario.name,
                    strategy.name,
                    frequency,
                    static_cast<std::uint64_t>(count),
                    mean,
                    standard_deviation,
                    std::sqrt(stats.squared_error_sum / n),
                    stats.absolute_error_sum / n,
                    tail.quantile,
                    tail.expected_shortfall_loss,
                    average(stats.turnover),
                    average(stats.costs),
                    standard_error,
                    mean - quantile * standard_error,
                    mean + quantile * standard_error});
            }
        }
    }
    return result;
}

void write_hedging_experiment_csv(const HedgingExperimentResult& result,
                                  const std::string& path) {
    std::ofstream output{path};
    if (!output) {
        throw std::runtime_error{"Unable to open hedging output CSV: " + path};
    }
    output << "record_type,market_scenario,hedge_model,hedge_every_steps,"
              "seed,path_index,path_count,initial_option_value,initial_delta,"
              "terminal_asset,terminal_cash,terminal_shares,terminal_payoff,"
              "hedging_pnl,turnover,transaction_costs,mean_hedging_error,"
              "standard_deviation,rmse,expected_absolute_error,"
              "fifth_percentile_pnl,expected_shortfall_5pct_loss,"
              "mean_turnover,mean_transaction_costs,mean_standard_error,"
              "confidence_interval_lower,confidence_interval_upper,"
              "event,step,time,spot,delta_before,delta_after,shares_traded,"
              "financing_cash_flow,dividend_cash_flow,transaction_cost,"
              "cash_after,shares_after,initial_option_value_error,"
              "initial_option_value_error_kind,time_step_count,maturity,"
              "strike,option_right,transaction_cost_rate,seed_count,"
              "confidence_level,initial_asset,initial_cash_after_hedge,"
              "initial_transaction_cost,initial_delta_pricer_converged,"
              "nonconverged_delta_rebalances,interval_delta_pricer_converged\n";
    const auto number = [](double value) {
        std::ostringstream stream;
        stream << std::setprecision(17) << value;
        return stream.str();
    };
    const auto emit = [&output](const std::vector<std::string>& fields) {
        for (std::size_t i = 0U; i < fields.size(); ++i) {
            if (i != 0U) {
                output << ',';
            }
            output << csv_escape(fields[i]);
        }
        output << '\n';
    };
    for (const auto& summary : result.summaries) {
        std::vector<std::string> row(54U);
        row[0] = "summary";
        row[1] = summary.market_scenario;
        row[2] = summary.hedge_model;
        row[3] = std::to_string(summary.hedge_every_steps);
        row[6] = std::to_string(summary.path_count);
        row[16] = number(summary.mean_hedging_error);
        row[17] = number(summary.standard_deviation);
        row[18] = number(summary.root_mean_square_error);
        row[19] = number(summary.expected_absolute_error);
        row[20] = number(summary.fifth_percentile_pnl);
        row[21] = number(summary.expected_shortfall_five_percent_loss);
        row[22] = number(summary.mean_turnover);
        row[23] = number(summary.mean_transaction_costs);
        row[24] = number(summary.mean_standard_error);
        row[25] = number(summary.confidence_interval_lower);
        row[26] = number(summary.confidence_interval_upper);
        row[41] = std::to_string(result.config.time_step_count);
        row[42] = number(result.config.maturity.value);
        row[43] = number(result.config.strike);
        row[44] = result.config.right == pricing::OptionRight::call
                      ? "call"
                      : "put";
        row[45] = number(result.config.proportional_transaction_cost);
        row[46] = std::to_string(result.config.seeds.size());
        row[47] = number(result.config.confidence_level);
        emit(row);
    }
    for (const auto& path_record : result.path_accounting) {
        std::vector<std::string> row(54U);
        row[0] = "path";
        row[1] = path_record.market_scenario;
        row[2] = path_record.hedge_model;
        row[3] = std::to_string(path_record.hedge_every_steps);
        row[4] = std::to_string(path_record.seed);
        row[5] = std::to_string(path_record.path_index);
        row[6] = std::to_string(result.config.path_count);
        row[7] = number(path_record.initial_option_value.value);
        row[8] = number(path_record.initial_delta);
        row[9] = number(path_record.terminal_asset);
        row[10] = number(path_record.terminal_cash);
        row[11] = number(path_record.terminal_shares);
        row[12] = number(path_record.terminal_payoff);
        row[13] = number(path_record.hedging_pnl);
        row[14] = number(path_record.turnover);
        row[15] = number(path_record.transaction_costs);
        row[27] = "terminal";
        row[39] = number(
            path_record.initial_option_value_numerical_error.value);
        row[40] = path_record.initial_option_value_error_kind;
        row[41] = std::to_string(result.config.time_step_count);
        row[42] = number(result.config.maturity.value);
        row[43] = number(result.config.strike);
        row[44] = result.config.right == pricing::OptionRight::call
                      ? "call"
                      : "put";
        row[45] = number(result.config.proportional_transaction_cost);
        row[46] = std::to_string(result.config.seeds.size());
        row[47] = number(result.config.confidence_level);
        row[48] = number(path_record.initial_asset.value);
        row[49] = number(path_record.initial_cash_after_hedge);
        row[50] = number(path_record.initial_transaction_cost);
        row[51] = path_record.initial_delta_pricer_converged ? "1" : "0";
        row[52] =
            std::to_string(path_record.nonconverged_delta_rebalances);
        emit(row);
        for (const auto& rebalance : path_record.interval_records) {
            std::vector<std::string> event(54U);
            event[0] = rebalance.rebalanced ? "rebalance" : "cash_flow";
            event[1] = path_record.market_scenario;
            event[2] = path_record.hedge_model;
            event[3] = std::to_string(path_record.hedge_every_steps);
            event[4] = std::to_string(path_record.seed);
            event[5] = std::to_string(path_record.path_index);
            event[6] = std::to_string(result.config.path_count);
            event[7] = number(path_record.initial_option_value.value);
            event[8] = number(path_record.initial_delta);
            event[9] = number(path_record.terminal_asset);
            event[10] = number(path_record.terminal_cash);
            event[11] = number(path_record.terminal_shares);
            event[12] = number(path_record.terminal_payoff);
            event[13] = number(path_record.hedging_pnl);
            event[14] = number(path_record.turnover);
            event[15] = number(path_record.transaction_costs);
            event[27] = rebalance.rebalanced ? "rebalance" : "cash_flow";
            event[28] = std::to_string(rebalance.step);
            event[29] = number(rebalance.time.value);
            event[30] = number(rebalance.spot.value);
            event[31] = number(rebalance.delta_before);
            event[32] = number(rebalance.delta_after);
            event[33] = number(rebalance.shares_traded);
            event[34] = number(rebalance.financing_cash_flow);
            event[35] = number(rebalance.dividend_cash_flow);
            event[36] = number(rebalance.transaction_cost);
            event[37] = number(rebalance.cash_after);
            event[38] = number(rebalance.shares_after);
            event[39] = number(
                path_record.initial_option_value_numerical_error.value);
            event[40] = path_record.initial_option_value_error_kind;
            event[41] = std::to_string(result.config.time_step_count);
            event[42] = number(result.config.maturity.value);
            event[43] = number(result.config.strike);
            event[44] = result.config.right == pricing::OptionRight::call
                            ? "call"
                            : "put";
            event[45] =
                number(result.config.proportional_transaction_cost);
            event[46] = std::to_string(result.config.seeds.size());
            event[47] = number(result.config.confidence_level);
            event[48] = number(path_record.initial_asset.value);
            event[49] = number(path_record.initial_cash_after_hedge);
            event[50] = number(path_record.initial_transaction_cost);
            event[51] = path_record.initial_delta_pricer_converged ? "1" : "0";
            event[52] =
                std::to_string(path_record.nonconverged_delta_rebalances);
            if (rebalance.delta_pricer_converged) {
                event[53] = *rebalance.delta_pricer_converged ? "1" : "0";
            }
            emit(event);
        }
    }
    if (!output) {
        throw std::runtime_error{"Failed while writing hedging output CSV"};
    }
}

const char* market_model_name(const MarketModel& model) noexcept {
    return std::visit(
        [](const auto& specification) {
            using T = std::decay_t<decltype(specification)>;
            if constexpr (std::is_same_v<T, BlackScholesModel>) {
                return "black_scholes";
            } else if constexpr (std::is_same_v<T, models::HestonParameters>) {
                return "heston";
            } else {
                return "rough_bergomi";
            }
        },
        model);
}

const char* hedge_model_name(const HedgeModel& model) noexcept {
    return std::visit(
        [](const auto& specification) {
            using T = std::decay_t<decltype(specification)>;
            if constexpr (std::is_same_v<T, BlackScholesModel>) {
                return "black_scholes";
            } else if constexpr (std::is_same_v<T, models::HestonParameters>) {
                return "heston";
            } else {
                return "rough_bergomi";
            }
        },
        model);
}

}  // namespace rvre::hedging
