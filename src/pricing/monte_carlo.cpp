#include "rvre/pricing/monte_carlo.hpp"

#include "rvre/simulation/time_grid.hpp"

#include <algorithm>
#include <cmath>
#include <limits>
#include <stdexcept>
#include <utility>

namespace rvre::pricing {
namespace {

constexpr std::uint64_t pilot_stream_tag = 0x50494C4F545F4356ULL;
constexpr std::uint64_t evaluation_stream_tag = 0x4556414C5F505249ULL;

struct RunningMoments {
    std::uint64_t count{0U};
    double mean{0.0};
    double sum_squares{0.0};

    void add(double observation) {
        if (!std::isfinite(observation)) {
            throw std::overflow_error{"Monte Carlo payoff is not finite"};
        }
        ++count;
        const double delta = observation - mean;
        mean += delta / static_cast<double>(count);
        sum_squares += delta * (observation - mean);
    }

    [[nodiscard]] double variance() const {
        if (count < 2U) {
            throw std::logic_error{"At least two independent observations are required"};
        }
        const double result = sum_squares / static_cast<double>(count - 1U);
        if (!std::isfinite(result) || result < 0.0) {
            throw std::overflow_error{"Monte Carlo sample variance is not finite"};
        }
        return result;
    }
};

struct RunningCovariance {
    std::uint64_t count{0U};
    double mean_x{0.0};
    double mean_y{0.0};
    double sum_y_squares{0.0};
    double sum_cross_products{0.0};

    void add(double x, double y) {
        if (!std::isfinite(x) || !std::isfinite(y)) {
            throw std::overflow_error{"Control-variate pilot observation is not finite"};
        }
        ++count;
        const double count_value = static_cast<double>(count);
        const double delta_x = x - mean_x;
        const double delta_y = y - mean_y;
        mean_x += delta_x / count_value;
        mean_y += delta_y / count_value;
        sum_y_squares += delta_y * (y - mean_y);
        sum_cross_products += delta_x * (y - mean_y);
    }

    [[nodiscard]] double beta() const {
        if (count < 2U || sum_y_squares <= 0.0) {
            throw std::runtime_error{
                "Control-variate pilot has insufficient control-variable variance"};
        }
        const double result = sum_cross_products / sum_y_squares;
        if (!std::isfinite(result)) {
            throw std::overflow_error{"Control-variate coefficient is not finite"};
        }
        return result;
    }
};

[[nodiscard]] double normal_quantile(double probability) {
    if (probability <= 0.0 || probability >= 1.0) {
        throw std::invalid_argument{"Confidence probability must be in (0, 1)"};
    }

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
    constexpr double lower_tail = 0.02425;
    constexpr double upper_tail = 1.0 - lower_tail;

    if (probability < lower_tail) {
        const double q = std::sqrt(-2.0 * std::log(probability));
        return (((((c1 * q + c2) * q + c3) * q + c4) * q + c5) * q + c6) /
               ((((d1 * q + d2) * q + d3) * q + d4) * q + 1.0);
    }
    if (probability > upper_tail) {
        const double q = std::sqrt(-2.0 * std::log1p(-probability));
        return -(((((c1 * q + c2) * q + c3) * q + c4) * q + c5) * q + c6) /
               ((((d1 * q + d2) * q + d3) * q + d4) * q + 1.0);
    }

    const double q = probability - 0.5;
    const double r = q * q;
    return (((((a1 * r + a2) * r + a3) * r + a4) * r + a5) * r + a6) * q /
           (((((b1 * r + b2) * r + b3) * r + b4) * r + b5) * r + 1.0);
}

[[nodiscard]] double payoff(double terminal_asset,
                            double strike,
                            OptionRight right) noexcept {
    return right == OptionRight::call
               ? std::max(terminal_asset - strike, 0.0)
               : std::max(strike - terminal_asset, 0.0);
}

[[nodiscard]] std::uint64_t mixed_stream_id(std::uint64_t seed,
                                            std::uint64_t tag,
                                            std::uint64_t maturity_index,
                                            std::uint64_t batch_index) noexcept {
    std::uint64_t value = seed ^ tag;
    value ^= maturity_index + 0x9E3779B97F4A7C15ULL + (value << 6U) +
             (value >> 2U);
    value ^= batch_index + 0x9E3779B97F4A7C15ULL + (value << 6U) +
             (value >> 2U);
    value = (value ^ (value >> 30U)) * 0xBF58476D1CE4E5B9ULL;
    value = (value ^ (value >> 27U)) * 0x94D049BB133111EBULL;
    return value ^ (value >> 31U);
}

struct MaturityWorkspace {
    simulation::TimeGrid grid;
    std::vector<std::size_t> instrument_indices;
    double risk_free_discount;
    double control_expectation;
};

[[nodiscard]] std::vector<std::uint64_t> make_batch_sizes(
    std::uint64_t total_paths,
    std::uint64_t target_batch_size,
    bool antithetic) {
    std::vector<std::uint64_t> sizes;
    if (antithetic) {
        for (std::uint64_t remaining = total_paths; remaining > 0U;) {
            const std::uint64_t batch = std::min(target_batch_size, remaining);
            sizes.push_back(batch);
            remaining -= batch;
        }
        return sizes;
    }

    const std::uint64_t batch_count =
        std::min(total_paths / 2U,
                 total_paths / target_batch_size +
                     (total_paths % target_batch_size == 0U ? 0U : 1U));
    const std::uint64_t base_size = total_paths / batch_count;
    const std::uint64_t larger_batch_count = total_paths % batch_count;
    sizes.reserve(static_cast<std::size_t>(batch_count));
    for (std::uint64_t batch = 0U; batch < batch_count; ++batch) {
        sizes.push_back(base_size + (batch < larger_batch_count ? 1U : 0U));
    }
    return sizes;
}

[[nodiscard]] std::vector<EuropeanOption> create_options(
    const std::vector<double>& strikes,
    const std::vector<core::Time>& maturities) {
    std::vector<EuropeanOption> options;
    options.reserve(strikes.size() * maturities.size() * 2U);
    for (const core::Time maturity : maturities) {
        for (const double strike : strikes) {
            options.push_back(EuropeanOption{strike, maturity, OptionRight::call});
            options.push_back(EuropeanOption{strike, maturity, OptionRight::put});
        }
    }
    return options;
}

[[nodiscard]] std::vector<MaturityWorkspace> create_workspaces(
    const models::RoughBergomiParameters& parameters,
    const MonteCarloPricingConfig& config,
    const std::vector<double>& strikes,
    const std::vector<core::Time>& maturities,
    const std::vector<VarianceReduction>& methods) {
    if (strikes.empty() || maturities.empty() || methods.empty()) {
        throw std::invalid_argument{
            "At least one strike, maturity, and variance-reduction method is required"};
    }
    if (config.path_count < 2U || config.batch_size < 2U ||
        config.batch_size > config.path_count ||
        config.path_count > std::numeric_limits<std::size_t>::max() ||
        config.batch_size > std::numeric_limits<std::size_t>::max() ||
        config.time_step_count == 0U ||
        config.time_step_count > std::numeric_limits<std::size_t>::max()) {
        throw std::invalid_argument{
            "Path, batch, and time-step counts must be valid and fit size_t"};
    }
    if (std::find(methods.begin(), methods.end(), VarianceReduction::antithetic) !=
            methods.end() &&
        config.path_count < 4U) {
        throw std::invalid_argument{
            "Antithetic pricing requires at least two independent path pairs"};
    }
    if (!std::isfinite(config.confidence_level) ||
        config.confidence_level <= 0.0 || config.confidence_level >= 1.0) {
        throw std::invalid_argument{"Confidence level must be finite and in (0, 1)"};
    }
    if (std::find(methods.begin(), methods.end(), VarianceReduction::antithetic) !=
            methods.end() &&
        (config.path_count % 2U != 0U || config.batch_size % 2U != 0U)) {
        throw std::invalid_argument{
            "Antithetic evaluation requires even path count and batch size"};
    }
    const bool uses_control =
        std::find(methods.begin(), methods.end(),
                  VarianceReduction::stock_control_variate) != methods.end();
    if (uses_control &&
        (config.pilot_path_count < 2U ||
         config.pilot_path_count > std::numeric_limits<std::size_t>::max())) {
        throw std::invalid_argument{
            "Control-variate pilot count must be at least two and fit size_t"};
    }
    for (std::size_t i = 0; i < methods.size(); ++i) {
        for (std::size_t j = i + 1U; j < methods.size(); ++j) {
            if (methods[i] == methods[j]) {
                throw std::invalid_argument{
                    "Variance-reduction methods must not contain duplicates"};
            }
        }
        if (methods[i] != VarianceReduction::standard &&
            methods[i] != VarianceReduction::antithetic &&
            methods[i] != VarianceReduction::stock_control_variate) {
            throw std::invalid_argument{"Unknown variance-reduction method"};
        }
    }
    for (const double strike : strikes) {
        if (!std::isfinite(strike) || strike <= 0.0) {
            throw std::invalid_argument{"Option strikes must be finite and positive"};
        }
    }
    for (const core::Time maturity : maturities) {
        if (!std::isfinite(maturity.value) || maturity.value <= 0.0) {
            throw std::invalid_argument{"Option maturities must be finite and positive"};
        }
    }
    models::validate_rough_bergomi_parameters(parameters);

    std::vector<MaturityWorkspace> workspaces;
    workspaces.reserve(maturities.size());
    std::size_t instrument_index = 0U;
    for (const core::Time maturity : maturities) {
        MaturityWorkspace workspace{
            .grid = simulation::TimeGrid::uniform(
                maturity, static_cast<std::size_t>(config.time_step_count)),
            .risk_free_discount =
                parameters.risk_free_discount_curve->discount_factor(maturity).value,
            .control_expectation =
                parameters.initial_asset.value *
                parameters.dividend_discount_curve->discount_factor(maturity).value};
        workspace.instrument_indices.reserve(strikes.size() * 2U);
        for (std::size_t strike_index = 0; strike_index < strikes.size();
             ++strike_index) {
            workspace.instrument_indices.push_back(instrument_index++);
            workspace.instrument_indices.push_back(instrument_index++);
        }
        workspaces.push_back(std::move(workspace));
    }
    return workspaces;
}

[[nodiscard]] PriceEstimate make_estimate(const EuropeanOption& option,
                                          VarianceReduction method,
                                          const RunningMoments& moments,
                                          std::uint64_t path_count,
                                          double beta,
                                          double confidence_level) {
    const double standard_error =
        std::sqrt(moments.variance() / static_cast<double>(moments.count));
    const double quantile = normal_quantile(0.5 + 0.5 * confidence_level);
    const double lower = moments.mean - quantile * standard_error;
    const double upper = moments.mean + quantile * standard_error;
    if (!std::isfinite(moments.mean) || !std::isfinite(standard_error) ||
        !std::isfinite(lower) || !std::isfinite(upper)) {
        throw std::overflow_error{"Price estimate or confidence interval is not finite"};
    }
    return PriceEstimate{
        .option = option,
        .method = method,
        .mean = core::OptionPrice{moments.mean},
        .standard_error = core::OptionPrice{standard_error},
        .confidence_interval =
            ConfidenceInterval{lower, upper},
        .path_count = path_count,
        .independent_observation_count = moments.count,
        .control_coefficient = beta};
}

}  // namespace

PricingResult price_rough_bergomi(
    const models::RoughBergomiParameters& parameters,
    const MonteCarloPricingConfig& config,
    const std::vector<double>& strikes,
    const std::vector<core::Time>& maturities,
    const std::vector<VarianceReduction>& methods) {
    auto workspaces =
        create_workspaces(parameters, config, strikes, maturities, methods);
    const auto options = create_options(strikes, maturities);
    const bool uses_control =
        std::find(methods.begin(), methods.end(),
                  VarianceReduction::stock_control_variate) != methods.end();

    PricingResult result{
        .confidence_level = config.confidence_level,
        .seed = config.seed,
        .path_count = config.path_count,
        .batch_size = config.batch_size,
        .pilot_path_count = uses_control ? config.pilot_path_count : 0U,
        .time_step_count = config.time_step_count};
    result.estimates.reserve(options.size() * methods.size());

    for (std::size_t maturity_index = 0; maturity_index < workspaces.size();
         ++maturity_index) {
        auto& workspace = workspaces[maturity_index];
        std::vector<double> coefficients(workspace.instrument_indices.size(), 0.0);

        if (uses_control) {
            std::vector<RunningCovariance> pilot(
                workspace.instrument_indices.size());
            const auto pilot_batches = make_batch_sizes(
                config.pilot_path_count, config.batch_size, false);
            for (std::size_t batch_index = 0; batch_index < pilot_batches.size();
                 ++batch_index) {
                const std::uint64_t batch_paths = pilot_batches[batch_index];
                models::RoughBergomiConfig model_config = config.model;
                model_config.path_count = static_cast<std::size_t>(batch_paths);
                model_config.seed = config.seed;
                model_config.stream_id = mixed_stream_id(
                    config.seed, pilot_stream_tag,
                    static_cast<std::uint64_t>(maturity_index),
                    static_cast<std::uint64_t>(batch_index));
                model_config.antithetic_pairs = false;
                const auto paths = models::simulate_rough_bergomi(
                    parameters, model_config, workspace.grid);
                const std::size_t final_point = workspace.grid.point_count() - 1U;
                for (std::size_t path = 0; path < paths.path_count(); ++path) {
                    const double terminal =
                        paths.asset_price(path, final_point).value;
                    const double discounted_asset =
                        workspace.risk_free_discount * terminal;
                    for (std::size_t local = 0;
                         local < workspace.instrument_indices.size(); ++local) {
                        const auto& option =
                            options[workspace.instrument_indices[local]];
                        const double discounted_payoff =
                            workspace.risk_free_discount *
                            payoff(terminal, option.strike, option.right);
                        pilot[local].add(discounted_payoff, discounted_asset);
                    }
                }
            }
            for (std::size_t local = 0; local < pilot.size(); ++local) {
                coefficients[local] = pilot[local].beta();
            }
        }

        for (const VarianceReduction method : methods) {
            const bool antithetic = method == VarianceReduction::antithetic;
            std::vector<RunningMoments> moments(
                workspace.instrument_indices.size());
            const auto evaluation_batches = make_batch_sizes(
                config.path_count, config.batch_size, antithetic);
            for (std::size_t batch_index = 0;
                 batch_index < evaluation_batches.size(); ++batch_index) {
                const std::uint64_t batch_paths = evaluation_batches[batch_index];
                models::RoughBergomiConfig model_config = config.model;
                model_config.path_count = static_cast<std::size_t>(batch_paths);
                model_config.seed = config.seed;
                model_config.stream_id = mixed_stream_id(
                    config.seed, evaluation_stream_tag,
                    static_cast<std::uint64_t>(maturity_index),
                    static_cast<std::uint64_t>(batch_index));
                model_config.antithetic_pairs = antithetic;
                const auto paths = models::simulate_rough_bergomi(
                    parameters, model_config, workspace.grid);
                const std::size_t final_point = workspace.grid.point_count() - 1U;

                if (antithetic) {
                    for (std::size_t path = 0; path < paths.path_count();
                         path += 2U) {
                        const double first =
                            paths.asset_price(path, final_point).value;
                        const double second =
                            paths.asset_price(path + 1U, final_point).value;
                        for (std::size_t local = 0;
                             local < workspace.instrument_indices.size(); ++local) {
                            const auto& option =
                                options[workspace.instrument_indices[local]];
                            moments[local].add(
                                0.5 * workspace.risk_free_discount *
                                (payoff(first, option.strike, option.right) +
                                 payoff(second, option.strike, option.right)));
                        }
                    }
                } else {
                    for (std::size_t path = 0; path < paths.path_count(); ++path) {
                        const double terminal =
                            paths.asset_price(path, final_point).value;
                        const double discounted_asset =
                            workspace.risk_free_discount * terminal;
                        for (std::size_t local = 0;
                             local < workspace.instrument_indices.size(); ++local) {
                            const auto& option =
                                options[workspace.instrument_indices[local]];
                            double observation =
                                workspace.risk_free_discount *
                                payoff(terminal, option.strike, option.right);
                            if (method ==
                                VarianceReduction::stock_control_variate) {
                                observation += coefficients[local] *
                                    (workspace.control_expectation -
                                     discounted_asset);
                            }
                            moments[local].add(observation);
                        }
                    }
                }
            }

            for (std::size_t local = 0;
                 local < workspace.instrument_indices.size(); ++local) {
                const std::size_t option_index =
                    workspace.instrument_indices[local];
                result.estimates.push_back(make_estimate(
                    options[option_index], method, moments[local],
                    config.path_count,
                    method == VarianceReduction::stock_control_variate
                        ? coefficients[local]
                        : 0.0,
                    config.confidence_level));
            }
        }
    }
    return result;
}

const char* to_string(OptionRight right) noexcept {
    return right == OptionRight::call ? "call" : "put";
}

const char* to_string(VarianceReduction method) noexcept {
    switch (method) {
        case VarianceReduction::standard:
            return "standard";
        case VarianceReduction::antithetic:
            return "antithetic";
        case VarianceReduction::stock_control_variate:
            return "stock_control_variate";
    }
    return "unknown";
}

}  // namespace rvre::pricing
