#include "rvre/models/rough_bergomi.hpp"

#include "rvre/simulation/brownian_motion.hpp"

#include <algorithm>
#include <cmath>
#include <limits>
#include <stdexcept>
#include <string>
#include <utility>

namespace rvre::models {
namespace {

constexpr std::uint64_t bridge_stream_tag = 0x5242554e49425247ULL;

void validate_knot_times(const std::vector<core::Time>& knot_times,
                         std::size_t value_count) {
    if (knot_times.empty() || knot_times.size() != value_count) {
        throw std::invalid_argument{"Curve knot times and values must have equal non-zero size"};
    }
    if (!std::isfinite(knot_times.front().value) || knot_times.front().value != 0.0) {
        throw std::invalid_argument{"Curve knots must start at zero years"};
    }
    for (std::size_t index = 0; index < knot_times.size(); ++index) {
        if (!std::isfinite(knot_times[index].value) ||
            (index > 0U && knot_times[index].value <= knot_times[index - 1U].value)) {
            throw std::invalid_argument{
                "Curve knot times must be finite and strictly increasing"};
        }
    }
}

[[nodiscard]] double interpolate_flat_extrapolation(
    const std::vector<double>& times,
    const std::vector<double>& values,
    core::Time time) {
    if (!std::isfinite(time.value) || time.value < 0.0) {
        throw std::invalid_argument{"Curve evaluation time must be finite and non-negative"};
    }
    if (time.value <= times.front()) {
        return values.front();
    }
    if (time.value >= times.back()) {
        return values.back();
    }
    const auto upper = std::upper_bound(times.begin(), times.end(), time.value);
    const auto upper_index = static_cast<std::size_t>(upper - times.begin());
    const std::size_t lower_index = upper_index - 1U;
    const double weight = (time.value - times[lower_index]) /
                          (times[upper_index] - times[lower_index]);
    return std::lerp(values[lower_index], values[upper_index], weight);
}

[[nodiscard]] std::size_t checked_path_storage(std::size_t paths,
                                               std::size_t points) {
    if (paths == 0U) {
        throw std::invalid_argument{"Rough Bergomi requires at least one path"};
    }
    if (points != 0U && paths > std::numeric_limits<std::size_t>::max() / points) {
        throw std::length_error{"Rough Bergomi path storage size overflows size_t"};
    }
    return paths * points;
}

void validate_model_curves_on_grid(const RoughBergomiParameters& parameters,
                                   const simulation::TimeGrid& grid,
                                   std::vector<double>& forward_variance,
                                   std::vector<double>& risk_free_discount,
                                   std::vector<double>& dividend_discount) {
    const double initial_risk_free_discount =
        parameters.risk_free_discount_curve->discount_factor(core::Time{0.0}).value;
    const double initial_dividend_discount =
        parameters.dividend_discount_curve->discount_factor(core::Time{0.0}).value;
    constexpr double origin_tolerance = 1.0e-12;
    if (!std::isfinite(initial_risk_free_discount) ||
        std::abs(initial_risk_free_discount - 1.0) > origin_tolerance) {
        throw std::invalid_argument{"Risk-free discount curve must equal one at time zero"};
    }
    if (!std::isfinite(initial_dividend_discount) ||
        std::abs(initial_dividend_discount - 1.0) > origin_tolerance) {
        throw std::invalid_argument{"Dividend discount curve must equal one at time zero"};
    }

    const std::size_t point_count = grid.point_count();
    forward_variance.resize(point_count);
    risk_free_discount.resize(point_count);
    dividend_discount.resize(point_count);
    for (std::size_t point = 0; point < point_count; ++point) {
        const core::Time time{grid.time_years(point)};
        const double variance =
            parameters.initial_forward_variance->value(time).value;
        const double risk_free =
            parameters.risk_free_discount_curve->discount_factor(time).value;
        const double dividend =
            parameters.dividend_discount_curve->discount_factor(time).value;
        if (!std::isfinite(variance) || variance <= 0.0) {
            throw std::invalid_argument{
                "Initial forward variance curve must be finite and positive on the simulation grid"};
        }
        if (!std::isfinite(risk_free) || risk_free <= 0.0) {
            throw std::invalid_argument{
                "Risk-free discount curve must be finite and positive on the simulation grid"};
        }
        if (!std::isfinite(dividend) || dividend <= 0.0) {
            throw std::invalid_argument{
                "Dividend discount curve must be finite and positive on the simulation grid"};
        }
        forward_variance[point] = variance;
        risk_free_discount[point] = risk_free;
        dividend_discount[point] = dividend;
    }
}

[[nodiscard]] double empirical_driver_correlation(
    const simulation::CorrelatedBrownianIncrements& correlated) {
    double first_mean = 0.0;
    double second_mean = 0.0;
    double first_square = 0.0;
    double second_square = 0.0;
    double cross_product = 0.0;
    std::size_t count = 0U;
    for (std::size_t path = 0; path < correlated.first.path_count(); ++path) {
        for (std::size_t step = 0; step < correlated.first.interval_count(); ++step) {
            const double scale = std::sqrt(correlated.first.grid().interval_years(step));
            const double first = correlated.first.increment(path, step) / scale;
            const double second = correlated.second.increment(path, step) / scale;
            ++count;
            const double count_value = static_cast<double>(count);
            first_mean += (first - first_mean) / count_value;
            second_mean += (second - second_mean) / count_value;
            first_square += (first - first_mean) * (first - first_mean);
            second_square += (second - second_mean) * (second - second_mean);
            cross_product += (first - first_mean) * (second - second_mean);
        }
    }
    const double covariance =
        cross_product / static_cast<double>(count - 1U);
    const double first_variance =
        first_square / static_cast<double>(count - 1U);
    const double second_variance =
        second_square / static_cast<double>(count - 1U);
    return covariance / std::sqrt(first_variance * second_variance);
}

[[nodiscard]] std::vector<VariancePathDiagnostic> diagnose_forward_variance(
    const std::vector<double>& variance_values,
    const std::vector<double>& forward_variance,
    const simulation::TimeGrid& grid,
    std::size_t path_count) {
    const std::size_t point_count = grid.point_count();
    std::vector<VariancePathDiagnostic> diagnostics;
    diagnostics.reserve(point_count);
    for (std::size_t point = 0; point < point_count; ++point) {
        double mean = 0.0;
        double sum_squares = 0.0;
        std::size_t count = 0U;
        for (std::size_t path = 0; path < path_count; ++path) {
            const double value = variance_values[path * point_count + point];
            ++count;
            const double count_value = static_cast<double>(count);
            const double delta = value - mean;
            mean += delta / count_value;
            sum_squares += delta * (value - mean);
        }
        const double standard_error =
            path_count > 1U
                ? std::sqrt((sum_squares / static_cast<double>(path_count - 1U)) /
                            static_cast<double>(path_count))
                : 0.0;
        diagnostics.push_back(VariancePathDiagnostic{
            .time = core::Time{grid.time_years(point)},
            .target_forward_variance = core::Variance{forward_variance[point]},
            .empirical_mean_variance = core::Variance{mean},
            .mean_standard_error = core::Variance{standard_error}});
    }
    return diagnostics;
}

}  // namespace

PiecewiseLinearForwardVarianceCurve::PiecewiseLinearForwardVarianceCurve(
    std::vector<core::Time> knot_times,
    std::vector<core::Variance> variances) {
    validate_knot_times(knot_times, variances.size());
    knot_times_.reserve(knot_times.size());
    variances_.reserve(variances.size());
    for (std::size_t index = 0; index < knot_times.size(); ++index) {
        if (!std::isfinite(variances[index].value) || variances[index].value <= 0.0) {
            throw std::invalid_argument{
                "Forward variance curve values must be finite and positive"};
        }
        knot_times_.push_back(knot_times[index].value);
        variances_.push_back(variances[index].value);
    }
}

core::Variance PiecewiseLinearForwardVarianceCurve::value(core::Time time) const {
    return core::Variance{
        interpolate_flat_extrapolation(knot_times_, variances_, time)};
}

PiecewiseLinearDiscountCurve::PiecewiseLinearDiscountCurve(
    std::vector<core::Time> knot_times,
    std::vector<core::DiscountFactor> discount_factors) {
    validate_knot_times(knot_times, discount_factors.size());
    knot_times_.reserve(knot_times.size());
    discount_factors_.reserve(discount_factors.size());
    for (std::size_t index = 0; index < knot_times.size(); ++index) {
        if (!std::isfinite(discount_factors[index].value) ||
            discount_factors[index].value <= 0.0) {
            throw std::invalid_argument{
                "Discount curve factors must be finite and positive"};
        }
        knot_times_.push_back(knot_times[index].value);
        discount_factors_.push_back(discount_factors[index].value);
    }
    if (std::abs(discount_factors_.front() - 1.0) > 1.0e-12) {
        throw std::invalid_argument{"Discount curve must equal one at time zero"};
    }
}

core::DiscountFactor PiecewiseLinearDiscountCurve::discount_factor(
    core::Time time) const {
    return core::DiscountFactor{
        interpolate_flat_extrapolation(knot_times_, discount_factors_, time)};
}

RoughBergomiPaths::RoughBergomiPaths(simulation::TimeGrid grid,
                                     simulation::VolterraPaths volterra,
                                     std::size_t path_count,
                                     std::vector<double> variance_values,
                                     std::vector<double> asset_values,
                                     RoughBergomiDiagnostics diagnostics)
    : grid_{std::move(grid)},
      volterra_{std::move(volterra)},
      path_count_{path_count},
      variance_values_{std::move(variance_values)},
      asset_values_{std::move(asset_values)},
      diagnostics_{std::move(diagnostics)} {
    const std::size_t expected = checked_path_storage(path_count_, grid_.point_count());
    if (variance_values_.size() != expected || asset_values_.size() != expected ||
        volterra_.path_count() != path_count_ ||
        volterra_.grid().points_years() != grid_.points_years()) {
        throw std::invalid_argument{"Rough Bergomi result paths have inconsistent dimensions"};
    }
    for (const double value : variance_values_) {
        if (!std::isfinite(value) || value <= 0.0) {
            throw std::invalid_argument{"Rough Bergomi variance paths must be finite and positive"};
        }
    }
    for (const double value : asset_values_) {
        if (!std::isfinite(value) || value <= 0.0) {
            throw std::invalid_argument{"Rough Bergomi asset paths must be finite and positive"};
        }
    }
}

double RoughBergomiPaths::variance(std::size_t path_index,
                                   std::size_t point_index) const {
    if (path_index >= path_count_ || point_index >= grid_.point_count()) {
        throw std::out_of_range{"Rough Bergomi variance path index is out of range"};
    }
    return variance_values_[path_index * grid_.point_count() + point_index];
}

core::AssetPrice RoughBergomiPaths::asset_price(std::size_t path_index,
                                                std::size_t point_index) const {
    if (path_index >= path_count_ || point_index >= grid_.point_count()) {
        throw std::out_of_range{"Rough Bergomi asset path index is out of range"};
    }
    return core::AssetPrice{
        asset_values_[path_index * grid_.point_count() + point_index]};
}

void validate_rough_bergomi_parameters(const RoughBergomiParameters& parameters) {
    if (!std::isfinite(parameters.initial_asset.value) ||
        parameters.initial_asset.value <= 0.0) {
        throw std::invalid_argument{
            "Rough Bergomi initial asset price must be finite and positive"};
    }
    if (!std::isfinite(parameters.hurst.value) ||
        parameters.hurst.value <= 0.0 || parameters.hurst.value >= 0.5) {
        throw std::invalid_argument{"Rough Bergomi H must be finite and in (0, 0.5)"};
    }
    if (!std::isfinite(parameters.eta.value) || parameters.eta.value <= 0.0) {
        throw std::invalid_argument{"Rough Bergomi eta must be finite and positive"};
    }
    if (!std::isfinite(parameters.correlation) ||
        parameters.correlation <= -1.0 || parameters.correlation >= 1.0) {
        throw std::invalid_argument{"Rough Bergomi rho must be finite and in (-1, 1)"};
    }
    if (!parameters.initial_forward_variance) {
        throw std::invalid_argument{"Initial forward variance curve is required"};
    }
    if (!parameters.risk_free_discount_curve) {
        throw std::invalid_argument{"Risk-free discount curve is required"};
    }
    if (!parameters.dividend_discount_curve) {
        throw std::invalid_argument{"Dividend discount curve is required"};
    }
}

RoughBergomiPaths simulate_rough_bergomi(
    const RoughBergomiParameters& parameters,
    const RoughBergomiConfig& config,
    const simulation::TimeGrid& grid) {
    validate_rough_bergomi_parameters(parameters);
    if (config.path_count == 0U) {
        throw std::invalid_argument{"Rough Bergomi path count must be positive"};
    }
    if (config.path_count < 2U) {
        throw std::invalid_argument{
            "Rough Bergomi requires at least two paths to estimate diagnostics"};
    }
    if (config.antithetic_pairs && config.path_count % 2U != 0U) {
        throw std::invalid_argument{
            "Rough Bergomi antithetic simulation requires an even path count"};
    }

    const std::size_t point_count = grid.point_count();
    const std::size_t storage_size = checked_path_storage(config.path_count, point_count);
    std::vector<double> forward_variance;
    std::vector<double> risk_free_discount;
    std::vector<double> dividend_discount;
    validate_model_curves_on_grid(
        parameters, grid, forward_variance, risk_free_discount, dividend_discount);

    const auto correlated_brownian =
        simulation::generate_correlated_brownian_increments(
            grid,
            config.path_count,
            parameters.correlation,
            random::RandomStream{config.seed, config.stream_id ^ 0x524256524542524FULL},
            config.antithetic_pairs);
    simulation::VolterraPaths volterra =
        config.volterra_scheme == simulation::VolterraScheme::cell_average_reference
            ? simulation::simulate_volterra(correlated_brownian.first,
                                            parameters.hurst.value)
            : simulation::simulate_volterra(
                  correlated_brownian.first,
                  parameters.hurst.value,
                  config.volterra_scheme,
                  random::RandomStream{config.seed,
                                       config.stream_id ^ bridge_stream_tag},
                  config.antithetic_pairs);

    std::vector<double> variance_values(storage_size);
    std::vector<double> asset_values(storage_size);
    std::vector<double> discrete_variance(point_count);
    for (std::size_t point = 0; point < point_count; ++point) {
        discrete_variance[point] = simulation::volterra_discrete_variance(
            grid,
            parameters.hurst.value,
            config.volterra_scheme,
            point);
    }
    for (std::size_t path = 0; path < config.path_count; ++path) {
        const std::size_t offset = path * point_count;
        asset_values[offset] = parameters.initial_asset.value;
        for (std::size_t point = 0; point < point_count; ++point) {
            const double driver = volterra.value(path, point);
            const double log_variance =
                parameters.eta.value * driver -
                0.5 * parameters.eta.value * parameters.eta.value *
                    discrete_variance[point];
            const double variance = forward_variance[point] * std::exp(log_variance);
            if (!std::isfinite(log_variance) || !std::isfinite(variance) ||
                variance <= 0.0) {
                throw std::overflow_error{
                    "Rough Bergomi produced a non-finite or non-positive variance at path " +
                    std::to_string(path) + ", point " + std::to_string(point)};
            }
            variance_values[offset + point] = variance;

            if (point + 1U < point_count) {
                const double time_step = grid.interval_years(point);
                const double integrated_rate =
                    std::log(risk_free_discount[point]) -
                    std::log(risk_free_discount[point + 1U]) -
                    std::log(dividend_discount[point]) +
                    std::log(dividend_discount[point + 1U]);
                const double asset_increment =
                    std::sqrt(variance) *
                    correlated_brownian.second.increment(path, point);
                const double log_asset =
                    std::log(asset_values[offset + point]) + integrated_rate -
                    0.5 * variance * time_step + asset_increment;
                const double next_asset = std::exp(log_asset);
                if (!std::isfinite(log_asset) || !std::isfinite(next_asset) ||
                    next_asset <= 0.0) {
                    throw std::overflow_error{
                        "Rough Bergomi produced a non-finite or non-positive asset at path " +
                        std::to_string(path) + ", point " + std::to_string(point + 1U)};
                }
                asset_values[offset + point + 1U] = next_asset;
            }
        }
    }

    const std::size_t final_point = point_count - 1U;
    const double final_numeraire = std::exp(
        std::log(risk_free_discount[final_point]) -
        std::log(dividend_discount[final_point]));
    const double initial_forward_target =
        parameters.initial_asset.value * dividend_discount.front() /
        risk_free_discount.front();
    double adjusted_mean = 0.0;
    double adjusted_sum_squares = 0.0;
    std::size_t count = 0U;
    for (std::size_t path = 0; path < config.path_count; ++path) {
        const double adjusted_asset =
            asset_values[path * point_count + final_point] * final_numeraire;
        if (!std::isfinite(adjusted_asset)) {
            throw std::overflow_error{"Discounted-forward diagnostic is non-finite"};
        }
        ++count;
        const double count_value = static_cast<double>(count);
        const double delta = adjusted_asset - adjusted_mean;
        adjusted_mean += delta / count_value;
        adjusted_sum_squares += delta * (adjusted_asset - adjusted_mean);
    }
    const double adjusted_standard_error =
        std::sqrt((adjusted_sum_squares /
                   static_cast<double>(config.path_count - 1U)) /
                  static_cast<double>(config.path_count));
    RoughBergomiDiagnostics diagnostics{
        .empirical_driver_correlation =
            empirical_driver_correlation(correlated_brownian),
        .forward_variance_diagnostics =
            diagnose_forward_variance(
                variance_values, forward_variance, grid, config.path_count),
        .discounted_forward_mean_asset_price = adjusted_mean,
        .discounted_forward_standard_error = adjusted_standard_error,
        .relative_discounted_forward_error =
            (adjusted_mean - initial_forward_target) / initial_forward_target};

    return RoughBergomiPaths{grid,
                             std::move(volterra),
                             config.path_count,
                             std::move(variance_values),
                             std::move(asset_values),
                             std::move(diagnostics)};
}

}  // namespace rvre::models
