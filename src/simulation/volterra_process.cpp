#include "rvre/simulation/volterra_process.hpp"

#include <array>
#include <cmath>
#include <limits>
#include <numbers>
#include <optional>
#include <stdexcept>
#include <utility>
#include <vector>

namespace rvre::simulation {
namespace {

constexpr std::size_t quadrature_order = 64U;
constexpr std::uint64_t bridge_substream_tag = 0x5256425249444745ULL;

[[nodiscard]] std::size_t checked_storage_size(std::size_t paths, std::size_t points) {
    if (points != 0U && paths > std::numeric_limits<std::size_t>::max() / points) {
        throw std::length_error{"Volterra path storage size overflows size_t"};
    }
    return paths * points;
}

void validate_hurst(double hurst) {
    if (!std::isfinite(hurst) || hurst <= 0.0 || hurst >= 0.5) {
        throw std::invalid_argument{"Volterra Hurst parameter must be in (0, 0.5)"};
    }
}

[[nodiscard]] double cell_average_weight(double cell_end_distance,
                                         double time_step,
                                         double hurst) {
    const double exponent = hurst + 0.5;
    double average_kernel = std::pow(cell_end_distance, exponent);
    if (time_step < cell_end_distance) {
        average_kernel *=
            -std::expm1(exponent * std::log1p(-time_step / cell_end_distance));
    }
    return std::sqrt(2.0 * hurst) * average_kernel /
           (exponent * time_step);
}

[[nodiscard]] double midpoint_weight(double midpoint_distance, double hurst) {
    const double exponent = hurst - 0.5;
    return std::sqrt(2.0 * hurst) * std::pow(midpoint_distance, exponent);
}

[[nodiscard]] std::array<double, quadrature_order> gauss_legendre_unit_nodes(
    std::array<double, quadrature_order>& weights) {
    std::array<double, quadrature_order> nodes{};
    constexpr double tolerance = 2.0e-15;
    const std::size_t half_order = quadrature_order / 2U;

    for (std::size_t index = 0; index < half_order; ++index) {
        double root = std::cos(std::numbers::pi_v<double> *
                               (static_cast<double>(index) + 0.75) /
                               (static_cast<double>(quadrature_order) + 0.5));
        double derivative = 0.0;
        for (std::size_t iteration = 0; iteration < 16U; ++iteration) {
            double previous = 1.0;
            double current = root;
            for (std::size_t degree = 2U; degree <= quadrature_order; ++degree) {
                const double degree_value = static_cast<double>(degree);
                const double next =
                    ((2.0 * degree_value - 1.0) * root * current -
                     (degree_value - 1.0) * previous) /
                    degree_value;
                previous = current;
                current = next;
            }
            derivative = static_cast<double>(quadrature_order) *
                        (root * current - previous) / (root * root - 1.0);
            const double correction = current / derivative;
            root -= correction;
            if (std::abs(correction) <= tolerance) {
                break;
            }
        }
        const double weight_on_minus_one_to_one =
            2.0 / ((1.0 - root * root) * derivative * derivative);
        const std::size_t mirror = quadrature_order - 1U - index;
        nodes[index] = 0.5 * (1.0 - root);
        nodes[mirror] = 0.5 * (1.0 + root);
        weights[index] = 0.5 * weight_on_minus_one_to_one;
        weights[mirror] = weights[index];
    }
    return nodes;
}

}  // namespace

VolterraPaths::VolterraPaths(TimeGrid grid,
                             double hurst,
                             VolterraScheme scheme,
                             std::size_t path_count,
                             std::vector<double> values)
    : grid_{std::move(grid)},
      hurst_{hurst},
      scheme_{scheme},
      path_count_{path_count},
      values_{std::move(values)} {
    validate_hurst(hurst_);
    const std::size_t expected = checked_storage_size(path_count_, grid_.point_count());
    if (path_count_ == 0U || values_.size() != expected) {
        throw std::invalid_argument{"Volterra path storage has an incorrect size"};
    }
    for (const double value : values_) {
        if (!std::isfinite(value)) {
            throw std::invalid_argument{"Volterra path values must be finite"};
        }
    }
}

double VolterraPaths::value(std::size_t path_index, std::size_t point_index) const {
    if (path_index >= path_count_ || point_index >= grid_.point_count()) {
        throw std::out_of_range{"Volterra path index is out of range"};
    }
    return values_[path_index * grid_.point_count() + point_index];
}

namespace {

VolterraPaths simulate_volterra_impl(const BrownianIncrements& brownian,
                                     double hurst,
                                     VolterraScheme scheme,
                                     random::RandomStream* bridge_random,
                                     bool antithetic_pairs) {
    validate_hurst(hurst);
    if (scheme != VolterraScheme::cell_average_reference &&
        scheme != VolterraScheme::hybrid_exact_diagonal) {
        throw std::invalid_argument{"Unsupported Volterra discretization scheme"};
    }
    const TimeGrid& grid = brownian.grid();
    const std::size_t path_count = brownian.path_count();
    const std::size_t point_count = grid.point_count();
    const std::size_t interval_count = grid.interval_count();
    if (antithetic_pairs && path_count % 2U != 0U) {
        throw std::invalid_argument{
            "Antithetic Volterra simulation requires an even path count"};
    }
    std::vector<double> values(checked_storage_size(path_count, point_count), 0.0);
    std::vector<double> reference_lag_weights;
    std::vector<double> hybrid_lag_weights;
    if (const auto uniform_step = grid.uniform_step_years()) {
        if (scheme == VolterraScheme::cell_average_reference) {
            reference_lag_weights.resize(interval_count + 1U);
            for (std::size_t lag = 1U; lag <= interval_count; ++lag) {
                reference_lag_weights[lag] =
                    cell_average_weight(static_cast<double>(lag) * *uniform_step,
                                        *uniform_step,
                                        hurst);
            }
        } else if (scheme == VolterraScheme::hybrid_exact_diagonal) {
            hybrid_lag_weights.resize(interval_count + 1U);
            for (std::size_t lag = 1U; lag <= interval_count; ++lag) {
                hybrid_lag_weights[lag] =
                    midpoint_weight((static_cast<double>(lag) - 0.5) * *uniform_step,
                                    hurst);
            }
        }
    }

    for (std::size_t path = 0; path < path_count; ++path) {
        const std::size_t output_offset = path * point_count;
        std::optional<random::RandomStream> path_bridge_stream;
        if (scheme == VolterraScheme::hybrid_exact_diagonal) {
            if (bridge_random == nullptr) {
                throw std::invalid_argument{
                    "Hybrid Volterra simulation requires an explicit bridge RNG stream"};
            }
            const std::size_t stream_path = antithetic_pairs ? path / 2U : path;
            path_bridge_stream.emplace(bridge_random->substream(bridge_substream_tag)
                                            .substream(
                                                static_cast<std::uint64_t>(stream_path)));
        }
        for (std::size_t point = 1U; point < point_count; ++point) {
            double sum = 0.0;
            const std::size_t diagonal_interval = point - 1U;
            const double diagonal_step = grid.interval_years(diagonal_interval);

            if (scheme == VolterraScheme::cell_average_reference) {
                for (std::size_t interval = 0; interval < point; ++interval) {
                    const double cell_end_distance =
                        grid.time_years(point) - grid.time_years(interval);
                    const std::size_t lag = point - interval;
                    const double weight = reference_lag_weights.empty()
                                              ? cell_average_weight(
                                                    cell_end_distance,
                                                    grid.interval_years(interval),
                                                    hurst)
                                              : reference_lag_weights[lag];
                    sum += weight * brownian.increment(path, interval);
                }
            } else if (scheme == VolterraScheme::hybrid_exact_diagonal) {
                for (std::size_t interval = 0; interval < diagonal_interval; ++interval) {
                    const double midpoint_distance =
                        grid.time_years(point) -
                        0.5 * (grid.time_years(interval) +
                               grid.time_years(interval + 1U));
                    const std::size_t lag = point - interval;
                    const double weight = hybrid_lag_weights.empty()
                                              ? midpoint_weight(midpoint_distance, hurst)
                                              : hybrid_lag_weights[lag];
                    sum += weight * brownian.increment(path, interval);
                }

                const double diagonal_exponent = hurst + 0.5;
                const double diagonal_coefficient =
                    std::sqrt(2.0 * hurst) *
                    std::pow(diagonal_step, hurst - 0.5) / diagonal_exponent;
                const double diagonal_residual_variance =
                    std::pow(diagonal_step, 2.0 * hurst) *
                    std::pow((0.5 - hurst) / diagonal_exponent, 2.0);
                if (diagonal_residual_variance < 0.0) {
                    throw std::runtime_error{
                        "Hybrid Volterra diagonal residual variance is negative"};
                }
                sum += diagonal_coefficient *
                           brownian.increment(path, diagonal_interval) +
                       std::sqrt(diagonal_residual_variance) *
                           (antithetic_pairs && path % 2U != 0U ? -1.0 : 1.0) *
                               path_bridge_stream->standard_normal();
            }
            values[output_offset + point] = sum;
        }
    }

    return VolterraPaths{grid, hurst, scheme, path_count, std::move(values)};
}

}  // namespace

VolterraPaths simulate_volterra(const BrownianIncrements& brownian, double hurst) {
    return simulate_volterra_impl(
        brownian, hurst, VolterraScheme::cell_average_reference, nullptr, false);
}

VolterraPaths simulate_volterra(const BrownianIncrements& brownian,
                                double hurst,
                                VolterraScheme scheme,
                                random::RandomStream bridge_random,
                                bool antithetic_pairs) {
    return simulate_volterra_impl(
        brownian, hurst, scheme, &bridge_random, antithetic_pairs);
}

double volterra_theoretical_covariance(double hurst,
                                       core::Time first_time,
                                       core::Time second_time) {
    validate_hurst(hurst);
    if (!std::isfinite(first_time.value) || !std::isfinite(second_time.value) ||
        first_time.value < 0.0 || second_time.value < 0.0) {
        throw std::invalid_argument{"Covariance times must be finite and non-negative"};
    }
    const double lower_time = std::min(first_time.value, second_time.value);
    const double upper_time = std::max(first_time.value, second_time.value);
    if (lower_time == 0.0) {
        return 0.0;
    }
    if (lower_time == upper_time) {
        return std::pow(lower_time, 2.0 * hurst);
    }

    const double kernel_exponent = hurst - 0.5;
    const double integration_power = 1.0 / (kernel_exponent + 1.0);
    const double separation = upper_time - lower_time;
    std::array<double, quadrature_order> weights{};
    const auto nodes = gauss_legendre_unit_nodes(weights);
    double integral = 0.0;
    for (std::size_t index = 0; index < quadrature_order; ++index) {
        const double transformed_time =
            lower_time * std::pow(nodes[index], integration_power);
        integral += weights[index] *
                    std::pow(separation + transformed_time, kernel_exponent);
    }
    integral *= std::pow(lower_time, kernel_exponent + 1.0) /
                (kernel_exponent + 1.0);
    return 2.0 * hurst * integral;
}

double volterra_discrete_variance(const TimeGrid& grid,
                                  double hurst,
                                  VolterraScheme scheme,
                                  std::size_t point_index) {
    validate_hurst(hurst);
    if (point_index >= grid.point_count()) {
        throw std::out_of_range{"Volterra variance point index is out of range"};
    }
    if (scheme != VolterraScheme::cell_average_reference &&
        scheme != VolterraScheme::hybrid_exact_diagonal) {
        throw std::invalid_argument{"Unsupported Volterra discretization scheme"};
    }
    if (point_index == 0U) {
        return 0.0;
    }

    double variance = 0.0;
    const std::size_t diagonal_interval = point_index - 1U;
    const double diagonal_step = grid.interval_years(diagonal_interval);
    const double exponent = hurst + 0.5;
    for (std::size_t interval = 0; interval < point_index; ++interval) {
        const double cell_end_distance =
            grid.time_years(point_index) - grid.time_years(interval);
        const double weight =
            scheme == VolterraScheme::cell_average_reference
                ? cell_average_weight(
                      cell_end_distance, grid.interval_years(interval), hurst)
                : interval == diagonal_interval
                      ? std::sqrt(2.0 * hurst) *
                            std::pow(diagonal_step, hurst - 0.5) / exponent
                      : midpoint_weight(
                            grid.time_years(point_index) -
                                0.5 * (grid.time_years(interval) +
                                       grid.time_years(interval + 1U)),
                            hurst);
        variance += weight * weight * grid.interval_years(interval);
    }
    if (scheme == VolterraScheme::hybrid_exact_diagonal) {
        variance += std::pow(diagonal_step, 2.0 * hurst) *
                    std::pow((0.5 - hurst) / exponent, 2.0);
    }
    if (!std::isfinite(variance) || variance < 0.0) {
        throw std::runtime_error{"Volterra discrete variance is invalid"};
    }
    return variance;
}

}  // namespace rvre::simulation
