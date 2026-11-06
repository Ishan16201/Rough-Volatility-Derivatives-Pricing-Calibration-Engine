#include "rvre/simulation/brownian_motion.hpp"

#include <cmath>
#include <limits>
#include <stdexcept>
#include <utility>

namespace rvre::simulation {
namespace {

[[nodiscard]] std::size_t checked_size(std::size_t paths, std::size_t points) {
    if (paths == 0U) {
        throw std::invalid_argument{"Brownian simulation requires at least one path"};
    }
    if (points != 0U && paths > std::numeric_limits<std::size_t>::max() / points) {
        throw std::length_error{"Brownian path storage size overflows size_t"};
    }
    return paths * points;
}

void validate_finite(const std::vector<double>& values, const char* message) {
    for (const double value : values) {
        if (!std::isfinite(value)) {
            throw std::invalid_argument{message};
        }
    }
}

}  // namespace

BrownianIncrements::BrownianIncrements(TimeGrid grid,
                                       std::size_t path_count,
                                       std::vector<double> increments)
    : grid_{std::move(grid)}, path_count_{path_count}, increments_{std::move(increments)} {
    const auto expected = checked_size(path_count_, grid_.interval_count());
    if (increments_.size() != expected) {
        throw std::invalid_argument{"Brownian increment storage has an incorrect size"};
    }
    validate_finite(increments_, "Brownian increments must be finite");
}

double BrownianIncrements::increment(std::size_t path_index,
                                     std::size_t interval_index) const {
    if (path_index >= path_count_ || interval_index >= grid_.interval_count()) {
        throw std::out_of_range{"Brownian increment index is out of range"};
    }
    return increments_[path_index * grid_.interval_count() + interval_index];
}

BrownianPaths::BrownianPaths(TimeGrid grid,
                             std::size_t path_count,
                             std::vector<double> values)
    : grid_{std::move(grid)}, path_count_{path_count}, values_{std::move(values)} {
    const auto expected = checked_size(path_count_, grid_.point_count());
    if (values_.size() != expected) {
        throw std::invalid_argument{"Brownian path storage has an incorrect size"};
    }
    validate_finite(values_, "Brownian path values must be finite");
}

double BrownianPaths::value(std::size_t path_index, std::size_t point_index) const {
    if (path_index >= path_count_ || point_index >= grid_.point_count()) {
        throw std::out_of_range{"Brownian path index is out of range"};
    }
    return values_[path_index * grid_.point_count() + point_index];
}

BrownianIncrements generate_brownian_increments(const TimeGrid& grid,
                                                std::size_t path_count,
                                                random::RandomStream base_stream) {
    const std::size_t interval_count = grid.interval_count();
    std::vector<double> increments(checked_size(path_count, interval_count));
    for (std::size_t path = 0; path < path_count; ++path) {
        auto stream = base_stream.substream(static_cast<std::uint64_t>(path));
        const std::size_t offset = path * interval_count;
        for (std::size_t step = 0; step < interval_count; ++step) {
            increments[offset + step] =
                std::sqrt(grid.interval_years(step)) * stream.standard_normal();
        }
    }
    return BrownianIncrements{grid, path_count, std::move(increments)};
}

BrownianPaths cumulative_brownian_paths(const BrownianIncrements& increments) {
    const TimeGrid& grid = increments.grid();
    const std::size_t points = grid.point_count();
    std::vector<double> values(checked_size(increments.path_count(), points), 0.0);
    for (std::size_t path = 0; path < increments.path_count(); ++path) {
        const std::size_t offset = path * points;
        for (std::size_t step = 0; step < grid.interval_count(); ++step) {
            values[offset + step + 1U] =
                values[offset + step] + increments.increment(path, step);
        }
    }
    return BrownianPaths{grid, increments.path_count(), std::move(values)};
}

CorrelatedBrownianIncrements generate_correlated_brownian_increments(
    const TimeGrid& grid,
    std::size_t path_count,
    double correlation,
    random::RandomStream base_stream,
    bool antithetic_pairs) {
    if (!std::isfinite(correlation) || correlation < -1.0 || correlation > 1.0) {
        throw std::invalid_argument{"Brownian correlation must be finite and in [-1, 1]"};
    }
    if (antithetic_pairs && path_count % 2U != 0U) {
        throw std::invalid_argument{
            "Antithetic Brownian simulation requires an even path count"};
    }
    const std::size_t interval_count = grid.interval_count();
    const std::size_t storage_size = checked_size(path_count, interval_count);
    std::vector<double> first(storage_size);
    std::vector<double> second(storage_size);
    const double orthogonal_scale = std::sqrt((1.0 - correlation) * (1.0 + correlation));

    for (std::size_t path = 0; path < path_count; ++path) {
        const std::size_t stream_path = antithetic_pairs ? path / 2U : path;
        const double sign = antithetic_pairs && path % 2U != 0U ? -1.0 : 1.0;
        auto first_stream =
            base_stream.substream(static_cast<std::uint64_t>(stream_path));
        auto orthogonal_stream =
            base_stream.substream(static_cast<std::uint64_t>(stream_path)).substream(1U);
        const std::size_t offset = path * interval_count;
        for (std::size_t step = 0; step < interval_count; ++step) {
            const double scale = std::sqrt(grid.interval_years(step));
            const double first_normal = first_stream.standard_normal();
            const double independent_normal = orthogonal_stream.standard_normal();
            first[offset + step] = sign * scale * first_normal;
            second[offset + step] =
                sign * scale *
                (correlation * first_normal + orthogonal_scale * independent_normal);
        }
    }

    return CorrelatedBrownianIncrements{
        BrownianIncrements{grid, path_count, std::move(first)},
        BrownianIncrements{grid, path_count, std::move(second)}};
}

}  // namespace rvre::simulation
