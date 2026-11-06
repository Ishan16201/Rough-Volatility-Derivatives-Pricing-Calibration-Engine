#include "rvre/simulation/time_grid.hpp"

#include <algorithm>
#include <cmath>
#include <limits>
#include <optional>
#include <stdexcept>
#include <utility>

namespace rvre::simulation {

TimeGrid::TimeGrid(std::vector<double> time_points_in_years,
                   std::optional<double> uniform_step_years)
    : times_in_years_{std::move(time_points_in_years)},
      uniform_step_years_{uniform_step_years} {}

TimeGrid TimeGrid::uniform(core::Time maturity, std::size_t interval_count) {
    if (!std::isfinite(maturity.value) || maturity.value <= 0.0) {
        throw std::invalid_argument{"Time-grid maturity must be finite and positive"};
    }
    if (interval_count == 0U) {
        throw std::invalid_argument{"Time grid must contain at least one interval"};
    }
    if (interval_count == std::numeric_limits<std::size_t>::max()) {
        throw std::length_error{"Time-grid point count overflows size_t"};
    }
    std::vector<double> points(interval_count + 1U);
    const double step = maturity.value / static_cast<double>(interval_count);
    for (std::size_t index = 0; index < interval_count; ++index) {
        points[index] = static_cast<double>(index) * step;
    }
    points.back() = maturity.value;
    return TimeGrid{std::move(points), step};
}

TimeGrid TimeGrid::from_points(std::vector<core::Time> time_points) {
    if (time_points.size() < 2U) {
        throw std::invalid_argument{"Time grid requires at least two points"};
    }
    std::vector<double> points;
    points.reserve(time_points.size());
    for (const core::Time time : time_points) {
        points.push_back(time.value);
    }
    if (points.front() != 0.0) {
        throw std::invalid_argument{"Time grid must start at zero years"};
    }
    for (std::size_t index = 0; index < points.size(); ++index) {
        if (!std::isfinite(points[index])) {
            throw std::invalid_argument{"Time-grid points must be finite"};
        }
        if (index > 0U && points[index] <= points[index - 1U]) {
            throw std::invalid_argument{"Time-grid points must be strictly increasing"};
        }
    }
    const double first_step = points[1] - points[0];
    bool is_uniform = true;
    for (std::size_t index = 1U; index + 1U < points.size(); ++index) {
        const double step = points[index + 1U] - points[index];
        const double scale = std::max(std::abs(first_step), std::abs(step));
        if (std::abs(step - first_step) >
            32.0 * std::numeric_limits<double>::epsilon() * scale) {
            is_uniform = false;
            break;
        }
    }
    return TimeGrid{std::move(points),
                    is_uniform ? std::optional<double>{first_step} : std::nullopt};
}

double TimeGrid::time_years(std::size_t point_index) const {
    return times_in_years_.at(point_index);
}

double TimeGrid::interval_years(std::size_t interval_index) const {
    return times_in_years_.at(interval_index + 1U) - times_in_years_.at(interval_index);
}

}  // namespace rvre::simulation
