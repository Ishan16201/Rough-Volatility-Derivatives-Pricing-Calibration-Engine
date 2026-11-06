#pragma once

#include "rvre/core/types.hpp"

#include <cstddef>
#include <optional>
#include <vector>

namespace rvre::simulation {

class TimeGrid {
public:
    [[nodiscard]] static TimeGrid uniform(core::Time maturity,
                                          std::size_t interval_count);
    [[nodiscard]] static TimeGrid from_points(std::vector<core::Time> time_points);

    [[nodiscard]] std::size_t interval_count() const noexcept {
        return times_in_years_.size() - 1U;
    }
    [[nodiscard]] std::size_t point_count() const noexcept { return times_in_years_.size(); }
    [[nodiscard]] double time_years(std::size_t point_index) const;
    [[nodiscard]] double interval_years(std::size_t interval_index) const;
    [[nodiscard]] const std::vector<double>& points_years() const noexcept {
        return times_in_years_;
    }
    [[nodiscard]] std::optional<double> uniform_step_years() const noexcept {
        return uniform_step_years_;
    }

private:
    explicit TimeGrid(std::vector<double> time_points_in_years,
                      std::optional<double> uniform_step_years);

    std::vector<double> times_in_years_;
    std::optional<double> uniform_step_years_;
};

}  // namespace rvre::simulation
