#pragma once

#include "rvre/random/random_stream.hpp"
#include "rvre/simulation/time_grid.hpp"

#include <cstddef>
#include <cstdint>
#include <vector>

namespace rvre::simulation {

class BrownianIncrements {
public:
    BrownianIncrements(TimeGrid grid,
                       std::size_t path_count,
                       std::vector<double> increments);

    [[nodiscard]] const TimeGrid& grid() const noexcept { return grid_; }
    [[nodiscard]] std::size_t path_count() const noexcept { return path_count_; }
    [[nodiscard]] std::size_t interval_count() const noexcept {
        return grid_.interval_count();
    }
    [[nodiscard]] double increment(std::size_t path_index,
                                   std::size_t interval_index) const;
    [[nodiscard]] const std::vector<double>& values() const noexcept { return increments_; }

private:
    TimeGrid grid_;
    std::size_t path_count_;
    std::vector<double> increments_;
};

class BrownianPaths {
public:
    BrownianPaths(TimeGrid grid, std::size_t path_count, std::vector<double> values);

    [[nodiscard]] const TimeGrid& grid() const noexcept { return grid_; }
    [[nodiscard]] std::size_t path_count() const noexcept { return path_count_; }
    [[nodiscard]] double value(std::size_t path_index, std::size_t point_index) const;
    [[nodiscard]] const std::vector<double>& values() const noexcept { return values_; }

private:
    TimeGrid grid_;
    std::size_t path_count_;
    std::vector<double> values_;
};

struct CorrelatedBrownianIncrements {
    BrownianIncrements first;
    BrownianIncrements second;
};

[[nodiscard]] BrownianIncrements generate_brownian_increments(
    const TimeGrid& grid,
    std::size_t path_count,
    random::RandomStream base_stream);

[[nodiscard]] BrownianPaths cumulative_brownian_paths(
    const BrownianIncrements& increments);

[[nodiscard]] CorrelatedBrownianIncrements generate_correlated_brownian_increments(
    const TimeGrid& grid,
    std::size_t path_count,
    double correlation,
    random::RandomStream base_stream,
    bool antithetic_pairs = false);

}  // namespace rvre::simulation
