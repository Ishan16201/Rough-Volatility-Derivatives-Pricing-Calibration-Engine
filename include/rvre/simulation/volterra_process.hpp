#pragma once

#include "rvre/core/types.hpp"
#include "rvre/random/random_stream.hpp"
#include "rvre/simulation/brownian_motion.hpp"

#include <cstddef>
#include <vector>

namespace rvre::simulation {

enum class VolterraScheme {
    cell_average_reference,
    hybrid_exact_diagonal
};

class VolterraPaths {
public:
    VolterraPaths(TimeGrid grid,
                  double hurst,
                  VolterraScheme scheme,
                  std::size_t path_count,
                  std::vector<double> values);

    [[nodiscard]] const TimeGrid& grid() const noexcept { return grid_; }
    [[nodiscard]] double hurst() const noexcept { return hurst_; }
    [[nodiscard]] VolterraScheme scheme() const noexcept { return scheme_; }
    [[nodiscard]] std::size_t path_count() const noexcept { return path_count_; }
    [[nodiscard]] double value(std::size_t path_index, std::size_t point_index) const;
    [[nodiscard]] const std::vector<double>& values() const noexcept { return values_; }

private:
    TimeGrid grid_;
    double hurst_;
    VolterraScheme scheme_;
    std::size_t path_count_;
    std::vector<double> values_;
};

[[nodiscard]] VolterraPaths simulate_volterra(
    const BrownianIncrements& brownian,
    double hurst);

[[nodiscard]] VolterraPaths simulate_volterra(
    const BrownianIncrements& brownian,
    double hurst,
    VolterraScheme scheme,
    random::RandomStream bridge_random,
    bool antithetic_pairs = false);

[[nodiscard]] double volterra_theoretical_covariance(double hurst,
                                                      core::Time first_time,
                                                      core::Time second_time);

[[nodiscard]] double volterra_discrete_variance(const TimeGrid& grid,
                                                 double hurst,
                                                 VolterraScheme scheme,
                                                 std::size_t point_index);

}  // namespace rvre::simulation
