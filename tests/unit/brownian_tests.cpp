#include "rvre/simulation/brownian_motion.hpp"

#include <catch2/catch_test_macros.hpp>

#include <cmath>
#include <vector>

TEST_CASE("Time grids validate zero origin and increasing finite times") {
    using rvre::core::Time;
    using rvre::simulation::TimeGrid;

    const auto uniform = TimeGrid::uniform(Time{1.0}, 4U);
    REQUIRE(uniform.point_count() == 5U);
    REQUIRE(uniform.interval_count() == 4U);
    REQUIRE(uniform.uniform_step_years().has_value());
    REQUIRE(uniform.interval_years(2U) == 0.25);

    const auto irregular =
        TimeGrid::from_points({Time{0.0}, Time{0.1}, Time{0.35}, Time{1.0}});
    REQUIRE_FALSE(irregular.uniform_step_years().has_value());
    REQUIRE(irregular.interval_years(2U) == 0.65);

    REQUIRE_THROWS_AS(TimeGrid::from_points({Time{0.1}, Time{1.0}}),
                      std::invalid_argument);
    REQUIRE_THROWS_AS(TimeGrid::from_points({Time{0.0}, Time{0.5}, Time{0.5}}),
                      std::invalid_argument);
    REQUIRE_THROWS_AS(TimeGrid::uniform(Time{1.0}, 0U), std::invalid_argument);
}

TEST_CASE("Brownian increments have the configured moments and independent steps") {
    using namespace rvre;
    constexpr std::size_t paths = 30000U;
    const auto grid = simulation::TimeGrid::from_points(
        {core::Time{0.0}, core::Time{0.2}, core::Time{0.7}});
    const auto increments = simulation::generate_brownian_increments(
        grid, paths, random::RandomStream{417U});

    for (std::size_t step = 0; step < grid.interval_count(); ++step) {
        const double expected_variance = grid.interval_years(step);
        double mean = 0.0;
        double second_moment = 0.0;
        for (std::size_t path = 0; path < paths; ++path) {
            const double value = increments.increment(path, step);
            mean += value;
            second_moment += value * value;
        }
        mean /= static_cast<double>(paths);
        const double variance =
            second_moment / static_cast<double>(paths - 1U) -
            static_cast<double>(paths) /
                static_cast<double>(paths - 1U) * mean * mean;
        const double mean_tolerance =
            6.0 * std::sqrt(expected_variance / static_cast<double>(paths));
        const double variance_tolerance =
            6.0 * expected_variance *
            std::sqrt(2.0 / static_cast<double>(paths - 1U));
        REQUIRE(std::abs(mean) < mean_tolerance);
        REQUIRE(std::abs(variance - expected_variance) < variance_tolerance);
    }

    double cross_product_mean = 0.0;
    double cross_path_product_mean = 0.0;
    std::size_t independent_path_pairs = 0U;
    for (std::size_t path = 0; path < paths; ++path) {
        cross_product_mean += increments.increment(path, 0U) *
                              increments.increment(path, 1U);
        if (path % 2U == 0U && path + 1U < paths) {
            cross_path_product_mean += increments.increment(path, 0U) *
                                       increments.increment(path + 1U, 0U);
            ++independent_path_pairs;
        }
    }
    cross_product_mean /= static_cast<double>(paths);
    const double covariance_standard_error =
        std::sqrt(grid.interval_years(0U) * grid.interval_years(1U) /
                  static_cast<double>(paths));
    REQUIRE(std::abs(cross_product_mean) < 6.0 * covariance_standard_error);
    cross_path_product_mean /= static_cast<double>(independent_path_pairs);
    const double cross_path_standard_error =
        grid.interval_years(0U) /
        std::sqrt(static_cast<double>(independent_path_pairs));
    REQUIRE(std::abs(cross_path_product_mean) < 6.0 * cross_path_standard_error);
}

TEST_CASE("Brownian paths are cumulative and begin at zero") {
    using namespace rvre;
    const auto grid = simulation::TimeGrid::uniform(core::Time{1.0}, 3U);
    const auto increments = simulation::generate_brownian_increments(
        grid, 32U, random::RandomStream{891U});
    const auto paths = simulation::cumulative_brownian_paths(increments);

    for (std::size_t path = 0; path < paths.path_count(); ++path) {
        REQUIRE(paths.value(path, 0U) == 0.0);
        for (std::size_t step = 0; step < grid.interval_count(); ++step) {
            REQUIRE(paths.value(path, step + 1U) ==
                    paths.value(path, step) + increments.increment(path, step));
        }
    }
}

TEST_CASE("Correlated Brownian increments preserve marginal variance and correlation") {
    using namespace rvre;
    constexpr std::size_t paths = 40000U;
    constexpr double target_correlation = -0.65;
    const auto grid = simulation::TimeGrid::uniform(core::Time{1.0}, 2U);
    const auto pair = simulation::generate_correlated_brownian_increments(
        grid, paths, target_correlation, random::RandomStream{19283U, 4U});

    double first_mean = 0.0;
    double second_mean = 0.0;
    double first_square = 0.0;
    double second_square = 0.0;
    double cross_product = 0.0;
    for (std::size_t path = 0; path < paths; ++path) {
        const double first = pair.first.increment(path, 0U);
        const double second = pair.second.increment(path, 0U);
        first_mean += first;
        second_mean += second;
        first_square += first * first;
        second_square += second * second;
        cross_product += first * second;
    }
    const double count = static_cast<double>(paths);
    first_mean /= count;
    second_mean /= count;
    const double first_variance = first_square / (count - 1.0) -
                                  count / (count - 1.0) * first_mean * first_mean;
    const double second_variance = second_square / (count - 1.0) -
                                   count / (count - 1.0) * second_mean * second_mean;
    const double covariance = cross_product / (count - 1.0) -
                              count / (count - 1.0) * first_mean * second_mean;
    const double correlation =
        covariance / std::sqrt(first_variance * second_variance);

    REQUIRE(std::abs(first_mean) < 6.0 / std::sqrt(count * 2.0));
    REQUIRE(std::abs(second_mean) < 6.0 / std::sqrt(count * 2.0));
    REQUIRE(std::abs(first_variance - 0.5) < 6.0 * std::sqrt(2.0 / (count - 1.0)) * 0.5);
    REQUIRE(std::abs(second_variance - 0.5) < 6.0 * std::sqrt(2.0 / (count - 1.0)) * 0.5);
    REQUIRE(std::abs(correlation - target_correlation) < 6.0 / std::sqrt(count));

    REQUIRE_THROWS_AS(
        simulation::generate_correlated_brownian_increments(
            grid, 1U, 1.01, random::RandomStream{1U}),
        std::invalid_argument);
}
