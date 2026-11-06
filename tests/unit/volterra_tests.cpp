#include "rvre/simulation/volterra_process.hpp"

#include <catch2/catch_test_macros.hpp>

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <vector>

namespace {

struct EmpiricalMoments {
    double variance;
    double second_variance;
    double covariance;
};

[[nodiscard]] EmpiricalMoments estimate_moments(
    const rvre::simulation::VolterraPaths& paths,
    std::size_t first_point,
    std::size_t second_point) {
    const double count = static_cast<double>(paths.path_count());
    double first_mean = 0.0;
    double second_mean = 0.0;
    double first_square = 0.0;
    double second_square = 0.0;
    double cross_product = 0.0;
    for (std::size_t path = 0; path < paths.path_count(); ++path) {
        const double first = paths.value(path, first_point);
        const double second = paths.value(path, second_point);
        first_mean += first;
        second_mean += second;
        first_square += first * first;
        second_square += second * second;
        cross_product += first * second;
    }
    first_mean /= count;
    second_mean /= count;
    return EmpiricalMoments{
        .variance = first_square / (count - 1.0) -
                    count / (count - 1.0) * first_mean * first_mean,
        .second_variance = second_square / (count - 1.0) -
                           count / (count - 1.0) * second_mean * second_mean,
        .covariance = cross_product / (count - 1.0) -
                      count / (count - 1.0) * first_mean * second_mean};
}

[[nodiscard]] double discrete_covariance(double hurst,
                                         std::size_t steps,
                                         std::size_t first_point,
                                         std::size_t second_point,
                                         rvre::simulation::VolterraScheme scheme) {
    const double dt = 1.0 / static_cast<double>(steps);
    const double exponent = hurst + 0.5;
    double covariance = 0.0;
    for (std::size_t interval = 0;
         interval < std::min(first_point, second_point);
         ++interval) {
        const auto weight = [hurst, exponent, dt, scheme](
                                std::size_t point, std::size_t cell) {
            const std::size_t lag = point - cell;
            if (scheme == rvre::simulation::VolterraScheme::hybrid_exact_diagonal) {
                if (lag == 1U) {
                    return std::sqrt(2.0 * hurst) *
                           std::pow(dt, hurst - 0.5) / exponent;
                }
                return std::sqrt(2.0 * hurst) *
                       std::pow((static_cast<double>(lag) - 0.5) * dt,
                                hurst - 0.5);
            }
            const double lag_value = static_cast<double>(lag);
            const double difference =
                lag == 1U
                    ? 1.0
                    : -std::pow(lag_value, exponent) *
                          std::expm1(exponent * std::log1p(-1.0 / lag_value));
            return std::sqrt(2.0 * hurst) * std::pow(dt, hurst - 0.5) *
                   difference / exponent;
        };
        covariance += weight(first_point, interval) *
                      weight(second_point, interval) * dt;
    }

    if (scheme == rvre::simulation::VolterraScheme::hybrid_exact_diagonal &&
        first_point == second_point) {
        covariance += std::pow(dt, 2.0 * hurst) *
                      std::pow((0.5 - hurst) / exponent, 2.0);
    }
    return covariance;
}

}  // namespace

TEST_CASE("Riemann-Liouville Volterra theoretical variance and covariance are correct") {
    using namespace rvre;
    constexpr double hurst = 0.2;
    const double time = 0.75;
    const double expected_variance = std::pow(time, 2.0 * hurst);
    REQUIRE(simulation::volterra_theoretical_covariance(
                hurst, core::Time{time}, core::Time{time}) ==
            expected_variance);

    const double covariance = simulation::volterra_theoretical_covariance(
        hurst, core::Time{0.5}, core::Time{0.75});
    REQUIRE(covariance > 0.0);
    REQUIRE(covariance < expected_variance);
    REQUIRE(covariance ==
            simulation::volterra_theoretical_covariance(
                hurst, core::Time{0.75}, core::Time{0.5}));
    REQUIRE_THROWS_AS(simulation::volterra_theoretical_covariance(
                          0.5, core::Time{1.0}, core::Time{1.0}),
                      std::invalid_argument);
}

TEST_CASE("Cell-average reference and hybrid Volterra paths match moments within sampling error") {
    using namespace rvre;
    constexpr double hurst = 0.2;
    constexpr std::size_t interval_count = 64U;
    constexpr std::size_t path_count = 30000U;
    constexpr std::size_t first_point = 32U;
    constexpr std::size_t second_point = 48U;
    const auto grid = simulation::TimeGrid::uniform(core::Time{1.0}, interval_count);
    const auto brownian = simulation::generate_brownian_increments(
        grid, path_count, random::RandomStream{564321U, 12U});

    const auto reference = simulation::simulate_volterra(brownian, hurst);
    const auto hybrid = simulation::simulate_volterra(
        brownian,
        hurst,
        simulation::VolterraScheme::hybrid_exact_diagonal,
        random::RandomStream{564321U, 99U});
    const double count = static_cast<double>(path_count);

    for (const auto* paths : {&reference, &hybrid}) {
        REQUIRE(paths->value(0U, 0U) == 0.0);
        const EmpiricalMoments moments =
            estimate_moments(*paths, first_point, second_point);
        const double expected_first_variance =
            discrete_covariance(hurst,
                                interval_count,
                                first_point,
                                first_point,
                                paths->scheme());
        const double expected_second_variance =
            discrete_covariance(hurst,
                                interval_count,
                                second_point,
                                second_point,
                                paths->scheme());
        const double expected_covariance =
            discrete_covariance(hurst,
                                interval_count,
                                first_point,
                                second_point,
                                paths->scheme());
        const double variance_sampling_tolerance =
            6.0 * expected_first_variance * std::sqrt(2.0 / (count - 1.0));
        const double covariance_sampling_tolerance =
            6.0 * std::sqrt(
                      (expected_first_variance * expected_second_variance +
                       expected_covariance * expected_covariance) /
                      (count - 1.0));
        REQUIRE(std::abs(moments.variance - expected_first_variance) <
                variance_sampling_tolerance);
        REQUIRE(std::abs(moments.second_variance - expected_second_variance) <
                6.0 * expected_second_variance *
                    std::sqrt(2.0 / (count - 1.0)));
        REQUIRE(std::abs(moments.covariance - expected_covariance) <
                covariance_sampling_tolerance);
    }
}

TEST_CASE("Volterra discretization errors decrease under time-grid refinement") {
    using namespace rvre;
    const auto deterministic_variance_error =
        [](double hurst, std::size_t steps, simulation::VolterraScheme scheme) {
            const double dt = 1.0 / static_cast<double>(steps);
            double variance = 0.0;
            const double exponent = hurst + 0.5;
            for (std::size_t lag = 1U; lag <= steps; ++lag) {
                const double k = static_cast<double>(lag);
                const double difference =
                    lag == 1U
                        ? 1.0
                        : -std::pow(k, exponent) *
                              std::expm1(exponent * std::log1p(-1.0 / k));
                const double direct_weight =
                    std::sqrt(2.0 * hurst) * std::pow(dt, hurst - 0.5) *
                    difference / exponent;
                double weight = direct_weight;
                if (scheme == simulation::VolterraScheme::hybrid_exact_diagonal) {
                    if (lag == 1U) {
                        weight = std::sqrt(2.0 * hurst) *
                                 std::pow(dt, hurst - 0.5) / exponent;
                        variance += std::pow(dt, 2.0 * hurst) *
                                    std::pow((0.5 - hurst) / exponent, 2.0);
                    } else {
                        weight = std::sqrt(2.0 * hurst) *
                                 std::pow((k - 0.5) * dt, hurst - 0.5);
                    }
                }
                variance += weight * weight * dt;
            }
            return std::abs(variance - 1.0);
        };

    for (const double hurst : {0.1, 0.2, 0.4}) {
        for (const auto scheme : {simulation::VolterraScheme::cell_average_reference,
                                  simulation::VolterraScheme::hybrid_exact_diagonal}) {
            const double coarse_error = deterministic_variance_error(hurst, 16U, scheme);
            const double medium_error = deterministic_variance_error(hurst, 64U, scheme);
            const double fine_error = deterministic_variance_error(hurst, 256U, scheme);
            REQUIRE(medium_error < coarse_error);
            REQUIRE(fine_error < medium_error);
        }

        const double continuous_covariance =
            simulation::volterra_theoretical_covariance(
                hurst, core::Time{0.5}, core::Time{0.75});
        for (const auto scheme : {simulation::VolterraScheme::cell_average_reference,
                                  simulation::VolterraScheme::hybrid_exact_diagonal}) {
            const auto covariance_error = [hurst, scheme, continuous_covariance](
                                              std::size_t steps) {
                const std::size_t first_point = steps / 2U;
                const std::size_t second_point = 3U * steps / 4U;
                const double approximation = discrete_covariance(
                    hurst, steps, first_point, second_point, scheme);
                return std::abs(approximation - continuous_covariance);
            };
            REQUIRE(covariance_error(64U) < covariance_error(16U));
            REQUIRE(covariance_error(256U) < covariance_error(64U));
        }
    }
}

TEST_CASE("Volterra discretizations accept nonuniform time grids") {
    using namespace rvre;
    const auto grid = simulation::TimeGrid::from_points(
        {core::Time{0.0}, core::Time{0.1}, core::Time{0.4}, core::Time{1.0}});
    const auto brownian = simulation::generate_brownian_increments(
        grid, 128U, random::RandomStream{816U});
    const auto paths = simulation::simulate_volterra(brownian, 0.25);
    REQUIRE(paths.grid().point_count() == 4U);
    REQUIRE(paths.value(3U, 0U) == 0.0);
}
