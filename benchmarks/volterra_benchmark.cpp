#include "rvre/random/random_stream.hpp"
#include "rvre/simulation/brownian_motion.hpp"
#include "rvre/simulation/volterra_process.hpp"

#include <charconv>
#include <chrono>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <iostream>
#include <limits>
#include <stdexcept>
#include <string>

namespace {

struct Measurement {
    double elapsed_seconds;
    double terminal_variance;
    double terminal_variance_standard_error;
    double sampled_covariance;
    double covariance_standard_error;
    double terminal_checksum;
    std::size_t estimated_payload_bytes;
};

[[nodiscard]] std::size_t checked_payload_bytes(std::size_t paths,
                                                std::size_t intervals) {
    constexpr std::size_t bytes_per_double = sizeof(double);
    constexpr std::size_t maximum = std::numeric_limits<std::size_t>::max();
    if (intervals > (std::numeric_limits<std::size_t>::max() - 1U) / 2U) {
        throw std::length_error{"Benchmark interval count overflows size_t"};
    }
    const std::size_t per_path_elements = 2U * intervals + 1U;
    if (per_path_elements > maximum / bytes_per_double ||
        paths > maximum / (per_path_elements * bytes_per_double)) {
        throw std::length_error{"Benchmark memory estimate overflows size_t"};
    }
    const std::size_t path_payload =
        paths * per_path_elements * bytes_per_double;
    if (intervals + 1U > maximum / bytes_per_double) {
        throw std::length_error{"Benchmark coefficient memory estimate overflows size_t"};
    }
    const std::size_t coefficient_payload =
        (intervals + 1U) * bytes_per_double;
    if (path_payload > maximum - coefficient_payload) {
        throw std::length_error{"Benchmark memory estimate overflows size_t"};
    }
    return path_payload + coefficient_payload;
}

[[nodiscard]] Measurement measure(rvre::simulation::VolterraScheme scheme,
                                  const rvre::simulation::BrownianIncrements& brownian,
                                  double hurst,
                                  std::size_t paths,
                                  std::size_t intervals) {
    const auto start = std::chrono::steady_clock::now();
    const auto result = scheme == rvre::simulation::VolterraScheme::cell_average_reference
                            ? rvre::simulation::simulate_volterra(brownian, hurst)
                            : rvre::simulation::simulate_volterra(
                                  brownian,
                                  hurst,
                                  scheme,
                                  rvre::random::RandomStream{8801U, 0x485942524944ULL});
    const auto stop = std::chrono::steady_clock::now();

    double mean = 0.0;
    double second_moment = 0.0;
    double checksum = 0.0;
    const std::size_t first_point = intervals / 2U;
    const std::size_t second_point = (3U * intervals) / 4U;
    double first_mean = 0.0;
    double second_mean = 0.0;
    double first_square = 0.0;
    double second_square = 0.0;
    double cross_product = 0.0;
    for (std::size_t path = 0; path < result.path_count(); ++path) {
        const double value = result.value(path, intervals);
        mean += value;
        second_moment += value * value;
        checksum += value;
        const double first = result.value(path, first_point);
        const double second = result.value(path, second_point);
        first_mean += first;
        second_mean += second;
        first_square += first * first;
        second_square += second * second;
        cross_product += first * second;
    }
    const double count = static_cast<double>(paths);
    mean /= count;
    const double sample_variance =
        (second_moment - static_cast<double>(paths) * mean * mean) /
        static_cast<double>(paths - 1U);
    first_mean /= count;
    second_mean /= count;
    const double first_variance =
        first_square / (count - 1.0) - count / (count - 1.0) *
                                             first_mean * first_mean;
    const double second_variance =
        second_square / (count - 1.0) - count / (count - 1.0) *
                                              second_mean * second_mean;
    const double covariance =
        cross_product / (count - 1.0) -
        count / (count - 1.0) * first_mean * second_mean;
    return Measurement{
        .elapsed_seconds =
            std::chrono::duration<double>(stop - start).count(),
        .terminal_variance = sample_variance,
        .terminal_variance_standard_error =
            sample_variance * std::sqrt(2.0 / static_cast<double>(paths - 1U)),
        .sampled_covariance = covariance,
        .covariance_standard_error =
            std::sqrt((first_variance * second_variance + covariance * covariance) /
                      (count - 1.0)),
        .terminal_checksum = checksum,
        .estimated_payload_bytes = checked_payload_bytes(paths, intervals)};
}

[[nodiscard]] std::size_t parse_size(const char* text, const char* name) {
    const std::string value{text};
    std::size_t parsed = 0U;
    const auto [end, error] =
        std::from_chars(value.data(), value.data() + value.size(), parsed);
    if (error != std::errc{} || end != value.data() + value.size() || parsed == 0U) {
        throw std::invalid_argument{std::string{name} + " must be a positive integer"};
    }
    return parsed;
}

}  // namespace

int main(int argc, char* argv[]) {
    if (argc > 4) {
        std::cerr << "Usage: rvre_volterra_benchmark [paths] [time_steps] [H]\n";
        return 2;
    }
    try {
        const std::size_t paths = argc > 1 ? parse_size(argv[1], "paths") : 2000U;
        const std::size_t intervals =
            argc > 2 ? parse_size(argv[2], "time_steps") : 128U;
        std::size_t consumed = 0U;
        const std::string hurst_input = argc > 3 ? argv[3] : "0.1";
        const double hurst = std::stod(hurst_input, &consumed);
        if (consumed != hurst_input.size()) {
            throw std::invalid_argument{"H must be a number"};
        }
        if (!std::isfinite(hurst) || hurst <= 0.0 || hurst >= 0.5) {
            throw std::invalid_argument{"H must be finite and in (0, 0.5)"};
        }
        if (paths < 2U) {
            throw std::invalid_argument{"Benchmark paths must be at least two"};
        }

        const double target_variance = std::pow(1.0, 2.0 * hurst);
        std::cout << "seed=20261106 paths=" << paths << " steps=" << intervals
                  << " H=" << hurst << " grid_maturity_years=1\n"
                  << "method,time_steps,brownian_seconds,elapsed_seconds,terminal_variance,"
                     "absolute_variance_error,variance_sampling_se,"
                     "covariance_t1_t2,absolute_covariance_error,covariance_sampling_se,"
                     "covariance_time_1_years,covariance_time_2_years,"
                     "estimated_payload_bytes,approximate_kernel_terms,checksum\n";
        const std::size_t step_counts[] = {
            intervals / 4U == 0U ? 1U : intervals / 4U,
            intervals / 2U == 0U ? 1U : intervals / 2U,
            intervals};
        std::size_t previous_steps = 0U;
        for (const std::size_t step_count : step_counts) {
            if (step_count == previous_steps) {
                continue;
            }
            previous_steps = step_count;
            const auto grid =
                rvre::simulation::TimeGrid::uniform(rvre::core::Time{1.0}, step_count);
            const auto brownian_start = std::chrono::steady_clock::now();
            const auto brownian = rvre::simulation::generate_brownian_increments(
                grid,
                paths,
                rvre::random::RandomStream{20261106U, 7U});
            const double brownian_seconds =
                std::chrono::duration<double>(std::chrono::steady_clock::now() -
                                              brownian_start)
                    .count();
            for (const auto scheme : {
                     rvre::simulation::VolterraScheme::cell_average_reference,
                     rvre::simulation::VolterraScheme::hybrid_exact_diagonal}) {
                const auto measurement =
                    measure(scheme, brownian, hurst, paths, step_count);
                const char* method =
                    scheme == rvre::simulation::VolterraScheme::cell_average_reference
                        ? "cell_average_reference"
                        : "hybrid_exact_diagonal";
                const double approximate_pair_count =
                    static_cast<double>(paths) *
                    static_cast<double>(step_count) *
                    static_cast<double>(step_count + 1U) * 0.5;
                const double first_time =
                    grid.time_years(step_count / 2U);
                const double second_time =
                    grid.time_years((3U * step_count) / 4U);
                const double target_covariance =
                    rvre::simulation::volterra_theoretical_covariance(
                        hurst,
                        rvre::core::Time{first_time},
                        rvre::core::Time{second_time});
                std::cout << method << ',' << step_count << ','
                          << brownian_seconds << ','
                          << measurement.elapsed_seconds << ','
                          << measurement.terminal_variance << ','
                          << std::abs(measurement.terminal_variance - target_variance)
                          << ',' << measurement.terminal_variance_standard_error << ','
                          << measurement.sampled_covariance << ','
                          << std::abs(measurement.sampled_covariance - target_covariance)
                          << ',' << measurement.covariance_standard_error << ','
                          << first_time << ',' << second_time << ','
                          << measurement.estimated_payload_bytes << ','
                          << approximate_pair_count << ','
                          << measurement.terminal_checksum << '\n';
            }
        }
        std::cout << "Memory is an estimated retained double-payload (Brownian input, "
                     "Volterra output, and coefficient arrays); allocator overhead "
                     "and RNG state are excluded.\n";
        std::cout << "Kernel terms are approximate O(paths*steps^2) operation counts; "
                     "elapsed time includes output allocation and hybrid bridge draws.\n";
        return 0;
    } catch (const std::exception& error) {
        std::cerr << "RVRE Volterra benchmark error: " << error.what() << '\n';
        return 1;
    }
}
