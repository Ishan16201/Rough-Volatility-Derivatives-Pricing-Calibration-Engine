#include "rvre/config/experiment_config.hpp"
#include "rvre/experiment.hpp"

#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_string.hpp>

#include <chrono>
#include <filesystem>
#include <fstream>
#include <stdexcept>
#include <string>

#ifndef RVRE_SOURCE_DIR
#error RVRE_SOURCE_DIR must be defined for configuration tests
#endif

namespace {

class TemporaryConfig {
public:
    explicit TemporaryConfig(const std::string& contents)
        : path_(std::filesystem::temp_directory_path() /
                ("rvre-invalid-test-config-" +
                 std::to_string(std::chrono::steady_clock::now()
                                    .time_since_epoch()
                                    .count()) +
                 ".rvre")) {
        std::ofstream output{path_};
        output << contents;
        output.close();
        if (!output) {
            throw std::runtime_error{"Unable to create temporary test configuration"};
        }
    }

    ~TemporaryConfig() {
        std::error_code error;
        std::filesystem::remove(path_, error);
    }

    [[nodiscard]] const std::filesystem::path& path() const noexcept { return path_; }

private:
    std::filesystem::path path_;
};

}  // namespace

TEST_CASE("Default configuration loads and placeholder reports no price") {
    const auto path =
        std::filesystem::path{RVRE_SOURCE_DIR} / "config" / "default.rvre";
    const auto config = rvre::config::load_experiment_config(path);
    const auto result = rvre::run_placeholder_experiment(config);

    REQUIRE(config.name == "foundation-smoke-test");
    REQUIRE(config.model.black_scholes.initial_asset.value == 100.0);
    REQUIRE(config.model.black_scholes.volatility.value == 0.2);
    REQUIRE(config.monte_carlo.seed == 20261106U);
    REQUIRE(result.status == rvre::core::SimulationStatus::configured_only);
    REQUIRE_FALSE(result.price_estimate.has_value());
}

TEST_CASE("Configuration rejects unknown keys") {
    const TemporaryConfig config{
        "experiment_name = invalid\n"
        "initial_asset = 100\n"
        "volatility = 0.2\n"
        "risk_free_rate = 0.01\n"
        "dividend_yield = 0\n"
        "maturity_years = 1\n"
        "seed = 1\n"
        "paths = 10\n"
        "time_steps = 2\n"
        "surprise = 7\n"};

    REQUIRE_THROWS_AS(rvre::config::load_experiment_config(config.path()),
                      rvre::config::ConfigError);
}

TEST_CASE("Configuration rejects invalid model and simulation domains") {
    const TemporaryConfig config{
        "experiment_name = invalid\n"
        "initial_asset = 100\n"
        "volatility = -0.2\n"
        "risk_free_rate = 0.01\n"
        "dividend_yield = 0\n"
        "maturity_years = 1\n"
        "seed = 1\n"
        "paths = 10\n"
        "time_steps = 2\n"};

    REQUIRE_THROWS_WITH(
        rvre::config::load_experiment_config(config.path()),
        Catch::Matchers::ContainsSubstring("volatility must be finite and non-negative"));
}

TEST_CASE("Configuration parser reports missing files") {
    const auto path = std::filesystem::temp_directory_path() /
                      "rvre-file-that-does-not-exist.rvre";
    REQUIRE_THROWS_AS(rvre::config::load_experiment_config(path),
                      rvre::config::ConfigError);
}
