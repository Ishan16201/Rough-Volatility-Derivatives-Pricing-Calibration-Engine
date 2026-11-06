#include "rvre/core/types.hpp"

#include <catch2/catch_test_macros.hpp>
#include <type_traits>

TEST_CASE("Financial quantities use distinct explicit types") {
    STATIC_REQUIRE(!std::is_convertible_v<double, rvre::core::AssetPrice>);
    STATIC_REQUIRE(!std::is_convertible_v<rvre::core::AssetPrice, rvre::core::Variance>);
    STATIC_REQUIRE(!std::is_convertible_v<rvre::core::Time, rvre::core::InterestRate>);

    const rvre::core::AssetPrice spot{100.0};
    const rvre::core::Volatility volatility{0.2};
    const rvre::core::Time maturity{1.0};
    REQUIRE(spot.value == 100.0);
    REQUIRE(volatility.value == 0.2);
    REQUIRE(maturity.value == 1.0);
}

TEST_CASE("Simulation results distinguish an unrun experiment from a price") {
    const rvre::core::SimulationResult result{
        .status = rvre::core::SimulationStatus::configured_only,
        .seed = 123U,
        .path_count = 64U,
        .time_step_count = 12U,
        .price_estimate = std::nullopt};

    REQUIRE(result.status == rvre::core::SimulationStatus::configured_only);
    REQUIRE_FALSE(result.price_estimate.has_value());
}
