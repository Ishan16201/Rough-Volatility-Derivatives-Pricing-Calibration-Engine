#include "rvre/random/random_stream.hpp"

#include <catch2/catch_test_macros.hpp>

TEST_CASE("Philox4x32-10 matches the Random123 zero-counter known-answer vector") {
    const auto output =
        rvre::random::philox4x32_10({0U, 0U, 0U, 0U}, {0U, 0U});

    REQUIRE(output[0] == 0x6627E8D5U);
    REQUIRE(output[1] == 0xE169C58DU);
    REQUIRE(output[2] == 0xBC57AC4CU);
    REQUIRE(output[3] == 0x9B00DBD8U);
}

TEST_CASE("Random streams replay and substreams are stable and distinct") {
    rvre::random::RandomStream first{7351U, 9U};
    rvre::random::RandomStream replay{7351U, 9U};
    rvre::random::RandomStream independent{7351U, 10U};

    for (int draw = 0; draw < 100; ++draw) {
        REQUIRE(first.standard_normal() == replay.standard_normal());
    }

    auto child_a = rvre::random::RandomStream{7351U, 9U}.substream(42U);
    auto child_b = rvre::random::RandomStream{7351U, 9U}.substream(42U);
    auto child_other = rvre::random::RandomStream{7351U, 9U}.substream(43U);
    REQUIRE(child_a.standard_normal() == child_b.standard_normal());
    REQUIRE(child_a.standard_normal() != child_other.standard_normal());
    REQUIRE(independent.standard_normal() !=
            rvre::random::RandomStream{7351U, 9U}.standard_normal());
}

TEST_CASE("Uniform conversion is strictly between zero and one") {
    rvre::random::RandomStream stream{20261106U};
    for (int draw = 0; draw < 10000; ++draw) {
        const double value = stream.uniform_open();
        REQUIRE(value > 0.0);
        REQUIRE(value < 1.0);
    }
}
