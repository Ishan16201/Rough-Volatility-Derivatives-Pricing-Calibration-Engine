#include "rvre/random/random_stream.hpp"

#include <cmath>
#include <numbers>

namespace rvre::random {
namespace {

constexpr std::uint32_t philox_multiplier0 = 0xD2511F53U;
constexpr std::uint32_t philox_multiplier1 = 0xCD9E8D57U;
constexpr std::uint32_t philox_weyl0 = 0x9E3779B9U;
constexpr std::uint32_t philox_weyl1 = 0xBB67AE85U;

[[nodiscard]] constexpr std::uint64_t splitmix64(std::uint64_t value) noexcept {
    value += 0x9E3779B97F4A7C15ULL;
    value = (value ^ (value >> 30U)) * 0xBF58476D1CE4E5B9ULL;
    value = (value ^ (value >> 27U)) * 0x94D049BB133111EBULL;
    return value ^ (value >> 31U);
}

[[nodiscard]] constexpr std::uint32_t high_product(std::uint32_t lhs,
                                                   std::uint32_t rhs) noexcept {
    return static_cast<std::uint32_t>((static_cast<std::uint64_t>(lhs) * rhs) >> 32U);
}

[[nodiscard]] constexpr std::uint32_t low_product(std::uint32_t lhs,
                                                  std::uint32_t rhs) noexcept {
    return static_cast<std::uint32_t>(static_cast<std::uint64_t>(lhs) * rhs);
}

}  // namespace

std::array<std::uint32_t, 4> philox4x32_10(
    std::array<std::uint32_t, 4> counter,
    std::array<std::uint32_t, 2> key) noexcept {
    std::uint32_t counter0 = counter[0];
    std::uint32_t counter1 = counter[1];
    std::uint32_t counter2 = counter[2];
    std::uint32_t counter3 = counter[3];
    std::uint32_t key0 = key[0];
    std::uint32_t key1 = key[1];

    for (int round = 0; round < 10; ++round) {
        const std::uint32_t product0_high = high_product(philox_multiplier0, counter0);
        const std::uint32_t product0_low = low_product(philox_multiplier0, counter0);
        const std::uint32_t product1_high = high_product(philox_multiplier1, counter2);
        const std::uint32_t product1_low = low_product(philox_multiplier1, counter2);

        const std::uint32_t next0 = product1_high ^ counter1 ^ key0;
        const std::uint32_t next1 = product1_low;
        const std::uint32_t next2 = product0_high ^ counter3 ^ key1;
        const std::uint32_t next3 = product0_low;
        counter0 = next0;
        counter1 = next1;
        counter2 = next2;
        counter3 = next3;

        if (round != 9) {
            key0 += philox_weyl0;
            key1 += philox_weyl1;
        }
    }

    return {counter0, counter1, counter2, counter3};
}

RandomStream::RandomStream(std::uint64_t seed, std::uint64_t stream_id) noexcept
    : seed_{seed}, stream_id_{stream_id} {}

RandomStream RandomStream::substream(std::uint64_t identifier) const noexcept {
    return RandomStream{seed_, splitmix64(stream_id_ ^ splitmix64(identifier))};
}

void RandomStream::refill() {
    const auto block = philox4x32_10(
        {static_cast<std::uint32_t>(block_counter_),
         static_cast<std::uint32_t>(block_counter_ >> 32U),
         static_cast<std::uint32_t>(stream_id_),
         static_cast<std::uint32_t>(stream_id_ >> 32U)},
        {static_cast<std::uint32_t>(seed_),
         static_cast<std::uint32_t>(seed_ >> 32U)});
    words_[0] = block[0];
    words_[1] = block[1];
    words_[2] = block[2];
    words_[3] = block[3];
    word_index_ = 0U;
    ++block_counter_;
}

std::uint32_t RandomStream::next_word() {
    if (word_index_ == 4U) {
        refill();
    }
    const std::uint32_t word = words_[word_index_];
    ++word_index_;
    return word;
}

double RandomStream::uniform_open() {
    const std::uint64_t upper = static_cast<std::uint64_t>(next_word()) << 20U;
    const std::uint64_t lower = static_cast<std::uint64_t>(next_word() >> 12U);
    const std::uint64_t mantissa = upper | lower;
    return (static_cast<double>(mantissa) + 0.5) * 0x1.0p-52;
}

double RandomStream::standard_normal() {
    if (has_spare_normal_) {
        has_spare_normal_ = false;
        return spare_normal_;
    }

    const double radius = std::sqrt(-2.0 * std::log(uniform_open()));
    const double angle = 2.0 * std::numbers::pi_v<double> * uniform_open();
    spare_normal_ = radius * std::sin(angle);
    has_spare_normal_ = true;
    return radius * std::cos(angle);
}

}  // namespace rvre::random
