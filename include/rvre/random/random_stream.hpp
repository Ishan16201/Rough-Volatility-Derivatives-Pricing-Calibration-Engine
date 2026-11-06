#pragma once

#include <array>
#include <cstdint>

namespace rvre::random {

[[nodiscard]] std::array<std::uint32_t, 4> philox4x32_10(
    std::array<std::uint32_t, 4> counter,
    std::array<std::uint32_t, 2> key) noexcept;

class RandomStream {
public:
    explicit RandomStream(std::uint64_t seed, std::uint64_t stream_id = 0U) noexcept;

    [[nodiscard]] std::uint64_t seed() const noexcept { return seed_; }
    [[nodiscard]] std::uint64_t stream_id() const noexcept { return stream_id_; }

    [[nodiscard]] RandomStream substream(std::uint64_t identifier) const noexcept;
    [[nodiscard]] double uniform_open();
    [[nodiscard]] double standard_normal();

private:
    void refill();
    [[nodiscard]] std::uint32_t next_word();

    std::uint64_t seed_;
    std::uint64_t stream_id_;
    std::uint64_t block_counter_{0U};
    std::uint32_t words_[4]{};
    std::uint32_t word_index_{4U};
    bool has_spare_normal_{false};
    double spare_normal_{0.0};
};

}  // namespace rvre::random
