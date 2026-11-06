#pragma once

#include "rvre/core/types.hpp"

#include <filesystem>
#include <stdexcept>

namespace rvre::config {

class ConfigError : public std::runtime_error {
public:
    using std::runtime_error::runtime_error;
};

[[nodiscard]] core::ExperimentConfig load_experiment_config(
    const std::filesystem::path& path);

void validate_experiment_config(const core::ExperimentConfig& config);

}  // namespace rvre::config
