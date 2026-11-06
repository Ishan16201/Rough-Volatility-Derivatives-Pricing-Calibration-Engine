#include "rvre/experiment.hpp"

#include "rvre/config/experiment_config.hpp"

namespace rvre {

core::SimulationResult run_placeholder_experiment(const core::ExperimentConfig& config) {
    config::validate_experiment_config(config);
    return core::SimulationResult{
        .status = core::SimulationStatus::configured_only,
        .seed = config.monte_carlo.seed,
        .path_count = config.monte_carlo.path_count,
        .time_step_count = config.monte_carlo.time_step_count,
        .price_estimate = std::nullopt};
}

}  // namespace rvre
