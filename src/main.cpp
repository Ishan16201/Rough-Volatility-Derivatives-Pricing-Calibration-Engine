#include "rvre/config/experiment_config.hpp"
#include "rvre/experiment.hpp"

#include <exception>
#include <filesystem>
#include <iostream>

#ifndef RVRE_DEFAULT_CONFIG_PATH
#define RVRE_DEFAULT_CONFIG_PATH "config/default.rvre"
#endif

int main(int argc, char* argv[]) {
    if (argc > 2) {
        std::cerr << "Usage: rvre_experiment [configuration-file]\n";
        return 2;
    }

    const std::filesystem::path config_path =
        argc == 2 ? std::filesystem::path{argv[1]}
                  : std::filesystem::path{RVRE_DEFAULT_CONFIG_PATH};
    try {
        const auto config = rvre::config::load_experiment_config(config_path);
        const auto result = rvre::run_placeholder_experiment(config);
        std::cout << "Experiment: " << config.name << '\n'
                  << "Model: Black-Scholes configuration (not priced)\n"
                  << "Seed: " << result.seed << '\n'
                  << "Paths: " << result.path_count << '\n'
                  << "Time steps: " << result.time_step_count << '\n'
                  << "Status: configuration validated; pricing is not implemented\n";
        return 0;
    } catch (const std::exception& error) {
        std::cerr << "RVRE experiment error: " << error.what() << '\n';
        return 1;
    }
}
