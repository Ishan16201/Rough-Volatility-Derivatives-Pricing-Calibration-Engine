#include "rvre/config/experiment_config.hpp"

#include <array>
#include <charconv>
#include <cmath>
#include <cstdint>
#include <fstream>
#include <locale>
#include <map>
#include <sstream>
#include <string>
#include <string_view>

namespace rvre::config {
namespace {

using Values = std::map<std::string, std::string>;

[[nodiscard]] std::string trim(std::string_view input) {
    constexpr std::string_view whitespace{" \t\r\n"};
    const auto first = input.find_first_not_of(whitespace);
    if (first == std::string_view::npos) {
        return {};
    }
    const auto last = input.find_last_not_of(whitespace);
    return std::string{input.substr(first, last - first + 1U)};
}

[[nodiscard]] Values read_values(const std::filesystem::path& path) {
    std::ifstream input{path};
    if (!input) {
        throw ConfigError{"Unable to open configuration file: " + path.string()};
    }

    Values values;
    std::string line;
    std::size_t line_number = 0;
    while (std::getline(input, line)) {
        ++line_number;
        const auto comment = line.find('#');
        const std::string content = trim(std::string_view{line}.substr(0, comment));
        if (content.empty()) {
            continue;
        }

        const auto separator = content.find('=');
        if (separator == std::string::npos) {
            throw ConfigError{"Expected key=value at line " + std::to_string(line_number)};
        }
        const auto key = trim(std::string_view{content}.substr(0, separator));
        const auto value = trim(std::string_view{content}.substr(separator + 1U));
        if (key.empty() || value.empty()) {
            throw ConfigError{"Empty key or value at line " + std::to_string(line_number)};
        }
        if (!values.emplace(key, value).second) {
            throw ConfigError{"Duplicate key '" + key + "' at line " +
                              std::to_string(line_number)};
        }
    }
    if (input.bad()) {
        throw ConfigError{"Failed while reading configuration file: " + path.string()};
    }
    return values;
}

[[nodiscard]] std::string required_value(const Values& values, const std::string& key) {
    const auto found = values.find(key);
    if (found == values.end()) {
        throw ConfigError{"Missing required configuration key: " + key};
    }
    return found->second;
}

[[nodiscard]] double parse_double(const Values& values, const std::string& key) {
    const std::string text = required_value(values, key);
    std::istringstream parser{text};
    parser.imbue(std::locale::classic());
    double value = 0.0;
    if (!(parser >> value) || !std::isfinite(value)) {
        throw ConfigError{"Configuration key '" + key + "' must be a finite number"};
    }
    parser >> std::ws;
    if (!parser.eof()) {
        throw ConfigError{"Configuration key '" + key + "' must be a finite number"};
    }
    return value;
}

[[nodiscard]] std::uint64_t parse_unsigned(const Values& values, const std::string& key) {
    const std::string text = required_value(values, key);
    std::uint64_t value = 0;
    const auto [end, error] = std::from_chars(text.data(), text.data() + text.size(), value);
    if (error != std::errc{} || end != text.data() + text.size()) {
        throw ConfigError{"Configuration key '" + key + "' must be an unsigned integer"};
    }
    return value;
}

}  // namespace

void validate_experiment_config(const core::ExperimentConfig& config) {
    const auto& model = config.model.black_scholes;
    const auto& monte_carlo = config.monte_carlo;

    if (config.name.empty()) {
        throw ConfigError{"Experiment name must not be empty"};
    }
    if (!std::isfinite(model.initial_asset.value) || model.initial_asset.value <= 0.0) {
        throw ConfigError{"initial_asset must be finite and greater than zero"};
    }
    if (!std::isfinite(model.volatility.value) || model.volatility.value < 0.0) {
        throw ConfigError{"volatility must be finite and non-negative"};
    }
    if (!std::isfinite(model.risk_free_rate.value) ||
        !std::isfinite(model.dividend_yield.value)) {
        throw ConfigError{"interest rates must be finite"};
    }
    if (!std::isfinite(model.maturity.value) || model.maturity.value < 0.0) {
        throw ConfigError{"maturity_years must be finite and non-negative"};
    }
    if (monte_carlo.path_count == 0U) {
        throw ConfigError{"paths must be greater than zero"};
    }
    if (monte_carlo.time_step_count == 0U) {
        throw ConfigError{"time_steps must be greater than zero"};
    }
}

core::ExperimentConfig load_experiment_config(const std::filesystem::path& path) {
    const Values values = read_values(path);
    constexpr std::array<std::string_view, 9> allowed_keys{
        "experiment_name",
        "initial_asset",
        "volatility",
        "risk_free_rate",
        "dividend_yield",
        "maturity_years",
        "seed",
        "paths",
        "time_steps"};

    for (const auto& [key, unused_value] : values) {
        (void)unused_value;
        bool known = false;
        for (const auto allowed_key : allowed_keys) {
            if (key == allowed_key) {
                known = true;
                break;
            }
        }
        if (!known) {
            throw ConfigError{"Unknown configuration key: " + key};
        }
    }

    core::ExperimentConfig config{
        .name = required_value(values, "experiment_name"),
        .model =
            core::ModelParameters{
                .black_scholes =
                    core::BlackScholesParameters{
                        .initial_asset = core::AssetPrice{parse_double(values, "initial_asset")},
                        .volatility = core::Volatility{parse_double(values, "volatility")},
                        .risk_free_rate =
                            core::InterestRate{parse_double(values, "risk_free_rate")},
                        .dividend_yield =
                            core::InterestRate{parse_double(values, "dividend_yield")},
                        .maturity = core::Time{parse_double(values, "maturity_years")}}},
        .monte_carlo =
            core::MonteCarloConfig{
                .seed = parse_unsigned(values, "seed"),
                .path_count = parse_unsigned(values, "paths"),
                .time_step_count = parse_unsigned(values, "time_steps")}};
    validate_experiment_config(config);
    return config;
}

}  // namespace rvre::config
