#include "rvre/calibration/market_data.hpp"

#include <algorithm>
#include <cmath>
#include <fstream>
#include <stdexcept>
#include <string>
#include <string_view>

namespace rvre::calibration {
namespace {

[[nodiscard]] std::vector<std::string> split_csv_row(const std::string& line) {
    std::vector<std::string> fields;
    std::string field;
    bool quoted = false;
    for (std::size_t i = 0; i < line.size(); ++i) {
        const char character = line[i];
        if (character == '"') {
            if (quoted && i + 1U < line.size() && line[i + 1U] == '"') {
                field.push_back('"');
                ++i;
            } else {
                quoted = !quoted;
            }
        } else if (character == ',' && !quoted) {
            fields.push_back(std::move(field));
            field.clear();
        } else {
            field.push_back(character);
        }
    }
    if (quoted) {
        throw std::invalid_argument{"CSV contains an unterminated quoted field"};
    }
    fields.push_back(std::move(field));
    return fields;
}

[[nodiscard]] double parse_required_double(const std::string& value,
                                          const std::string& field_name) {
    std::size_t consumed = 0U;
    double parsed = 0.0;
    try {
        parsed = std::stod(value, &consumed);
    } catch (const std::exception&) {
        throw std::invalid_argument{"Invalid numeric CSV field: " + field_name};
    }
    if (consumed != value.size() || !std::isfinite(parsed)) {
        throw std::invalid_argument{"Invalid finite numeric CSV field: " +
                                    field_name};
    }
    return parsed;
}

}  // namespace

void validate_option_observations(
    const std::vector<OptionObservation>& observations) {
    if (observations.empty()) {
        throw std::invalid_argument{"Calibration requires at least one observation"};
    }
    for (std::size_t i = 0; i < observations.size(); ++i) {
        const auto& observation = observations[i];
        if ((observation.right != pricing::OptionRight::call &&
             observation.right != pricing::OptionRight::put) ||
            !std::isfinite(observation.strike) || observation.strike <= 0.0 ||
            !std::isfinite(observation.maturity.value) ||
            observation.maturity.value < 0.0 ||
            !std::isfinite(observation.market_price.value) ||
            observation.market_price.value < 0.0 ||
            !std::isfinite(observation.underlying_price.value) ||
            observation.underlying_price.value <= 0.0 ||
            !std::isfinite(observation.risk_free_discount.value) ||
            observation.risk_free_discount.value <= 0.0 ||
            !std::isfinite(observation.dividend_discount.value) ||
            observation.dividend_discount.value <= 0.0 ||
            !std::isfinite(observation.reliability) ||
            observation.reliability < 0.0 || observation.reliability > 1.0) {
            throw std::invalid_argument{
                "Invalid option observation at index " + std::to_string(i)};
        }
        if (observation.bid.has_value() != observation.ask.has_value()) {
            throw std::invalid_argument{
                "Bid and ask must either both be present or both be absent"};
        }
        if (observation.bid &&
            (!std::isfinite(observation.bid->value) ||
             !std::isfinite(observation.ask->value) ||
             observation.bid->value < 0.0 ||
             observation.ask->value < observation.bid->value)) {
            throw std::invalid_argument{
                "Option observation has invalid or crossed bid/ask quotes"};
        }
        if (observation.implied_volatility &&
            (!std::isfinite(observation.implied_volatility->value) ||
             observation.implied_volatility->value < 0.0)) {
            throw std::invalid_argument{
                "Option observation implied volatility must be finite and nonnegative"};
        }
        if (observation.maturity.value == 0.0 &&
            (std::abs(observation.risk_free_discount.value - 1.0) > 1.0e-12 ||
             std::abs(observation.dividend_discount.value - 1.0) > 1.0e-12)) {
            throw std::invalid_argument{
                "Zero-maturity observation discount factors must equal one"};
        }
    }
}

std::vector<OptionObservation> read_option_observations_csv(
    const std::string& path) {
    std::ifstream input{path};
    if (!input) {
        throw std::runtime_error{"Unable to open option observations CSV: " + path};
    }
    std::string header_line;
    if (!std::getline(input, header_line)) {
        throw std::invalid_argument{"Option observations CSV is empty"};
    }
    const auto headers = split_csv_row(header_line);
    const std::vector<std::string> required_headers{
        "right", "strike", "maturity", "market_price", "underlying_price",
        "risk_free_discount", "dividend_discount"};
    for (const auto& required : required_headers) {
        if (std::find(headers.begin(), headers.end(), required) == headers.end()) {
            throw std::invalid_argument{"CSV missing required column: " + required};
        }
    }
    const auto column = [&headers](const std::string& name) {
        const auto found = std::find(headers.begin(), headers.end(), name);
        return static_cast<std::size_t>(std::distance(headers.begin(), found));
    };
    std::vector<OptionObservation> observations;
    std::string line;
    std::size_t line_number = 1U;
    while (std::getline(input, line)) {
        ++line_number;
        if (line.empty()) {
            continue;
        }
        const auto fields = split_csv_row(line);
        if (fields.size() != headers.size()) {
            throw std::invalid_argument{"CSV column count mismatch on line " +
                                        std::to_string(line_number)};
        }
        const auto field = [&fields, line_number](std::size_t index) -> const std::string& {
            if (index >= fields.size()) {
                throw std::invalid_argument{"Invalid CSV column on line " +
                                            std::to_string(line_number)};
            }
            return fields[index];
        };
        const std::string& right_text = field(column("right"));
        pricing::OptionRight right;
        if (right_text == "call") {
            right = pricing::OptionRight::call;
        } else if (right_text == "put") {
            right = pricing::OptionRight::put;
        } else {
            throw std::invalid_argument{"Option right must be call or put on line " +
                                        std::to_string(line_number)};
        }
        const auto optional_column = [&headers, &fields](
                                         const std::string& name)
            -> std::optional<std::string> {
            const auto found = std::find(headers.begin(), headers.end(), name);
            if (found == headers.end()) {
                return std::nullopt;
            }
            const auto index =
                static_cast<std::size_t>(std::distance(headers.begin(), found));
            if (fields[index].empty()) {
                return std::nullopt;
            }
            return fields[index];
        };
        const auto bid = optional_column("bid");
        const auto ask = optional_column("ask");
        const auto implied_volatility = optional_column("implied_volatility");
        const auto reliability = optional_column("reliability");
        observations.push_back(OptionObservation{
            .instrument_id = optional_column("instrument_id").value_or(""),
            .right = right,
            .strike = parse_required_double(field(column("strike")), "strike"),
            .maturity = core::Time{
                parse_required_double(field(column("maturity")), "maturity")},
            .market_price = core::OptionPrice{parse_required_double(
                field(column("market_price")), "market_price")},
            .bid = bid ? std::optional<core::OptionPrice>{core::OptionPrice{
                             parse_required_double(*bid, "bid")}}
                       : std::nullopt,
            .ask = ask ? std::optional<core::OptionPrice>{core::OptionPrice{
                             parse_required_double(*ask, "ask")}}
                       : std::nullopt,
            .implied_volatility =
                implied_volatility
                    ? std::optional<core::Volatility>{core::Volatility{
                          parse_required_double(*implied_volatility,
                                                "implied_volatility")}}
                    : std::nullopt,
            .underlying_price = core::AssetPrice{parse_required_double(
                field(column("underlying_price")), "underlying_price")},
            .risk_free_discount = core::DiscountFactor{parse_required_double(
                field(column("risk_free_discount")), "risk_free_discount")},
            .dividend_discount = core::DiscountFactor{parse_required_double(
                field(column("dividend_discount")), "dividend_discount")},
            .reliability = reliability && !reliability->empty()
                               ? parse_required_double(*reliability, "reliability")
                               : 1.0});
    }
    if (!input.eof()) {
        throw std::runtime_error{"Error while reading option observations CSV: " +
                                 path};
    }
    validate_option_observations(observations);
    return observations;
}

}  // namespace rvre::calibration
