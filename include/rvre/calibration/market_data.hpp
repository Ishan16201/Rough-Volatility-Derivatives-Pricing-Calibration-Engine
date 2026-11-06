#pragma once

#include "rvre/core/types.hpp"
#include "rvre/pricing/option.hpp"

#include <optional>
#include <string>
#include <vector>

namespace rvre::calibration {

struct OptionObservation {
    std::string instrument_id;
    pricing::OptionRight right;
    double strike;
    core::Time maturity;
    core::OptionPrice market_price;
    std::optional<core::OptionPrice> bid;
    std::optional<core::OptionPrice> ask;
    std::optional<core::Volatility> implied_volatility;
    core::AssetPrice underlying_price;
    core::DiscountFactor risk_free_discount;
    core::DiscountFactor dividend_discount;
    double reliability{1.0};
};

struct HistoricalOptionQuery {
    std::string underlying;
    std::string observation_date;
};

class HistoricalOptionDataSource {
public:
    virtual ~HistoricalOptionDataSource() = default;
    [[nodiscard]] virtual std::vector<OptionObservation> load(
        const HistoricalOptionQuery& query) const = 0;
};

void validate_option_observations(
    const std::vector<OptionObservation>& observations);

}  // namespace rvre::calibration
