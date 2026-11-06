#pragma once

#include "rvre/core/types.hpp"
#include "rvre/simulation/time_grid.hpp"
#include "rvre/simulation/volterra_process.hpp"

#include <cstddef>
#include <cstdint>
#include <memory>
#include <vector>

namespace rvre::models {

struct HurstExponent {
    explicit constexpr HurstExponent(double value) noexcept : value{value} {}
    double value;
};

struct VolatilityOfVolatility {
    explicit constexpr VolatilityOfVolatility(double value_per_year_power_h) noexcept
        : value{value_per_year_power_h} {}
    double value;
};

class ForwardVarianceCurve {
public:
    virtual ~ForwardVarianceCurve() = default;
    [[nodiscard]] virtual core::Variance value(core::Time time) const = 0;
};

class DiscountCurve {
public:
    virtual ~DiscountCurve() = default;
    [[nodiscard]] virtual core::DiscountFactor discount_factor(core::Time time) const = 0;
};

class PiecewiseLinearForwardVarianceCurve final : public ForwardVarianceCurve {
public:
    PiecewiseLinearForwardVarianceCurve(std::vector<core::Time> knot_times,
                                        std::vector<core::Variance> variances);
    [[nodiscard]] core::Variance value(core::Time time) const override;

private:
    std::vector<double> knot_times_;
    std::vector<double> variances_;
};

class PiecewiseLinearDiscountCurve final : public DiscountCurve {
public:
    PiecewiseLinearDiscountCurve(std::vector<core::Time> knot_times,
                                 std::vector<core::DiscountFactor> discount_factors);
    [[nodiscard]] core::DiscountFactor discount_factor(core::Time time) const override;

private:
    std::vector<double> knot_times_;
    std::vector<double> discount_factors_;
};

struct RoughBergomiParameters {
    core::AssetPrice initial_asset;
    HurstExponent hurst;
    VolatilityOfVolatility eta;
    double correlation;
    std::shared_ptr<const ForwardVarianceCurve> initial_forward_variance;
    std::shared_ptr<const DiscountCurve> risk_free_discount_curve;
    std::shared_ptr<const DiscountCurve> dividend_discount_curve;
};

struct RoughBergomiConfig {
    std::size_t path_count;
    std::uint64_t seed;
    simulation::VolterraScheme volterra_scheme{
        simulation::VolterraScheme::cell_average_reference};
    std::uint64_t stream_id{0U};
    bool antithetic_pairs{false};
};

struct VariancePathDiagnostic {
    core::Time time;
    core::Variance target_forward_variance;
    core::Variance empirical_mean_variance;
    core::Variance mean_standard_error;
};

struct RoughBergomiDiagnostics {
    double empirical_driver_correlation;
    std::vector<VariancePathDiagnostic> forward_variance_diagnostics;
    double discounted_forward_mean_asset_price;
    double discounted_forward_standard_error;
    double relative_discounted_forward_error;
};

class RoughBergomiPaths {
public:
    RoughBergomiPaths(simulation::TimeGrid grid,
                      simulation::VolterraPaths volterra,
                      std::size_t path_count,
                      std::vector<double> variance_values,
                      std::vector<double> asset_values,
                      RoughBergomiDiagnostics diagnostics);

    [[nodiscard]] const simulation::TimeGrid& grid() const noexcept { return grid_; }
    [[nodiscard]] const simulation::VolterraPaths& volterra() const noexcept {
        return volterra_;
    }
    [[nodiscard]] std::size_t path_count() const noexcept { return path_count_; }
    [[nodiscard]] double variance(std::size_t path_index,
                                  std::size_t point_index) const;
    [[nodiscard]] core::AssetPrice asset_price(std::size_t path_index,
                                               std::size_t point_index) const;
    [[nodiscard]] const std::vector<double>& variance_values() const noexcept {
        return variance_values_;
    }
    [[nodiscard]] const std::vector<double>& asset_values() const noexcept {
        return asset_values_;
    }
    [[nodiscard]] const RoughBergomiDiagnostics& diagnostics() const noexcept {
        return diagnostics_;
    }

private:
    simulation::TimeGrid grid_;
    simulation::VolterraPaths volterra_;
    std::size_t path_count_;
    std::vector<double> variance_values_;
    std::vector<double> asset_values_;
    RoughBergomiDiagnostics diagnostics_;
};

void validate_rough_bergomi_parameters(const RoughBergomiParameters& parameters);

[[nodiscard]] RoughBergomiPaths simulate_rough_bergomi(
    const RoughBergomiParameters& parameters,
    const RoughBergomiConfig& config,
    const simulation::TimeGrid& grid);

}  // namespace rvre::models
