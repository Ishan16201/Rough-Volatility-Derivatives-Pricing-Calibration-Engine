#pragma once

#include "rvre/calibration/market_data.hpp"
#include "rvre/models/heston.hpp"
#include "rvre/pricing/black_scholes.hpp"
#include "rvre/pricing/monte_carlo.hpp"
#include "rvre/simulation/volterra_process.hpp"

#include <cstddef>
#include <cstdint>
#include <optional>
#include <string>
#include <variant>
#include <vector>

namespace rvre::calibration {

struct ParameterControl {
    double initial;
    double lower_bound;
    double upper_bound;
    bool optimize{true};
};

struct BlackScholesCalibrationModel {
    ParameterControl volatility;
};

struct HestonCalibrationModel {
    ParameterControl initial_variance;
    ParameterControl mean_reversion_speed;
    ParameterControl long_run_variance;
    ParameterControl volatility_of_variance;
    ParameterControl correlation;
};

struct ForwardVarianceControl {
    core::Time knot_time;
    ParameterControl variance;
};

struct RoughBergomiCalibrationModel {
    ParameterControl hurst;
    ParameterControl eta;
    ParameterControl correlation;
    std::vector<ForwardVarianceControl> forward_variance_curve;
};

using CalibrationModel =
    std::variant<BlackScholesCalibrationModel, HestonCalibrationModel,
                 RoughBergomiCalibrationModel>;

enum class ObjectiveKind {
    option_price,
    implied_volatility,
    bid_ask_aware
};

enum class WeightingScheme {
    equal,
    reliability,
    inverse_bid_ask_width
};

struct ObjectiveConfig {
    ObjectiveKind kind{ObjectiveKind::option_price};
    WeightingScheme weighting{WeightingScheme::reliability};
    double minimum_reliability{0.0};
    double maximum_relative_bid_ask_width{0.25};
    double minimum_relative_bid_ask_width{1.0e-4};
    double maximum_inverse_width_weight{100.0};
    double rough_variance_curve_smoothness{0.0};
};

struct CalibrationOptimizerConfig {
    std::uint32_t maximum_iterations{1000U};
    std::uint32_t starting_points{4U};
    double objective_tolerance{1.0e-8};
    double parameter_tolerance{1.0e-5};
    double initial_simplex_step{0.08};
    std::uint64_t seed{20261106U};
};

struct CalibrationNumerics {
    std::uint64_t rough_path_count{2048U};
    std::uint64_t rough_batch_size{512U};
    std::uint64_t rough_time_step_count{32U};
    std::uint64_t common_random_seed{731927U};
    std::uint64_t common_random_stream_id{0U};
    double heston_price_error_budget{1.0e-4};
    simulation::VolterraScheme rough_volterra_scheme{
        simulation::VolterraScheme::cell_average_reference};
    models::HestonCharacteristicFunctionConfig heston_characteristic_function{};
};

struct CalibrationRequest {
    CalibrationModel model;
    std::vector<OptionObservation> observations;
    ObjectiveConfig objective;
    CalibrationOptimizerConfig optimizer;
    CalibrationNumerics numerics;
};

enum class CalibrationStatus {
    converged,
    maximum_iterations
};

struct ParameterEstimate {
    std::string name;
    double value;
    double lower_bound;
    double upper_bound;
    bool optimized;
};

struct CalibrationResidual {
    std::string instrument_id;
    pricing::OptionRight right;
    double strike;
    core::Time maturity;
    core::OptionPrice market_price;
    core::OptionPrice model_price;
    core::OptionPrice model_price_standard_error;
    core::OptionPrice model_price_numerical_error;
    std::optional<bool> heston_quadrature_converged;
    double price_residual;
    std::optional<double> market_implied_volatility;
    std::optional<double> model_implied_volatility;
    double objective_residual;
    double applied_weight;
    bool included;
};

struct StartingPointResult {
    std::uint32_t starting_point_index;
    double initial_objective;
    double final_objective;
    std::uint32_t iterations;
    CalibrationStatus status;
};

struct ParameterSensitivity {
    std::string parameter_name;
    double perturbation;
    double objective_up;
    double objective_down;
    double central_objective_derivative;
};

struct CalibrationReport {
    CalibrationStatus status;
    std::vector<ParameterEstimate> parameters;
    double objective_value;
    double initial_objective_value;
    std::uint32_t iterations;
    std::uint64_t objective_evaluations;
    std::uint64_t common_random_seed;
    std::uint64_t common_random_stream_id;
    std::uint64_t validation_random_seed;
    std::uint64_t validation_random_stream_id;
    std::uint64_t rough_path_count;
    std::uint64_t rough_time_step_count;
    double heston_price_error_budget;
    simulation::VolterraScheme rough_volterra_scheme;
    double validation_objective_value;
    double estimated_objective_sampling_error;
    double maximum_model_price_standard_error;
    double maximum_model_price_numerical_error;
    double duration_seconds;
    std::vector<StartingPointResult> starting_points;
    std::vector<CalibrationResidual> residuals;
    std::vector<ParameterSensitivity> sensitivities;
};

[[nodiscard]] CalibrationReport calibrate(const CalibrationRequest& request);

[[nodiscard]] std::vector<OptionObservation> read_option_observations_csv(
    const std::string& path);

void write_calibration_report_csv(const CalibrationReport& report,
                                  const std::string& path);

}  // namespace rvre::calibration
