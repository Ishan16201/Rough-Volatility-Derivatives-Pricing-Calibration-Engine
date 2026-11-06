# RVRE architecture

RVRE is a C++20 research library and a set of reproducible synthetic
experiments. Its primary dependency direction is

```text
core types
  -> random streams -> simulation (time grids, Brownian, Volterra)
  -> models (Rough Bergomi, Heston)
  -> pricing (Black--Scholes, Monte Carlo, Heston characteristic function)
  -> calibration, risk sensitivities, hedging experiments
  -> examples and benchmarks
```

Public interfaces live under `include/rvre/`; implementations live under
`src/`. `rvre::rvre` is the CMake library target. Catch2 is fetched at a pinned
revision for tests; the pricing library has no Python runtime dependency.
Python is not required for the CLI or tests. Data interchange is currently
CSV.

## Modules and data flow

* `core/`: unit-distinct asset price, variance, volatility, time, rates,
  discounts and option-price types.
* `random/`: explicit, seedable Philox stream ownership; path and stream IDs
  isolate random sequences.
* `simulation/`: uniform/nonuniform year-fraction grids, Brownian increments,
  correlated Brownian drivers and two Volterra discretizations.
* `models/`: Rough Bergomi path simulation and Heston dynamics/pricing.
  Forward variance, risk-free discounting and dividends are curve interfaces.
* `pricing/`: analytical Black--Scholes prices/Greeks/implied volatility;
  European Monte Carlo pricing and variance reduction.
* `calibration/`: validated option observations, a historical-data adapter
  interface, objectives, constrained multi-start fitting and CSV diagnostics.
* `risk/` and `hedging/`: finite-difference sensitivity tools, discrete
  hedge accounting, and paired synthetic experiments.
* `examples/`, `benchmarks/`, `tests/`: executable workflows, optional
  measurement programs and automated checks.

Calibration accepts normalized option observations. A model pricer evaluates
each candidate parameter vector; objective evaluation, optimizer diagnostics,
residuals and numerical error estimates are retained in a report. A hedge
experiment simulates a market scenario once per seed, then applies each
strategy to those same paths. Strategies in a scenario use a common
scenario-generated initial option mark.

## Reproducibility and failure behavior

Seeds, stream namespaces, path count, time grid, discretization scheme,
optimizer settings and quote weights are explicit inputs. Fixed seed and
configuration replay within the same supported build/runtime. Platform-wide
bitwise equality of floating-point transcendental functions is not promised.
Invalid inputs, failed numerical convergence and nonfinite simulation state
are surfaced as errors or diagnostics; a price with a non-convergence flag is
not represented as a validated result.

Model approximations are selectable/configurable rather than hidden:
cell-average reference and hybrid Volterra simulation; full-truncation
Euler Heston Monte Carlo; characteristic-function quadrature; finite-step
asset integration; and frozen-implied-volatility Rough Bergomi hedge proxy.
Sampling error, discretization error and numerical solver diagnostics are
distinct quantities.

## Build configurations

`CMakePresets.json` supplies `debug`, `release` and
`sanitizers` (AddressSanitizer and UndefinedBehaviorSanitizer where supported).
See [validation](./validation.md) for a clean-checkout reproduction sequence.
The exact CMake and compiler versions are environment-specific and should be
recorded alongside any research result.
