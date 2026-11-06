# Rough Volatility Research Engine (RVRE)

RVRE is a C++20 research project for stochastic-volatility simulation,
derivatives pricing, calibration, risk and hedging studies. The current module
specifications are collected under [`docs/`](docs/); the initial architecture
proposal is retained in [RVRE_DESIGN.md](RVRE_DESIGN.md).

## Foundation status

The current code provides unit-distinct core types, an explicit Philox RNG,
Brownian/Volterra simulation, Rough Bergomi and Heston models, European option
pricing, implied volatility, synthetic calibration, Greeks and discrete
hedging.
The `rvre_experiment` executable remains a configuration placeholder; use
`rvre_price` for pricing runs.

## Build and test

Requires CMake 3.25+, a C++20 compiler, and Git access to fetch the pinned
Catch2 3.8.1 test dependency on first configure.

```sh
cmake --preset debug
cmake --build --preset debug
ctest --preset debug
./build/debug/rvre_experiment
```

Use the `release` preset for optimized builds and `sanitizers` for
AddressSanitizer/UndefinedBehaviorSanitizer builds on supported Clang/GNU
toolchains. See [docs/FOUNDATION.md](docs/FOUNDATION.md) for the legacy foundation
configuration format and additional details.

The stochastic-process equations, numerical methods, limitations, and
benchmark usage are documented in
[docs/STOCHASTIC_SIMULATION.md](docs/STOCHASTIC_SIMULATION.md).
The Rough Bergomi parameters, curve conventions, discrete variance correction,
asset update, and validation are documented in
[docs/ROUGH_BERGOMI.md](docs/ROUGH_BERGOMI.md).
Pricing estimators, uncertainty reporting, variance reduction, CLI options, and
benchmark conventions are documented in [docs/PRICING.md](docs/PRICING.md).
Black--Scholes analytics, Greeks, implied volatility, Heston characteristic
function and full-truncation Monte Carlo pricing, and cross-model benchmark
assumptions are documented in [docs/CONVENTIONAL_MODELS.md](docs/CONVENTIONAL_MODELS.md).
Calibration data conventions, objectives, optimization limits, uncertainty
diagnostics, and reproducible workflows are documented in
[docs/calibration.md](docs/calibration.md) and its
[detailed reference](docs/calibration_details.md).
Risk-sensitivity estimators, hedge-account conventions, model-misspecification
experiments, and limitations are documented in
[docs/RISK_AND_HEDGING.md](docs/RISK_AND_HEDGING.md).

The research release audit, mathematical conventions, validation, and
benchmarking limitations are summarized in
[docs/research_report.md](docs/research_report.md),
[docs/mathematical_foundations.md](docs/mathematical_foundations.md),
[docs/validation.md](docs/validation.md), and
[docs/benchmarking.md](docs/benchmarking.md).
