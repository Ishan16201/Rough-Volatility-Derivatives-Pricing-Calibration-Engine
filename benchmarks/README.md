# Benchmarks

With `-DRVRE_BUILD_BENCHMARKS=ON`, RVRE provides:

* `rvre_volterra_benchmark`: Brownian generation and Volterra reference/hybrid
  runtime, grid-refinement variance/covariance errors, sampling uncertainty
  and estimated retained memory.
* `rvre_pricing_benchmark`: standard, antithetic and stock-control
  Rough Bergomi pricing estimates, confidence intervals and runtime.
* `rvre_model_comparison`: Black--Scholes, Heston and Rough Bergomi synthetic
  benchmark prices with sampling/quadrature diagnostics.

Usage and interpretation are documented in
[benchmarking](../docs/benchmarking.md),
[stochastic simulation](../docs/STOCHASTIC_SIMULATION.md), and
[conventional models](../docs/CONVENTIONAL_MODELS.md).

There are no dedicated microbenchmark targets yet for calibration, Greeks or
hedging. Future benchmarks must report numerical error and convergence
alongside runtime, and record compiler, build mode, hardware, dependencies,
and experiment settings.
