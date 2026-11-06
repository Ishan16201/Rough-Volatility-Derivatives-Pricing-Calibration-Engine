# Benchmarking and performance evidence

Optional CMake targets are enabled with `-DRVRE_BUILD_BENCHMARKS=ON`.
`rvre_volterra_benchmark` reports runtime, discretization variance/covariance
error, sampling standard errors, approximate kernel work and estimated
payload memory across grids. `rvre_pricing_benchmark` reports pricing
estimators, prices, confidence intervals and runtime. `rvre_model_comparison`
reports analytical/characteristic-function prices and Monte Carlo estimates
under declared common synthetic inputs.

Example:

```sh
cmake --preset release -DRVRE_BUILD_BENCHMARKS=ON
cmake --build --preset release --target rvre_volterra_benchmark \
  rvre_pricing_benchmark rvre_model_comparison
./build/release/rvre_volterra_benchmark 1000 64 0.1
./build/release/rvre_pricing_benchmark
./build/release/rvre_model_comparison
```

Calibration, sensitivities and hedging are not currently dedicated
microbenchmark suites; wall-clock smoke-run time is end-to-end CLI time, not
isolated kernel profiling. No measured bottleneck or speedup is claimed from
those runs. Dedicated microbenchmarks, repeated warmups,
hardware/compiler metadata and profiling counters are required before
performance optimization claims. Volterra memory is estimated payload bytes,
not process RSS.

Record compiler/version, CMake version, build mode, CPU/OS, seed, path count,
grid, schemes, model parameters and error diagnostics with any new benchmark.
Compare runtime only alongside numerical error and convergence; deterministic
replay and independent parallel-stream testing remain required for any
parallelization.

## Release smoke measurements

One non-repeated Release run was made on macOS 27.2 (arm64, Apple M5), Apple
Clang 21 and CMake 4.4.2. These are execution checks, not stable performance
estimates:

| Workload | Configuration | Observed time |
|---|---:|---:|
| Brownian generation | 1,000 paths, 64 intervals | 0.00167 s |
| Cell-average Volterra | 1,000 paths, 64 intervals, H=0.1 | 0.00739 s |
| Hybrid Volterra | Same paths/grid/H | 0.00939 s |
| Rough Bergomi pricing, standard | 5,000 paths, 32 steps | 0.01484 s |
| Rough Bergomi pricing, antithetic | Same | 0.01453 s |
| Rough Bergomi pricing, stock control | Same, includes 2,000 pilot paths | 0.02010 s |

At the 64-step Volterra setting the estimated double-payload memory was
1,032,520 bytes. The sample's reference terminal-variance error was 0.2039
(sampling SE 0.0356); hybrid error was 0.0091 (sampling SE 0.0452). The
reference/hybrid covariance sampling diagnostics are also emitted by the
benchmark. Coarse stochastic error in this small benchmark is not a validated
convergence conclusion. A synthetic Black--Scholes calibration with 40 observations and four starts
recovered volatility 0.27 in 13 iterations (reported objective zero; 0.000293
s on this run). A larger hedge benchmark was not performed: the paired hedge
CLI smoke used only 8 paths, 8 steps and one seed (0.38 s), and emitted 36
Heston characteristic-function non-convergence flags. The flag is a
substantive numerical warning, not a
performance result.

The sensitivity test completed 40 assertions in less than the timer's
0.01-second resolution; that measures test overhead as much as the estimator.
There is no isolated calibration/Greek/hedging profiler, and these one-shot
wall times are not substitutes for one. No optimization or speedup is claimed.
