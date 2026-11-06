# Validation and reproducibility

## What is checked

Automated tests cover typed domains; Philox known-answer/replay; Brownian
increment moments and correlation; the RL Volterra variance/covariance and
grid behavior; discrete Rough Bergomi variance correction, positivity,
correlation, discounted-forward diagnostics and replay; Black--Scholes
analytic benchmarks and implied-volatility inversion; Heston characteristic
function and full-truncation Monte Carlo; Monte Carlo estimator uncertainty
and variance reduction; synthetic calibration recovery/objectives; Greeks
refinement; and hedge cash accounting and common-mark behavior.

Statistical tests use fixed seeds and tolerances related to expected sampling
variation. They are regression evidence under finite samples, not proofs of
stochastic identities. Model consistency tests are distinguished from market
validation.

On the audit machine (macOS arm64, Apple M5, Apple Clang 21, CMake 4.4.2),
all 53 Catch2 cases passed in Debug (17.71 s), Release (3.66 s) and
AddressSanitizer/UndefinedBehaviorSanitizer (37.44 s). These are environment-
specific test-suite times, not benchmark claims.

## Reproduce from a fresh checkout

Requires CMake 3.25+, a C++20 compiler and Git connectivity on first
configuration to fetch pinned Catch2.

```sh
cmake --preset debug
cmake --build --preset debug
ctest --preset debug --output-on-failure
./build/debug/rvre_rough_bergomi_example /tmp/rvre_paths.csv 128 32
./build/debug/rvre_price --paths 500 --batch-size 100 --steps 8 \
  --seed 20261106 --strikes 100 --maturities 0.5 \
  --methods standard,antithetic,control
./build/debug/rvre_calibrate --model black_scholes --seed 20261106 \
  --starts 4 --output /tmp/rvre_calibration.csv
./build/debug/rvre_hedge --paths 8 --steps 8 --seed 731927 \
  --seed-count 1 --cost 0.001 --output /tmp/rvre_hedging.csv
```

For release and sanitizers:

```sh
cmake --preset release && cmake --build --preset release
ctest --preset release --output-on-failure
cmake --preset sanitizers && cmake --build --preset sanitizers
ctest --preset sanitizers --output-on-failure
```

The default `rvre_experiment` validates the foundational `.rvre` config and is
not a pricing run.

## Audit findings and unresolved validation

The requested RL integral and standard-fBm covariance are mathematically
incompatible for general \(H\); the implementation follows the RL integral
and tests its covariance integral. Hybrid Volterra preserves the diagonal
integral's marginal/covariance but approximates its full joint law.
Characteristic-function tail diagnostics are not rigorous error bounds.
Heston Euler and Rough Bergomi left-point stock updates have discretization
bias; confidence intervals do not capture it.

There is no genuine historical option dataset or vendor adapter in the
repository. Empirical quote validation, calibration, out-of-sample pricing,
and market hedging therefore remain incomplete. No model superiority or
economic benefit is established.
