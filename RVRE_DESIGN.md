# Rough Volatility Research Engine (RVRE)

## 1. Purpose and research standard

RVRE is a C++20 research platform for simulating stochastic-volatility models,
pricing derivatives, calibrating to supplied option observations, and testing
hedging hypotheses. Its initial focus is a controlled comparison of Black–Scholes,
Heston, and rough Bergomi. Rough Heston is a later extension, not part of the
initial implementation.

The primary question is whether rough volatility improves pricing and hedging
relative to conventional stochastic volatility, after accounting for sampling
uncertainty, discretization error, calibration choices, and model complexity.
The engine must report what was measured, on what data, with what numerical
settings. It must not present synthetic fixtures as market data or claim
empirical evidence until genuine, documented market observations have been used.

This document is a design specification only. No engine implementation or
empirical result is included.

## 2. Mathematical conventions and model specifications

### 2.1 Common conventions

- Pricing is under a risk-neutral measure \(Q\). Initially, \(r\) and \(q\)
  are constant continuously compounded risk-free and dividend yields. The
  stock drift is therefore \(r-q\); the martingale equation in the project
  brief is the diffusion component.
- Time is measured in years. Rates and instantaneous variances have units
  year\(^{-1}\). Brownian increments have units year\(^{1/2}\).
- \(S_t\) is a positive asset price, \(K>0\) a strike, and \(T\geq0\) maturity.
  Monetary values use the currency/unit of \(S\). Model parameters are typed
  quantities where practical, rather than unlabelled generic doubles.
- The initial implementation prices European vanilla options. Early exercise,
  jumps, stochastic rates, and transaction costs are out of scope unless an
  experiment explicitly adds them.
- Curves are deterministic inputs with documented interpolation and extrapolation.
  The initial constant-rate / constant-curve cases remain supported.

### 2.2 Black–Scholes

Under \(Q\):

\[
dS_t/S_t=(r-q)\,dt+\sigma\,dW_t,\qquad \sigma\geq0.
\]

Inputs satisfy \(S_0>0\), \(K>0\), \(T\geq0\), and finite \(r,q\). Volatility
\(\sigma\) has units year\(^{-1/2}\). The analytic European price is the
reference implementation; a log-Euler simulator is useful for generic Monte
Carlo interfaces but must not replace the analytic benchmark.

### 2.3 Heston

\[
\begin{aligned}
dS_t/S_t &= (r-q)\,dt+\sqrt{v_t}\,dW_t^S,\\
dv_t &= \kappa(\theta-v_t)\,dt+\xi\sqrt{v_t}\,dW_t^v,\\
d\langle W^S,W^v\rangle_t &= \rho\,dt.
\end{aligned}
\]

Use \(v_0\geq0,\ \kappa>0,\ \theta>0,\ \xi\geq0,\ \rho\in[-1,1]\);
\(v,\theta,v_0\) have units year\(^{-1}\), \(\kappa\) year\(^{-1}\), and
\(\xi\) year\(^{-1/2}\). The Feller condition \(2\kappa\theta\geq\xi^2\)
is recorded as a useful boundary-attainment condition, not silently imposed as
a universal calibration constraint.

The initial simulation scheme is full-truncation Euler. With
\(v_i^+=\max(v_i,0)\), update the raw variance state as
\[
v_{i+1}=v_i+\kappa(\theta-v_i^+)\Delta t+
\xi\sqrt{v_i^+}\Delta W_i^v.
\]
The raw state may be negative; only its positive part is used for the next
drift/diffusion and for the stock step
\[
\log(S_{i+1}/S_i)=(r-q-\tfrac12v_i^+)\Delta t+
\sqrt{v_i^+}\Delta W_i^S.
\]
Use \(\Delta W_i^S=\rho\Delta W_i^v+
\sqrt{1-\rho^2}\Delta W_i^\perp\). This scheme is not claimed to be
unbiased or exact. Add a quadratic-exponential scheme only as a separately
validated alternative. State representation and non-negative variance used
by the stock process must not be conflated.
Analytic Heston European prices (where the selected characteristic-function
implementation is validated) provide an independent pricing benchmark.

### 2.4 Rough Bergomi

The initial model uses the specified Volterra Gaussian driver:

\[
Y_t=\sqrt{2H}\int_0^t(t-s)^{H-\frac12}\,dW_s,\qquad
V_t=\xi_0(t)\exp\left(\eta Y_t-\tfrac12\eta^2t^{2H}\right),
\]
\[
dS_t/S_t=(r-q)\,dt+\sqrt{V_t}\,dW_t^S,\qquad
d\langle W^S,W\rangle_t=\rho\,dt.
\]

Here \(H\in(0,\frac12)\), \(\eta\geq0\), \(\rho\in[-1,1]\), and
\(\xi_0(t)>0\) is a deterministic initial forward-variance curve. \(\xi_0\)
and \(V_t\) have units year\(^{-1}\); \(Y_t\) has units year\(^{H}\);
\(\eta\) has units year\(^{-H}\). Thus the exponent is dimensionless.
The normalization gives
\(\operatorname{Var}(Y_t)=2H\int_0^t(t-s)^{2H-1}ds=t^{2H}\), so
\(\mathbb E[V_t]=\xi_0(t)\) for deterministic \(\xi_0\). This is the reason
for the \(-\frac12\eta^2t^{2H}\) correction; it is not an extra fitted drift.
At \(t=0\), \(Y_0=0\) and \(V_0=\xi_0(0)\).

For increments \(\Delta W_j=W_{t_{j+1}}-W_{t_j}\), correlation is constructed
without changing the Volterra driver:

\[
\Delta W^S_j=\rho\,\Delta W_j+\sqrt{1-\rho^2}\,\Delta W^\perp_j,
\]

where \(W^\perp\) is independent. This yields the requested instantaneous
covariation. At \(\rho=\pm1\), the independent component vanishes.

### 2.5 Discretization contract

The simulation grid is \(t_i=i\Delta t\), \(i=0,\ldots,N\), with
\(\Delta t=T/N\), unless a caller supplies a documented nonuniform grid.
Brownian increments are independent \(N(0,\Delta t)\) variates. The initial
stock update for stochastic-volatility models is the left-point log step

\[
\log(S_{i+1}/S_i)=(r-q-\tfrac12V_i)\Delta t+
\sqrt{V_i}\Delta W^S_i.
\]

This is a time-discretization approximation; reporting a Monte Carlo confidence
interval does not account for its bias. Grid refinement is therefore a
separate requirement.

For rough Bergomi, distinguish two methods:

1. **Direct convolution (baseline).** Approximate the stochastic convolution
   using cell-averaged kernel weights and Brownian increments. For a uniform
   grid, one usable weight convention is
   \[
   a_{i,j}=\frac{\sqrt{2H}}{\Delta t(H+\frac12)}
   \left[(t_i-t_j)^{H+\frac12}-(t_i-t_{j+1})^{H+\frac12}\right],
   \quad j<i,
   \qquad Y_{t_i}\approx\sum_{j<i}a_{i,j}\Delta W_j.
   \]
   This is the average of the kernel over each cell, not an exact integral
   against Brownian motion. It costs \(O(N^2)\) operations per path with a
   straightforward implementation. The singularity at the most recent cell
   makes discretization studies essential, especially at small \(H\).
2. **Hybrid scheme (validated alternative).** Treat the near-diagonal,
   singular-kernel contribution with a scheme that preserves its Gaussian
   variance and covariance with the interval Brownian increment, and use
   coarser kernel approximations away from the diagonal. The implementation
   must document its exact coefficient construction and cite its method
   provenance in code/docs. It is not accepted merely because it is faster:
   compare its covariance and prices against direct convolution and a
   fine-grid reference over the declared \(H\) range.

An exact-grid Gaussian reference can be formed from the covariance
\[
C(t,s)=2H\int_0^{\min(t,s)}(t-u)^{H-\frac12}(s-u)^{H-\frac12}\,du
\]
and a stable factorization. It is expensive (dense storage/factorization) and
is for small-grid verification, not production paths. Optional Markovian
sum-of-exponentials approximations require an explicit approximation interval,
error tolerance, and convergence study; they are not the baseline and must
not be described as exact.

## 3. Numerical and reproducibility policy

### 3.1 Randomness and experiment identity

- A `RandomStream` uses counter-based Philox4x32-10 with fixed-width,
  versioned seed/key and path/time counter mapping, plus a project-defined
  normal transform. Include published/reference known-answer tests. Do not
  assume that `std::normal_distribution` produces identical values across
  standard-library vendors.
- Derive independent streams by stable experiment/path identifiers, not by
  thread scheduling. Record seed, RNG version, path count, grid, model
  parameters, compiler/build mode, and relevant dependency versions.
- Promise deterministic replay for a fixed supported build/runtime. Treat
  cross-platform bitwise equality as a separately tested capability; floating
  point and transcendental implementations can differ.
- Store experiment configuration and results in a machine-readable format
  alongside a human-readable summary. Do not overwrite inputs or results.

### 3.2 Error reporting

Every Monte Carlo price includes sample count, sample mean, standard error,
confidence interval and, where applicable, the variance-reduction method and
its estimated variance reduction. Confidence intervals quantify sampling
uncertainty, not model error or time-discretization bias. Report grid and path
refinement independently. Numerical failures (non-finite values, invalid
inputs, failed root finding, covariance/factorization failure) are explicit
typed errors with context; do not silently clamp or return a plausible price.

Use double precision by default. Any clipping, truncation, fallback scheme, or
regularization must be an explicit configuration recorded with results.

### 3.3 Variance reduction and Greeks

Implement antithetic variates and a control variate with a known expectation
(initially the Black–Scholes price under a clearly stated matching setup).
Estimate control coefficients from an independent pilot or use a documented
estimator that accounts for coefficient estimation. Common random numbers
support paired model comparisons and finite-difference Greeks.

Greeks expose method and error limitations. Initially use pathwise estimators
where payoff/model differentiability permits and common-random-number finite
differences as a general fallback. Likelihood-ratio estimators may be added
for discontinuous payoffs. Delta hedging uses documented re-hedge dates and
self-financing cash-account conventions; results include hedge P&L summaries
and uncertainty, not only a single path.

### 3.4 Implied volatility

Invert European Black–Scholes prices with a bracketed method (Brent or
bisection as robust baseline); Newton steps may accelerate only while staying
inside a valid bracket. Enforce and report no-arbitrage price bounds for the
given \(S_0,K,T,r,q\). Handle \(T=0\), zero volatility, intrinsic-boundary
prices, and out-of-bounds inputs explicitly. Return convergence status,
iterations, final price residual, and bracket; no silent zero-volatility
fallback.

### 3.5 Calibration

Calibration accepts an explicit set of option observations with fields for
instrument, strike, maturity, bid/ask or selected target, quote timestamp,
source, and quality flags. The initial objective is a configurable weighted
least-squares residual, with weights and price/volatility residual choice
recorded. Do not treat bid/ask midpoints as error-free truth.

Enforce parameter domains by documented transforms or constrained optimization.
Use multiple declared starting points, report convergence diagnostics,
residuals by maturity/strike, parameter-bound hits, and sensitivity to starts
and weighting. Separate optimizer convergence from identifiability and from
out-of-sample performance. Synthetic-surface fitting is an integration test
and calibration-recovery experiment, not evidence of market fit.

## 4. Architecture and repository layout

Keep the public API narrow, value-oriented, and independent of a command-line
front end. Suggested initial tree:

```text
/
  CMakeLists.txt
  CMakePresets.json
  README.md
  RVRE_DESIGN.md
  cmake/
  include/rvre/
    core/              # dates/time grids, units, errors, result types
    random/            # RNG interface, streams, normal variates
    simulation/        # Brownian and correlated Brownian paths
    rough/             # Volterra process and rough Bergomi
    models/            # Black–Scholes and Heston
    pricing/            # payoffs, Monte Carlo, variance reduction
    numerics/           # root finding, quadrature, linear algebra adapters
    calibration/        # quote sets, objectives, optimizers/adapters
    risk/               # Greeks and delta hedging
    data/               # market-data interfaces and fixture formats
    analytics/          # diagnostics and research statistics
    benchmark/          # benchmark inputs and reporting interfaces
  src/                  # implementation matching public headers
  apps/                 # experiment/benchmark CLI (later milestone)
  tests/
    unit/
    integration/
    regression/
    fixtures/           # deterministic synthetic fixtures, labeled synthetic
  experiments/          # versioned configurations; no untracked claims/data
```

Core types include `TimeGrid`, `BrownianPath`, `CorrelatedBrownianPaths`,
`OptionContract`, `MarketQuote`, `ModelParameters`, `PriceEstimate`,
`CalibrationReport`, `GreekEstimate`, and a typed `Result<T, Error>`. Define
ownership, units, and valid domains at construction boundaries. Avoid exposing
Eigen types as a long-term public ABI if an internal linear algebra dependency
is used.

### Module responsibilities

| Module | Responsibility / boundary |
|---|---|
| Core mathematical types | Units, time grids, contracts, validated parameters, error/result and estimate types. No model-specific simulation. |
| RNG | Seeded, versioned streams and independent substreams; no hidden global generator. |
| Brownian motion | Generate increments and paths on a supplied grid; deterministic replay. |
| Correlated Brownian motion | Correlate streams by explicit matrix/factorization; validate positive semidefiniteness and report failures. |
| Volterra Gaussian process | Build rough driver from increments using named discretization; expose method metadata and small-grid covariance checks. |
| Rough Bergomi | Map driver to variance, correlate stock Brownian motion, evolve stock, expose the forward-variance curve. |
| Black–Scholes | Analytic price/Greeks and optional path simulator serving as reference/control variate. |
| Heston | Parameter validation, selected simulation schemes, and validated analytic price reference. |
| Monte Carlo pricing | Payoffs, path aggregation, uncertainty estimates, batching and deterministic parallel stream assignment. |
| Variance reduction | Antithetic and control variates with estimator diagnostics; composable with pricing. |
| Implied volatility | Bounds, robust inversion, residual and convergence reporting. |
| Calibration | Data-independent objective/optimizer interface and model-specific parameter adapters; report diagnostics. |
| Greeks / hedging | Greek estimators and self-financing discrete hedge accounting with uncertainty. |
| Market data ingestion | Import adapter boundary, schema, provenance, normalization, validation; no network dependency in core. |
| Research analytics | Paired comparisons, confidence intervals, residual and hedge-P&L diagnostics; no causal claims from summary metrics alone. |
| Benchmarking | Reproducible workloads and accuracy-vs-cost reports; compare numerical error as well as runtime. |

### Dependency graph

```text
core
 ├── random
 ├── numerics
 ├── data
 └── simulation
      ├── brownian
      ├── correlated_brownian
      └── volterra
           └── rough_bergomi
core + simulation ──> black_scholes / heston / rough_bergomi
models + pricing ──> variance_reduction
models + pricing ──> implied_volatility
data + models + pricing + numerics ──> calibration
models + simulation + pricing ──> risk / hedging
pricing + calibration + risk + data ──> research_analytics
all measured components ──> benchmark / experiment applications
```

Dependencies point inward to core abstractions. Core must not depend on models,
data vendors, command-line code, or a particular optimizer. A market-data
adapter can depend on the data schema, never the reverse.

## 5. Build and dependency decisions

- Require CMake 3.25 or newer and C++20 (the preset schema version requires
  3.25); enable compiler warnings and
  `BUILD_TESTING` through CTest. Export a library target and keep applications
  optional.
- Prefer standard-library facilities for core types, filesystem, and
  threading. Use Eigen 3 for internal dense linear algebra/factorization if
  needed; pin or constrain the supported version and hide it behind internal
  adapters. Do not require a large finance framework.
- Use Catch2 3 for tests, pinned in the build configuration. Builds must support
  offline use when dependencies are already provisioned; dependency fetching
  must be explicit, not a surprise in a normal build.
- Use CI configurations for at least one supported compiler on Linux and
  macOS, with sanitizers in a dedicated configuration. Release benchmarks
  record compiler, flags, CPU, and dependency versions.
- Keep serialization format/version explicit (for example, versioned JSON
  experiment configurations); avoid making a particular data vendor or
  database a core dependency.

## 6. Market data and synthetic fixtures

Define `MarketDataSource`/import adapters that produce a normalized,
provenance-carrying `OptionQuoteSet`; adapters do not price, calibrate, or
silently repair quotes. Preserve raw quote fields and attach normalization,
filtering, timezone, underlying, rate/dividend, and corporate-action metadata.
Validate quote identity, timestamp, units, strike/maturity, bid/ask ordering,
and static arbitrage checks where enough market context exists.

Initially supply deterministic, hand-specified synthetic fixtures for unit and
integration tests. Mark every fixture as synthetic in file names, schema, and
reports. Synthetic surfaces may be generated from known parameters plus
explicitly specified noise; record generator, seed, and noise model. They are
not genuine historical data, do not establish calibration quality, and must
never be reported as market evidence. Real-data experiments require licensed
or otherwise authorized data, source attribution, timestamp and cleaning
policy, and a reproducible snapshot identifier; do not commit restricted data.

## 7. Milestones and acceptance tests

Each milestone is accepted only when tests and experiment metadata are
reproducible. Statistical tests use confidence thresholds and fixed seeds;
avoid brittle assertions on one stochastic sample.

1. **Foundation and build.** CMake library, core types, typed errors, RNG,
   configuration serialization.
   - Accept: clean Debug/Release builds; invalid domains fail explicitly;
     fixed seed and stream reproduce paths; separate streams differ; no global
     RNG state.
2. **Brownian simulation.** Standard and correlated paths.
   - Accept: increment means/variances and cross-correlation fall within
     predeclared statistical tolerances over sufficiently many seeded samples;
     endpoints and grid conventions are tested; invalid correlation matrices
     return errors.
3. **Black–Scholes and pricing.** Analytic price, payoff, Monte Carlo, CI,
   antithetic/control variates.
   - Accept: analytic prices satisfy put-call parity and bounds; Monte Carlo
     agrees with analytic values within a predeclared multiple of standard
     error across cases; variance-reduction estimators remain unbiased within
     test tolerance and report measured variance ratios.
4. **Implied volatility.** Robust inversion.
   - Accept: round trips over a grid of \(S,K,T,r,q,\sigma\); explicit
     boundary/out-of-bounds handling; residual and convergence status verified.
5. **Heston.** Simulation and analytic reference.
   - Accept: variance non-negativity under selected scheme; correlation and
     moment diagnostics; convergence toward analytic prices under grid/path
     refinement; Feller and non-Feller cases are both covered. The raw full-
     truncation state may be negative, but every variance passed to the stock
     evolution or square root must be non-negative.
6. **Volterra and rough Bergomi.** Baseline convolution then separately
   validated hybrid scheme.
   - Accept: covariance/variance of generated \(Y\) agrees with theoretical
     values within statistical tolerances; \(E[V_{t_i}]\) agrees with
     \(\xi_0(t_i)\) within sampling error and discretization limitations;
     correlation is validated; prices are compared under joint path/grid
     refinement to a declared reference. Include small and moderate \(H\),
     \(\rho\) near endpoints, and multiple grids.
7. **Calibration.** Synthetic quote sets, objective, bounds, diagnostics.
   - Accept: known synthetic parameters are recovered within declared
     tolerances for identifiable cases; ill-conditioned/non-identifiable
     cases are flagged rather than passed; multiple starts and residuals are
     included in reports. This does not count as market calibration.
8. **Greeks and hedging.** Pathwise / finite-difference Greek interfaces and
   discrete delta hedge.
   - Accept: BS Greeks match analytic values within estimated uncertainty;
     finite differences are stable over a reported bump study; hedge cashflows
     reconcile and P&L summaries include uncertainty and re-hedge schedule.
9. **Data, analytics, benchmarks.** Adapter schema, paired research reports,
   reproducible benchmark harness.
   - Accept: synthetic data provenance is visible end-to-end; import
     validation rejects malformed quotes; paired model comparisons use
     identical market cases and report uncertainty; performance reports pair
     runtime with accuracy/error and full environment metadata.

## 8. Numerical risks and mitigations

| Risk | Why it matters | Required mitigation / diagnostic |
|---|---|---|
| Rough-kernel singularity and low \(H\) | Naive left-point quadrature can materially bias \(Y\), variance and option prices. | Cell-average baseline; hybrid validation; covariance checks; grid refinement across \(H\). |
| Time-step bias mistaken for sampling error | More paths narrow the CI without fixing discretization bias. | Independent path-count and time-grid refinement; report both. |
| Correlation broken by independent approximations | Wrong leverage effect can change skew and hedge outcomes. | Construct correlated increments explicitly; verify empirical covariance and edge cases. |
| Heston variance negativity / bias | Naive Euler can create invalid square roots or biased prices. | Named scheme, explicit truncation policy, non-negativity tests, analytic benchmark. |
| Overflow/underflow in lognormal variance | Exponentiation can produce non-finite paths, especially for stress parameters. | Validate domains, detect non-finite intermediates, report context; no silent clipping. |
| Variance-reduction bias | Fitted control coefficients or dependent estimators can distort reported errors. | Document estimator, independent pilot or valid correction, benchmark against plain MC. |
| Calibration local minima / weak identification | A converged optimizer does not imply unique or stable parameters. | Multiple starts, bounds and residual diagnostics, sensitivity and profile checks. |
| Misleading implied-volatility inversion | Prices near intrinsic/upper bounds can be ill-conditioned. | Arbitrage bounds, bracketed solver, residual and conditioning diagnostics. |
| Parallel nondeterminism | Thread-dependent random consumption can invalidate replay/comparisons. | Stable path-indexed streams and deterministic reduction option; record execution mode. |
| Synthetic-to-market overclaim | A designed surface omits market microstructure and data issues. | Machine-readable provenance and explicit synthetic labels; no empirical claims absent real data. |
| Runtime-only benchmarking | A fast approximation can be inaccurately presented as superior. | Report accuracy/convergence alongside time and memory; include reference method. |

## 9. Research experiment plan

### Question A: pricing

1. Generate deterministic synthetic vanilla surfaces from declared parameter
   sets for Black–Scholes, Heston, and rough Bergomi. Include matched and
   intentionally misspecified generating models, documented strikes and
   maturities, and noise scenarios only when their magnitude/model is stated.
2. Fit each eligible model on the same in-sample quotes, using declared
   objective, weights, constraints, and starts. Keep synthetic truth and
   calibration output distinguishable.
3. Evaluate on held-out strikes/maturities and quote conditions. Report
   price/volatility errors by moneyness and maturity, calibration residuals,
   parameter stability, runtime, and uncertainty. Use paired resampling or
   another declared method for uncertainty in differences.
4. Repeat over paths and time grids to separate Monte Carlo noise and
   discretization bias from apparent model differences. Compare rough
   discretization schemes against the small-grid covariance/reference tests.

### Question B: hedging

1. Predeclare option set, initial calibration information, hedge frequency,
   re-hedging convention, rates/dividends, and whether transaction costs are
   excluded (initially yes).
2. Simulate common underlying scenarios from each declared data-generating
   model and hedge with each fitted model. Use paired scenarios/common random
   numbers where meaningful; do not calibrate on the hedge evaluation sample.
3. Report mean, standard deviation, quantiles and confidence intervals of
   terminal hedge error and cumulative P&L, stratified by maturity and
   moneyness. Include sensitivity to hedge frequency and numerical settings.
4. Distinguish “best against a synthetic generator” from real-world hedging
   evidence. Later genuine historical evaluation needs a timestamp-safe
   calibration/evaluation split, data provenance, survivorship/selection
   discussion, and realistic quote/hedging constraints.

### Reporting and decision rule

Pre-register comparison metrics and tolerances before inspecting outcomes.
Report all model specifications and failed/unstable runs, not only favorable
cases. A conclusion of improvement requires out-of-sample improvement with
uncertainty quantified, numerical error controlled below the observed effect,
and the data source clearly identified. Synthetic experiments can validate
methodology and reveal model behavior but cannot answer empirical market
performance on their own.

## 10. Explicit initial exclusions

Rough Heston, American options, jumps, stochastic interest rates, live market
connectors, vendor-specific formats, distributed simulation, GPU kernels,
transaction-cost optimization, and automatic model-selection claims are
extensions. Interfaces should not preclude them, but the initial engine must
not imply they are implemented or validated.
