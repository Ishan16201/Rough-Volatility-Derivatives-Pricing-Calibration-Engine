# Numerical methods and error controls

## Randomness and simulation

RVRE uses explicit Philox4x32-10 random streams and a fixed seed/stream
identity. Brownian increments are scaled normal draws. The Box--Muller
transform relies on floating-point transcendental functions, so the exact
replay contract is for a fixed build/runtime, not all platforms.

Brownian and Volterra output use contiguous path-major storage. Direct
cell-average Volterra simulation is the transparent reference implementation.
Its work is \(O(PN^2)\), kernel storage \(O(N)\), and retained path storage
\(O(PN)\). The hybrid diagonal candidate is also \(O(PN^2)\) in this version;
it has no presumed complexity advantage. Tests/benchmark report numerical
variance and covariance errors separately from elapsed time.

Rough Bergomi uses a left-point variance in the log-stock step and a
scheme-specific marginal variance correction. Heston uses full-truncation
Euler. Refining the time grid probes these approximation errors; increasing
paths does not remove them.

## Sampling and variance reduction

Plain Monte Carlo uses sample variance over paths. Antithetic simulation
reverses Gaussian shocks in pairs and treats a pair-average payoff as one
independent observation. The stock control variate estimates its coefficient
on a separately seeded pilot sample and freezes it for evaluation. Each
method remains separately selectable. A reported confidence interval is a
normal approximation based on independent estimator observations; small
samples/skewed payoffs may have different coverage. Pilot coefficient
uncertainty is not included in the conditional evaluation standard error.

Sampling error, time-step error, Volterra approximation error, and model
specification error are not interchangeable. Price convergence should vary
path count at fixed grid and vary grid at adequate path count, reporting both.

## Conventional pricing and inversion

Black--Scholes prices and Greeks are analytical. Implied volatility uses a
no-arbitrage bound check, volatility bracket expansion and safeguarded
Brent-style root solving. Prices at the upper bound do not yield a finite IV;
expiration/zero-variance cases are handled explicitly.

Heston characteristic-function evaluation uses adaptive quadrature and a
doubled-frequency-cutoff comparison. The reported difference is a numerical
diagnostic, not a rigorous tail bound. Non-convergence flags must be retained.
Heston characteristic-function prices and full-truncation MC prices represent
different numerical objects: the former targets the continuous model, the
latter has time discretization and sampling error.

## Calibration and sensitivities

Calibration uses scaled bounded parameters and deterministic multi-start
Nelder--Mead. Rough Bergomi evaluations hold common random numbers fixed
across candidates. The independent validation evaluation is distinct from
the fit objective. Model price Monte Carlo standard errors and deterministic
numerical diagnostics are reported separately. An objective-level propagated
sampling diagnostic is approximate; it is not parameter confidence.

Finite-difference sensitivity estimates use common random numbers and a
sequence of reduced bumps. Their refinement disagreement is not a complete
error bound; users must check bump, path, grid and solver convergence.
Black--Scholes analytic Greeks and Heston characteristic-function finite
differences are preferred where available. Rough Bergomi Greeks remain
finite-difference estimates of a discretized stochastic price.

## Failure and stability

Inputs and generated values are checked for invalid domains and nonfinite
values. Kernel clipping is not used. Heston CF non-convergence and calibration
optimizer termination status remain observable. Extreme exponential states
may exceed floating-point range; those cases fail explicitly instead of
silently clipping.
