# Conventional benchmark models

The conventional-model module provides Black--Scholes and Heston benchmarks.
Calibration and hedging interfaces are described in their respective docs.
All times are year fractions; rates, variance, and volatility use annualized
units.

## Black--Scholes

With constant volatility \(\sigma\), constant continuously compounded rates
\(r,q\), and maturity \(T\), discounted spot and strike are
\(S_0e^{-qT}\) and \(Ke^{-rT}\). Call/put prices use the standard \(d_1,d_2\)
formula and satisfy put-call parity. Analytical Greeks are:

* Delta is with respect to spot.
* Gamma is the second derivative with respect to spot, common to calls and
  puts.
* Vega is the derivative with respect to volatility expressed per 1.0 change
  in annualized volatility, not per volatility point.
* Theta is calendar-time decay per year: the derivative with respect to
  passage of time, i.e. minus the derivative with respect to remaining
  maturity, holding spot, rates, dividend yield, and volatility fixed.

Greeks require positive maturity and volatility. The price function supports
zero maturity and zero total variance using discounted intrinsic value. The
implied-volatility solver checks no-arbitrage bounds, expands a volatility
bracket, then uses a safeguarded Brent-style method with explicit price and
volatility tolerances and iteration limits. A finite implied volatility is not
returned at the upper no-arbitrage bound.

Black--Scholes terminal Monte Carlo samples the exact lognormal terminal
distribution. It has sampling error but no time-grid discretization bias, and
is used as an analytic implementation benchmark.

## Heston model

The risk-neutral dynamics are
\[
\frac{dS_t}{S_t}=(r-q)\,dt+\sqrt{v_t}\,dW_t^S,\qquad
dv_t=\kappa(\theta-v_t)\,dt+\sigma_v\sqrt{v_t}\,dW_t^v,\qquad
d\langle W^S,W^v\rangle_t=\rho\,dt.
\]
Inputs require \(S_0>0,\ v_0,\theta\ge0,\ \kappa>0,\ \sigma_v\ge0\),
\(-1\le\rho\le1\), and finite constant rates. \(v_0,\theta\) have units
year\(^{-1}\), \(\kappa\) year\(^{-1}\), \(\sigma_v\) year\(^{-1/2}\), and
\(\rho\) is dimensionless. The Feller condition is not required: it determines
whether zero is attainable for the continuous CIR process, not whether the
model is mathematically defined.

### Full-truncation Euler Monte Carlo

The implementation stores a raw Euler variance \(v_n\), then uses
\(v_n^+=\max(v_n,0)\) in the drift, diffusion, and asset update:
\[
v_{n+1}=v_n+\kappa(\theta-v_n^+)\Delta t+
 \sigma_v\sqrt{v_n^+}\sqrt{\Delta t}\,Z_n^v,
\]
\[
\log(S_{n+1}/S_n)=(r-q-\tfrac12v_n^+)\Delta t+
\sqrt{v_n^+}\sqrt{\Delta t}
(\rho Z_n^v+\sqrt{1-\rho^2}Z_n^\perp).
\]
The next raw variance is not projected or clipped. If it is negative, the
next step uses zero instantaneous variance but retains the raw state, allowing
the mean-reversion drift to move it back above zero. This explicit boundary
rule avoids a square root of a negative number, but introduces finite-step
bias and can yield negative raw states. Diagnostics include the fraction of
updates where \(v_n\le0\), the observed raw variance range, and the
Monte Carlo discounted-forward mean/error. The scheme is not claimed to be
high-order or unbiased for coarse steps.

### Characteristic-function pricing

European prices use \(C=S_0e^{-qT}P_1-Ke^{-rT}P_2\), with the Heston
characteristic function evaluated in the stable \(g=(b-d)/(b+d)\) form and
adaptive Simpson quadrature of the two probability integrals. Put values follow
from parity. The continuous model's characteristic function is distinct from
the full-truncation Euler transition law, so exact numerical equality with
finite-step Monte Carlo is not expected.

Diagnostics report adaptive quadrature error estimates, integrand evaluation
counts, the final frequency cutoff, and convergence status. To expose
frequency-tail sensitivity, each integral is evaluated at the configured
cutoff and twice that cutoff; their price difference is included in the error
estimate and convergence flag. This is a practical sensitivity check, not a
rigorous tail bound. A non-converged flag must be treated as an explicit
numerical warning even though a price estimate is returned. Zero vol-of-
variance is treated analytically as deterministic variance and reduces to
Black--Scholes with integrated variance
\[
\int_0^T v_t\,dt=\theta T+(v_0-\theta)(1-e^{-\kappa T})/\kappa.
\]

## Comparing the three models

Black--Scholes fixes variance; Heston variance is a continuous mean-reverting
CIR diffusion with square-root innovations; rough Bergomi specifies a
lognormal forward variance driven by a non-Markovian Volterra Gaussian
process. Heston's state is Markovian in \((S,v)\); rough Bergomi retains
history-dependent covariance through its Volterra driver. Black--Scholes is a
useful analytic control but has no stochastic-volatility smile dynamics.
Heston can produce skew and smile, but its variance paths are not rough in the
fractional-kernel sense. Rough Bergomi directly models rough log-forward
variance but depends on a chosen Volterra discretization and is more
computationally demanding.

`rvre_model_comparison` uses the same spot, strike, maturity, rates, dividend
yield, initial variance, path count, and seed for each model, and the same
time-step count for Heston and rough-Bergomi simulation. Its illustrative
Heston values \((\kappa,\theta,\sigma_v,\rho)\) and rough-Bergomi
\((H,\eta,\rho)\) are intentionally different model-specific dynamics.
Black--Scholes volatility is set to the square root of the common initial
variance. The analytic/characteristic-function prices are model benchmarks;
MC intervals describe sampling error only, not time bias or quadrature-tail
error. This synthetic comparison is not empirical evidence or calibration.

Build and run the comparative fixture:

```sh
cmake --preset release -DRVRE_BUILD_BENCHMARKS=ON
cmake --build --preset release --target rvre_model_comparison
./build/release/rvre_model_comparison
```

Tests use known Black--Scholes prices and Greeks, an independent canonical
Heston characteristic-function benchmark, exact deterministic-variance
Heston/Black--Scholes reduction, and MC comparisons. Path-count changes
measure sampling error; Heston time-step changes expose the separate
full-truncation discretization effect.
