# Rough Bergomi simulation

This document specifies Rough Bergomi path simulation. Pricing and calibration
are separate modules; no empirical market-data conclusions are made.

## Model and parameter domains

On a risk-neutral measure, use
\[
Y_t=\sqrt{2H}\int_0^t(t-s)^{H-\frac12}\,dW_s,\qquad
V_t=\xi_0(t)\exp\left(\eta Y_t-\tfrac12\eta^2t^{2H}\right),
\]
\[
dS_t/S_t=(r_t-q_t)dt+\sqrt{V_t}\,dW_t^S,\qquad
d\langle W^S,W\rangle_t=\rho dt.
\]
Inputs satisfy \(0<H<1/2,\ \eta>0,\ -1<\rho<1,\ S_0>0\), and
\(\xi_0(t)>0\). Time and maturity are in years, \(V,\xi_0\) in
year\(^{-1}\), \(Y\) in year\(^{H}\), and \(\eta\) in year\(^{-H}\);
\(\rho,H\) are dimensionless. Curves are evaluated at every simulation
gridpoint and invalid, nonfinite or nonpositive values are errors.

`ForwardVarianceCurve` returns instantaneous initial forward variance.
`DiscountCurve` returns \(D(0,t)\), required to be positive and normalized to
one at zero (within \(10^{-12}\)). The provided piecewise-linear curve types
interpolate between knots and extrapolate flat beyond the first/last knots;
all supplied forward variances and discount factors must be positive. Users
may provide other curve implementations, which are validated at every
simulation point. The dividend curve is also represented as a normalized
discount-factor curve \(D_q(0,t)\). It need not be monotone: negative rates or
dividend yields are not silently prohibited.

## Volterra and variance correction

The simulator consumes the existing `VolterraPaths` implementations:
`cell_average_reference` or the explicitly selected
`hybrid_exact_diagonal`. The reference method's cell-average coefficients and
the hybrid method's midpoint coefficients, diagonal Brownian contribution,
and residual variance are described in
[STOCHASTIC_SIMULATION.md](./STOCHASTIC_SIMULATION.md).

For either scheme let \(c_{i,j}\) denote the coefficient of Brownian increment
\(\Delta W_j\) in its simulated \(Y_{t_i}\), and let \(b_i\) be the variance
of any independent diagonal residual (zero for the reference scheme). The
actual marginal variance of the discrete Gaussian driver is
\[
\nu_i=\operatorname{Var}(Y_{t_i}^{\rm scheme})
=\sum_{j<i}c_{i,j}^2\Delta t_j+b_i.
\]
At zero \(\nu_0=0\). Compute \(\nu_i\) from the selected scheme using
`volterra_discrete_variance`; do not substitute the continuous \(t_i^{2H}\)
except when defining/validating the continuous model. Set
\[
V_i=\xi_0(t_i)\exp\left(\eta Y_i^{\rm scheme}
-\tfrac12\eta^2\nu_i\right).
\]
Since \(Y_i^{\rm scheme}\) is centered Gaussian with variance \(\nu_i\), this
ensures \(\mathbb E[V_i]=\xi_0(t_i)\) for either approximate discretization.
It does not remove all discretization error in cross-time covariance or in
asset-option behavior. In the hybrid method, the marginal correction is
consistent with its exact diagonal plus midpoint older-cell approximation;
the method's approximate joint covariance limitation remains.

## Correlation and stock update

For every interval generate independent normals \(Z_i,Z_i^\perp\) and set
\[
\Delta W_i=\sqrt{\Delta t_i}Z_i,\qquad
\Delta W_i^S=\sqrt{\Delta t_i}
(\rho Z_i+\sqrt{(1-\rho)(1+\rho)}Z_i^\perp).
\]
The two increments have correlation \(\rho\). The Volterra driver uses
\(\Delta W_i\); the asset uses \(\Delta W_i^S\). The hybrid bridge residual
uses a separate domain-separated seeded stream.

The discrete update is
\[
\log(S_{i+1}/S_i)=
\log\!\frac{D_r(0,t_i)}{D_r(0,t_{i+1})}
-\log\!\frac{D_q(0,t_i)}{D_q(0,t_{i+1})}
-\tfrac12V_i\Delta t_i+\sqrt{V_i}\Delta W_i^S.
\]
The curve contribution is the integrated interval rate difference, not a
pointwise rate approximation. Under these assumptions \(V_i\) is measurable
before the interval-\(i\) Brownian increment; the lognormal step therefore
satisfies
\[
\mathbb E[S_{i+1}\mid\mathcal F_{t_i}]
=S_i\frac{D_r(0,t_i)}{D_r(0,t_{i+1})}
  \frac{D_q(0,t_{i+1})}{D_q(0,t_i)}.
\]
Thus \(\mathbb E[S_TD_r(0,T)/D_q(0,T)]
=S_0D_q(0,0)/D_r(0,0)=S_0\) for normalized curves. Simulation diagnostics
report the empirical driver correlation, the mean and standard error of
\(V_i\) against \(\xi_0(t_i)\), and discounted-forward mean/error. These
diagnostics are estimates from the requested finite path sample, not guarantees
that a finite realization exactly equals its target.

Every nonfinite exponent, variance, asset, or discounted-forward observation
fails explicitly with path/time context where available. No clipping is used.
The output contains path-major asset, variance, and Volterra-driver values.

## Example

Build and run:

```sh
cmake --preset debug
cmake --build --preset debug --target rvre_rough_bergomi_example
./build/debug/rvre_rough_bergomi_example paths.csv 128 32
```

The example writes CSV data to the requested file, diagnostics to stderr, and
uses a fixed seed and synthetic constant curves. It is a reproducible software
demonstration, not an empirical calibration or market-data result.

## Validation

Tests check strict parameter/curve domains, path positivity, deterministic
replay, empirical driver correlation, discrete forward-variance consistency,
discounted-forward martingale consistency within Monte Carlo standard errors,
and convergence of the discrete variance correction to \(t^{2H}\) for both
Volterra schemes over representative Hurst values. Time-step and finite-sample
errors are distinct: the exponential correction targets the actual discrete
Gaussian variance, while the grid-refinement tests assess its difference from
the continuous limit.
