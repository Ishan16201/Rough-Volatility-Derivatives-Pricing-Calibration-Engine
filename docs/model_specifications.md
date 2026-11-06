# Model specifications

## Black--Scholes

For constant annualized volatility \(\sigma\), risk-free rate \(r\), dividend
yield \(q\), and maturity \(T\), the risk-neutral stock is geometric Brownian
motion. Analytical call/put prices use discounted spot \(S_0e^{-qT}\) and
strike \(Ke^{-rT}\). The implementation provides delta, gamma, vega, theta
(calendar passage-of-time convention), rho, and bracketed implied volatility.
Terminal Monte Carlo samples the exact lognormal distribution; there is no
time-grid bias in that estimator.

## Heston

The risk-neutral dynamics and parameter units are in
[mathematical foundations](./mathematical_foundations.md). Parameters require
\(S_0>0,v_0,\theta\ge0,\kappa>0,\sigma_v\ge0,-1\le\rho\le1\); the Feller
condition is not required. Monte Carlo uses full-truncation Euler and retains
the raw next variance, while the stock uses the positive-part current
variance. Zero-boundary treatment is explicit but finite-step biased.
Characteristic-function pricing with adaptive quadrature is provided for
European options and has convergence diagnostics.

## Rough Bergomi

Rough Bergomi uses the RL Volterra driver with \(0<H<1/2,\eta>0\),
\(-1<\rho<1\), and strictly positive initial forward variance. It models
lognormal forward variance with a non-Markovian Gaussian driver. The selected
Volterra kernel scheme and its actual discrete variance enter the variance
correction. Asset paths use correlated Brownian increments and curve-based
discount/dividend inputs.

## Interpretation and limitations

Black--Scholes fixes variance; Heston is Markovian in \((S,v)\) with
mean-reverting CIR-type variance; Rough Bergomi has a history-dependent
Volterra driver and a specified initial forward-variance curve. These
structural differences do not by themselves imply better price fit or hedge
performance. Comparisons must hold market inputs and experiment design fixed
and report numerical uncertainty. Model comparison examples are synthetic;
no empirical superiority is established. Rough Heston is not implemented.
