# Calibration framework (detailed reference)

The calibration module fits Black--Scholes volatility, the five standard
Heston parameters, or a Rough Bergomi parameterization to a validated set of
European option observations. It supports synthetic fixtures and user-supplied
CSV; it does not fetch historical market data or imply that any real-market
calibration has been performed.

## Observations and data provenance

`OptionObservation` stores call/put right, strike, year-fraction maturity,
market price, optional bid/ask and implied volatility, spot, risk-free and
dividend discount factors, reliability in `[0,1]`, and an optional instrument
identifier. Inputs are checked for finite positive spot/strike, nonnegative
maturity/prices, valid discounts and quote ordering. CSV ingestion requires
`right,strike,maturity,market_price,underlying_price,risk_free_discount,dividend_discount`;
`instrument_id,bid,ask,implied_volatility,reliability` are optional columns.
Discount fields are discount factors, not rates. An integration may implement
`HistoricalOptionDataSource`; no historical vendor or feed is bundled.

All observations in a Rough Bergomi fit must share a spot, and discounts at
the same maturity must agree. The fitted forward-variance curve is a
piecewise-linear curve in variance, represented by 2--5 user-selected knots.
The knots bound the fitted curve and should be selected with the observation
maturities in mind; adding knots does not create information in the data.
Optional regularization penalizes squared changes in adjacent slopes of
log-forward variance. It is zero by default, is scaled in the optimizer's
dimensionless objective, and is a modeling choice rather than a statistical
prior.

## Objectives and quote weighting

The objective is a mean weighted squared residual:

* Price residuals are divided by spot.
* Implied-volatility residuals are divided by 0.20 (20 volatility points).
  Provided market implied volatilities take precedence; otherwise they are
  inverted from market prices.
* Bid/ask-aware price residuals are zero inside the quoted interval and
  otherwise measure distance to the nearest side, divided by the greater of
  half-spread and the configured minimum relative width times spot.

Weights can be equal, reliability-based, or capped inverse relative bid/ask
width multiplied by reliability. Quotes below the reliability threshold or
above the maximum relative spread are excluded. Weights are normalized to
mean one over included observations, so a large quote count does not itself
scale the reported objective. These heuristics reduce the influence of weak
quotes but are not a substitute for cleaning, timestamp alignment, or a
market-specific liquidity model. Bid/ask mode requires at least one quote;
observations without bid/ask are excluded from that objective. Implied
volatilities are calculated only for an implied-volatility objective; other
objective reports leave the IV columns blank rather than treating a noisy
Monte Carlo price outside Black--Scholes bounds as a valid volatility.

## Optimization and numerical error

Parameters have explicit initial values, lower and upper bounds, and an
`optimize` flag. Free parameters are scaled to `[0,1]`. A deterministic
multi-start design uses shifted Halton points initialized by the optimizer
seed; each start uses bounded Nelder--Mead. The result includes every start,
the best objective, termination status, finite-difference local sensitivities,
and per-observation residuals. Convergence is numerical termination, not proof
of global optimality or parameter identifiability. Inspect alternative starts,
bound hits, sensitivities, and parameter stability before interpreting a fit.

Black--Scholes prices are analytic. Heston uses the characteristic-function
pricer; each observation's discount factors are represented by equivalent
continuously compounded rates for that maturity. The configured Heston
quadrature error estimate is checked against a price-error budget and is
reported separately from sampling error. This adaptive quadrature estimate is
a numerical diagnostic, not a rigorous tail bound. Each residual also carries
the characteristic-function convergence flag; a false flag is preserved in
the CSV and reported by the CLI even when its estimated error remains within
the configured budget.

Rough Bergomi uses Monte Carlo with a fixed seed, stream namespace, batch
layout, path count, time-step count, and Volterra scheme across candidate
parameters (common random numbers). This reduces noise in objective
differences but does not remove sampling error, time-discretization bias, or
Volterra approximation error. A second deterministic seed/stream pair,
recorded in the report, evaluates the selected parameters independently.
Reports distinguish each model price's Monte Carlo standard error from its
deterministic numerical-error estimate. The estimated objective sampling
error is a first-order diagnostic propagated from price standard errors; it
is not a confidence interval for parameters or a guarantee that a small
objective improvement is significant. Increase paths and refine the time
grid to check stability, keeping those two error sources distinct.

## Reproducible workflow

Generate and fit a synthetic surface:

```sh
cmake --preset debug
cmake --build --preset debug --target rvre_calibrate
./build/debug/rvre_calibrate --model black_scholes --seed 20261106 \
  --starts 4 --output calibration_report.csv
```

For Rough Bergomi, additionally choose `--paths` and `--steps`; the synthetic
surface and fit use the same discretization settings and seed. Supply a CSV
with the required columns using `--input observations.csv`. Choose
`--objective price`, `--objective implied_volatility`, or
`--objective bid_ask`. The CSV report includes parameters, start diagnostics,
residuals, local sensitivities, random stream metadata, numerical diagnostics,
and an independent validation objective. Synthetic parameter recovery is a
software test only; it is not evidence of empirical fit or economic value.

The current optimizer is intentionally lightweight and derivative-free.
Nelder--Mead may converge slowly or to local minima in high-dimensional or
weakly identified fits. No automatic confidence intervals, profile
likelihoods, or global-optimality claims are made. Calibrate progressively,
fixing or bounding weakly identified parameters, and assess recovery under
multiple synthetic parameter sets before research conclusions.
