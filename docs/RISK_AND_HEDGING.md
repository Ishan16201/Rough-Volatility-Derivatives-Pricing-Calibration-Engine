# Risk sensitivities and discrete delta hedging

## Sensitivities

`risk::estimate_sensitivities` accepts a model-specific option-pricing
callback. Its state carries spot, strike, year-fraction maturity, Black-style
volatility, continuously compounded rates, roughness \(H\), vol-of-vol \(\eta\),
and an initial forward-variance curve. The factors are:

* Delta: \(\partial V/\partial S_0\).
* Gamma: \(\partial^2 V/\partial S_0^2\).
* Vega: \(\partial V/\partial \sigma\), per 1.0 change in annualized
  volatility (not per volatility point).
* Theta: calendar passage-of-time decay, \(-\partial V/\partial T\), with
  other inputs held fixed.
* Rho: \(\partial V/\partial r\), per 1.0 annual-rate change.
* \(H,\eta,\xi_0(t_i)\): derivatives with respect to the selected model
  parameter or forward-variance knot.

The generic estimator uses central first differences and a central second
difference for gamma. It halves the configured initial bump repeatedly and
reports the finest estimate, final bump, relative disagreement between the
last two refinements, and a convergence flag. It never silently claims
convergence: users should tighten/loosen tolerances only after checking a
step-size sequence and stability for their pricing callback. Perturbations
share the same seed and domain-separated stream for both sides, which is the
common-random-number (CRN) coupling needed to reduce Monte Carlo noise in
differences. Gamma also reuses that stream for its center value. Do not compare
the resulting numerical error estimate with Monte Carlo confidence intervals
as if they were interchangeable.

The existing analytic Black--Scholes delta, gamma, vega and theta are the
preferred estimators for that model; `black_scholes_rho` supplies analytic
rho. For Heston, finite differences of characteristic-function prices avoid
pathwise noise and price the specified continuous-time model, subject to
quadrature diagnostics. For Rough Bergomi and other callback-priced models,
finite differences are model-agnostic, but need CRN and grid/path refinement:
they inherit Monte Carlo and time-discretization bias. The caller must map
state shocks to the target model correctly; e.g., a model's volatility-of-
volatility parameter is not interchangeable with Black--Scholes vega.

Pathwise differentiation is valid when the payoff and path map are
differentiable almost surely and the derivative is integrable. In lognormal
Black--Scholes, pathwise delta and vega meet this condition for vanilla
options because the strike kink has probability zero; the library's analytic
Greeks are lower-noise and already available. Pathwise gamma for a vanilla
payoff requires differentiating the exercise indicator and is not an ordinary
pathwise derivative. Likelihood-ratio/score estimators can handle such payoff
kinks and can estimate gamma, but tail weights can have high variance; none is
silently substituted for the reported estimators. For rough volatility,
likelihood-ratio estimators would need a carefully derived score for the
discretized joint Volterra/asset law. This implementation uses transparent
finite differences instead.

## Hedge model estimators and limitations

The hedge module uses:

* Black--Scholes analytic delta.
* Heston characteristic-function prices with a central spot difference.
* Rough Bergomi Monte Carlo price at inception, inverted to a Black--Scholes
  implied volatility, then a **frozen implied-volatility Black--Scholes delta**
  through time. This is an explicitly approximate hedge, not a conditional
  Rough Bergomi value/Greek. The model's rough path history is not conditioned
  on at later hedge times, so the proxy should not be read as a true
  rough-Bergomi delta.
  The CSV separately flags whether the Heston characteristic-function
  calculations used for the initial delta and each subsequent hedge delta
  reported convergence. Non-converged estimates are retained but are not
  silently treated as validated; the CLI emits a warning and users should
  inspect or rerun those configurations with refined quadrature.

The generic sensitivity API offers multi-step finite-difference convergence
diagnostics for user-supplied Heston/Rough-Bergomi pricing callbacks. The
built-in hedging Heston spot bump and Rough Bergomi frozen-IV proxy are chosen
for an executable paired study; conduct bump and pricing-grid convergence
checks before treating them as production hedge Greeks.
The output stores the scenario market mark's numerical-error diagnostic and
labels whether it is analytic, a Heston characteristic-function error
estimate (including a non-convergence label), or a Rough Bergomi Monte Carlo
standard error. That quote uncertainty is not folded into the hedge-error
confidence interval.

## Accounting convention

The experiment sells one European option at the market scenario's initial
option mark and buys the hedge model's delta in shares. Every hedge strategy
within that scenario/seed uses the same mark. Initial cash is option premium
minus stock purchase and initial proportional transaction costs. At each
market-grid interval, cash earns the scenario risk-free financing growth and
the stock position receives a discrete dividend cash flow. At configured
rebalance times, cash pays for the share change and proportional costs.
At maturity, P&L is cash plus stock value minus the contractual payoff. No
borrow constraints, margin, discrete dividends, bid/ask asymmetry, market
impact, or funding spread are modeled. Continuous dividend yields/discount
curves are represented as interval cash flows using the beginning-of-interval
spot; this convention is a discrete approximation and is shared by all hedge
models within a scenario.

Turnover includes the initial delta purchase and subsequent absolute share
trades, each marked at that trade's spot. Transaction costs are the configured
proportional rate times marked turnover. Financing and dividend flows are
applied each simulation interval, even when the hedge is not rebalanced.
Rebalancing stops before the terminal payoff; there is no spurious terminal
trade.

## Paired out-of-sample experiment

`rvre_hedge` builds synthetic Black--Scholes, Heston, and Rough Bergomi
risk-neutral market scenarios and applies all three hedge specifications to
each same path set. This full cross-product makes market-model misspecification
visible while holding option terms, seeds, time grids, rebalance schedules,
and costs fixed. these are deterministic synthetic model fixtures, not calibration
observations or empirical outcomes. There is no calibration step: parameters
are declared ex ante, and the Rough Bergomi market mark and frozen-IV proxy
use a domain-separated pricing stream from the hedge evaluation paths. No
same-observation calibration/evaluation result is reported.

The experiment runs multiple explicit seeds. For each market/hedge/frequency
cell it reports pooled path-level mean hedge error, sample standard deviation,
RMSE, mean absolute error, lower 5th-percentile P&L, worst-5%-tail expected
shortfall loss (positive values are losses), mean turnover, mean transaction
cost, standard error of mean hedge error, and a normal-approximation
confidence interval. Tail summaries use empirical order statistics. Paths
share a seed and market scenario across hedge strategies, so comparisons are
paired; strategy differences should be analyzed pathwise rather than by
treating model/hedge cells as independent samples. Confidence intervals are
sampling uncertainty only and do not include model, discretization, parameter,
or Greek-estimator uncertainty.

Build and run:

```sh
cmake --preset release
cmake --build --preset release --target rvre_hedge
./build/release/rvre_hedge --paths 8 --steps 8 --seed 731927 \
  --seed-count 1 --cost 0.001 --output hedging_experiment.csv
```

The CSV contains summary rows, one terminal accounting row per path, and
separate rows for every interval's financing/dividend flow, labeled as either
`cash_flow` or `rebalance`. Increase paths and seeds for
stable tail estimates, increase the grid to test simulation bias, and compare
multiple rebalance frequencies. The Rough Bergomi market uses the reference
Volterra discretization; Heston uses full-truncation Euler. Different market
models are distinct risk-neutral scenarios, not a claim that identical
Brownian realizations define the same economic outcome across models.
Each row also records maturity, strike, option right, cost rate, simulation
grid, seed count, and confidence level.
