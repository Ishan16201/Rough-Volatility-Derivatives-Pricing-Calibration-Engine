# RVRE research and release report

## Motivation and background

Black--Scholes is the analytic constant-volatility benchmark; Heston adds
Markovian mean-reverting stochastic variance; rough-volatility research
studies models whose volatility exhibits rough temporal behavior. Foundational
references include Black and Scholes, “The Pricing of Options and Corporate
Liabilities” (1973); Heston, “A Closed-Form Solution for Options with
Stochastic Volatility with Applications to Bond and Currency Options” (1993);
Gatheral, Jaisson and Rosenbaum, “Volatility is rough” (2018); and Bayer,
Friz and Gatheral, “Pricing under rough volatility” (2016). These works
motivate model classes; this repository does not reproduce or claim their
empirical results.

The research question is whether rough volatility improves option fit and
hedging economically over conventional stochastic volatility. The current
repository has the software to run synthetic comparisons, but not the
verified market data or validated Rough Bergomi conditional Greeks needed to
answer that empirical question.

## Mathematical formulation and assumptions

Black--Scholes uses constant volatility. Heston uses a CIR-type variance
diffusion with correlated stock/variance Brownian motions. Rough Bergomi uses
\[
Y_t=\sqrt{2H}\int_0^t(t-s)^{H-\frac12}dW_s,\quad
V_t=\xi_0(t)e^{\eta Y_t-\frac12\eta^2t^{2H}},\quad
dS_t/S_t=(r_t-q_t)dt+\sqrt{V_t}dW_t^S,
\]
with \(d\langle W^S,W\rangle_t=\rho dt\). The simulation uses the discrete
scheme's marginal driver variance in the exponential correction. This is the
Riemann--Liouville driver; it is not standard fractional Brownian motion and
does not have the standard fBm off-diagonal covariance. All model dynamics
are risk-neutral specifications with year-fraction time, positive prices,
explicit discount/dividend curves and explicit parameter constraints.
Detailed formulas and numerical assumptions are in
[mathematical foundations](./mathematical_foundations.md) and
[model specifications](./model_specifications.md).

## Numerical methodology

RVRE implements Black--Scholes analytical pricing, Heston characteristic
function and full-truncation Euler Monte Carlo, and Rough Bergomi simulation
with reference cell-average and optional hybrid Volterra discretization.
Pricing uses discounted European payoffs, optional antithetic variates and a
pilot-fitted stock control variate. Calibration uses bounded multi-start
optimization and synthetic surfaces. Greeks use analytical formulas where
available and common-random-number finite differences otherwise. Hedging is
discrete, with explicit financing, dividends and transaction costs. Detailed
equations and caveats are in the linked module documents.

## Calibration methodology and experiment design

Calibration supports price, implied-volatility and bid/ask-aware objectives
with declared quote weights and exclusion rules. Parameters are bounded and
scaled; deterministic multi-start Nelder--Mead is used. Rough Bergomi uses
common random numbers across candidate parameters and a compact 2--5 knot
forward-variance curve. A separate seed evaluates the selected calibration.
Synthetic parameter-recovery tests constrain one parameter while holding
others fixed; this is not evidence that a full market calibration is
identified. See [calibration](./calibration.md).

Model comparisons use common spot, strike, maturity, rates, dividend, initial
variance, path count and seed where those inputs are shared. Model-specific
parameters remain different by construction. Hedging strategies share
scenario paths and the scenario's initial option mark, but use their own
model deltas and hedge schedules. This creates a controlled synthetic
comparison, not a held-out market evaluation. Detailed experimental
conventions and limits appear in [Greeks and hedging](./greeks_and_hedging.md)
and [benchmarking](./benchmarking.md).

## Results actually obtained

The 53-case unit suite passed in Debug, Release and AddressSanitizer/
UndefinedBehaviorSanitizer builds. The focused hedge test verifies that all
strategies share a scenario-generated Black--Scholes market mark while
retaining different model deltas. Exact command outcomes and toolchain are in
[validation](./validation.md).

One Release cross-model smoke (spot/strike 100, maturity 1 year,
\(r=3\%\), \(q=1\%\), initial variance 0.04, 6,000 paths, 32 stochastic-model
steps, seed 261106) returned:

| Model/pricer | Price | Reported uncertainty |
|---|---:|---:|
| Black--Scholes analytic | 8.82732 | exact formula |
| Black--Scholes terminal MC | 8.94630 | SE 0.17749; 95% CI [8.59842, 9.29418] |
| Heston characteristic function | 8.25285 | quadrature diagnostic 1.34e-7 |
| Heston full-truncation MC | 8.29947 | SE 0.12923; 95% CI [8.04617, 8.55276] |
| Rough Bergomi reference MC | 7.87194 | SE 0.12822; 95% CI [7.62064, 8.12324] |

This is one synthetic parameterization and finite-sample run; Heston and
Rough Bergomi have their own dynamics despite common spot, rates, dividend,
initial variance and correlation. The intervals exclude time-discretization
and model error. These numbers support code-path smoke checking only, not a
model ranking.

The paired hedge smoke used 8 paths, 8 steps and one seed. It emitted 36
Heston characteristic-function delta non-convergence flags. We preserve those
flags and treat Heston hedging convergence as unresolved rather than claiming
the flagged values are reliable. The low path count also makes its P&L/tail
statistics unsuitable for inference.

No experiment in this repository uses verified historical option quotes.
Consequently, there are no empirical calibration, out-of-sample market
pricing, realized hedge performance, or economic-improvement results to
report. Synthetic parameter recovery and model cross-products are software
validation only.

## Limitations and future work

The requested Riemann--Liouville driver is not standard fractional Brownian
motion and does not have the standard fBm off-diagonal covariance. The hybrid
scheme approximates joint bridge covariance. Heston Euler and Rough Bergomi
time discretization bias remain. Heston quadrature error estimates are not
rigorous tail bounds. Rough Bergomi hedging uses a frozen-IV delta proxy.
Calibration confidence/identifiability diagnostics are incomplete for
research claims. No genuine market-data source, quote-cleaning policy,
calendar alignment, or empirical validation is included. Parallel
performance, Greeks convergence across parameter ranges, and detailed
module-level profiling remain future work.

## Evidence classification

* **Software correctness:** unit/integration behavior covered by the reported
  builds and tests.
* **Mathematical validation:** equation and discretization checks documented
  in [mathematical foundations](./mathematical_foundations.md), including the
  RL/fBm inconsistency.
* **Synthetic validation:** reproducible simulation, pricing, calibration and
  hedge experiments only.
* **Empirical validation:** incomplete; no genuine option data used.
* **Unresolved:** discretization bias, numerical diagnostics noted above,
  approximation quality over broad parameter domains, and economic value.
