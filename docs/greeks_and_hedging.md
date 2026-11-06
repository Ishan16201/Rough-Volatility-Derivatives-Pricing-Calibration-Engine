# Greeks and discrete-time hedging

## Sensitivities

Black--Scholes provides analytic delta, gamma, vega, theta and rho. For
non-analytic quantities RVRE exposes central finite differences, common
random-number coupling and repeated bump refinement. Parameters include spot,
volatility, rates, Hurst exponent, vol-of-vol and explicit initial
forward-variance curve knots. The reported difference between final bumps is
a diagnostic, not a bound. Monte Carlo sensitivities additionally require
path-count and grid convergence checks.

Pathwise derivatives are valid only under differentiability/integrability
conditions. Vanilla Black--Scholes pathwise delta/vega are admissible almost
surely; gamma at the payoff kink is not a plain pathwise derivative.
Likelihood-ratio estimators may handle payoff discontinuities but can have
high variance. RVRE does not silently substitute either method for the
documented finite differences.

## Hedge accounting

An experiment sells one European option at the market scenario's initial
option mark. **All hedge strategies for that scenario, seed and contract use
that same mark**; each strategy computes its own initial and rebalance deltas.
Initial cash is premium less stock purchase and transaction costs. Each
interval applies scenario financing growth and discrete dividend cash flow;
configured rebalance dates debit/credit traded shares and proportional
transaction cost. Terminal P&L is cash plus terminal stock value less payoff.
Turnover includes initial and subsequent trades.

The market mark is analytic for Black--Scholes, characteristic-function based
for Heston, and Monte Carlo for Rough Bergomi. The CSV records its numerical
error label; this uncertainty is not part of hedge-P&L intervals. No margin,
borrow constraints, market impact, asymmetric bid/ask, or funding spread is
modeled. Continuous-dividend cash flow uses the beginning-of-interval spot,
an explicit discrete convention.

The Black--Scholes strategy uses analytic delta; Heston uses a central spot
bump of characteristic-function prices; Rough Bergomi estimates its initial
implied volatility and uses a frozen Black--Scholes delta. The Rough Bergomi
proxy is not the true conditional Rough Bergomi Greek. Heston non-convergence
flags are retained. Results should not be treated as production hedges.

The cross-product experiment applies all strategies to the same simulated
paths within each synthetic market scenario and seed. Reported mean, sample
standard deviation, RMSE, mean absolute error, lower-tail loss, turnover,
costs and mean-error confidence interval distinguish sampling uncertainty
from model/grid/Greek uncertainty. The workflow uses declared parameters
ex ante; it is not a calibrated out-of-sample market evaluation. See
[risk and hedging details](./RISK_AND_HEDGING.md) for executable instructions.
