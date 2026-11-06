# Monte Carlo European option pricing

The pricing module consumes Rough Bergomi paths and prices European calls and
puts. It does not calibrate parameters or read market data.

For strike \(K\), maturity \(T\), and terminal asset \(S_T\), the reported
payoff sample is discounted by the risk-free discount curve:
\[
\widehat C=D_r(0,T)\frac1N\sum_{i=1}^N(S_T^{(i)}-K)^+,\qquad
\widehat P=D_r(0,T)\frac1N\sum_{i=1}^N(K-S_T^{(i)})^+.
\]
One path simulation per maturity and batch is reused across requested strikes
and both payoff rights. The requested time-step count defines a uniform
year-fraction grid for each maturity. The stock update and variance
discretization are those documented in [ROUGH_BERGOMI.md](./ROUGH_BERGOMI.md).
The result is therefore subject to time-discretization bias in addition to
Monte Carlo sampling error; the confidence interval describes only the latter.

## Estimators and uncertainty

The standard error is the sample standard deviation of independent estimator
observations divided by the square root of their count. A two-sided normal
confidence interval uses the configured confidence level and a normal-quantile
approximation; for small samples or very skewed payoffs, coverage can differ
from nominal.

* **Standard:** each path payoff is one observation.
* **Antithetic:** the simulator generates Brownian paths in sign-reversed
  pairs. Each pair's average discounted payoff is one observation, so the
  standard error uses \(N/2\) observations. Both total paths and batch size
  must be even. No claim is made that antithetics always reduce error.
* **Stock control variate:** for each option, a pilot sample estimates
  \(\beta=\operatorname{Cov}(X,Y)/\operatorname{Var}(Y)\), where \(X\) is the
  discounted payoff and \(Y=D_r(0,T)S_T\). Evaluation uses
  \(X+\widehat\beta(S_0D_q(0,T)-Y)\). The pilot has domain-separated seeded
  streams from evaluation and is never reused to estimate the price. Each
  option receives its own pilot coefficient. The reported evaluation sample
  standard error does not include uncertainty in the estimated coefficient;
  conditional on the independent pilot, it is the standard error of the
  frozen-coefficient evaluation estimator.

Random streams are reproducible for a fixed seed, batch layout, maturity
ordering, and model configuration. Each maturity/batch and pilot/evaluation
purpose receives a deterministic domain-separated stream. Changing batch size
changes stream partitioning and is not expected to preserve the same path
sequence. The batch size is a target per simulation call. Non-antithetic paths are
distributed as evenly as possible across the required number of batches so
each simulation has at least two paths for path diagnostics. For example, a
three-path run with target batch size two uses a single three-path batch.
Results report the seed, path count, target batch size, pilot count, time
steps, method, independent observation count, coefficient, estimate, and
confidence interval.

The control expectation \(S_0D_q(0,T)\) follows from the model's normalized
discount curves and risk-neutral drift. It is not inferred from the simulated
evaluation sample. The estimator reports finite results explicitly and throws
on invalid inputs, nonfinite payoffs, or a degenerate pilot control variance.

## Executable

Build `rvre_price` and run, for example:

```sh
cmake --preset release
cmake --build --preset release --target rvre_price
./build/release/rvre_price --paths 500 --batch-size 100 --steps 8 \
  --seed 20261106 --strikes 100 --maturities 0.5 \
  --methods standard,antithetic,control
```

The executable writes metadata as a comment line and estimate records as CSV.
It uses deterministic synthetic flat forward variance and discount curves.
These fixtures are examples only and are not market observations or empirical
calibration results. `--help` lists the model and run options.

## Numerical validation

Tests cover analytical Black--Scholes parity and deterministic boundaries,
constant-variance convergence of the Monte Carlo estimate to Black--Scholes,
replay, antithetic independent-observation counts, control-pilot isolation,
and invalid inputs. Use independent grid-refinement experiments to investigate
time bias; increasing path count narrows sampling error but does not remove
time-discretization error. The optional pricing benchmark reports estimates,
standard errors, intervals, and elapsed time for each variance-reduction
method. Its timings are environment-specific and no universal speedup or
empirical pricing improvement is implied.
