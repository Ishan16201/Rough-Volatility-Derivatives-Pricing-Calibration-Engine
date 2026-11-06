# Calibration workflow and interpretation

The implementation accepts option observations with call/put, strike,
maturity, price, spot, risk-free/dividend discounts, and optional bid, ask,
implied volatility, reliability and instrument identifier. CSV ingestion is
available. A historical-data source interface is defined, but no real
historical market feed is bundled or used.

Three objective families are available: price residuals normalized by spot,
implied-volatility residuals normalized by 0.20, and bid/ask-aware price
residuals that are zero inside the spread. Reliability and capped inverse
spread weights can be selected; low-reliability and excessively wide quotes
can be excluded. These weighting rules mitigate quote-quality effects but
are not a complete market microstructure model.

Parameters have initial values, explicit bounds and optimize/fix switches.
Free dimensions are scaled to \([0,1]\). Deterministic shifted Halton starts
feed bounded Nelder--Mead. Reports include start outcomes, optimizer status,
residuals, local sensitivities, elapsed time, objective and numerical
diagnostics. Local optimization does not guarantee global optimality,
identifiability or statistical confidence.

Rough Bergomi's initial forward-variance curve is represented by 2--5
piecewise-linear knots. Optional regularization penalizes curvature in
log-forward variance. Calibration uses common random numbers with fixed
streams/settings; a separate evaluation stream assesses the selected
parameters. Model price Monte Carlo standard errors and Heston quadrature
diagnostics are distinct; the propagated objective sampling-error estimate
is only a first-order diagnostic.

Reproduce a synthetic Black--Scholes workflow:

```sh
cmake --preset release
cmake --build --preset release --target rvre_calibrate
./build/release/rvre_calibrate --model black_scholes --seed 20261106 \
  --starts 4 --output calibration_report.csv
```

See [calibration details](./calibration_details.md) for CSV columns, quote
filtering, parameter examples and further numerical limitations. Synthetic
recovery is a software check, not market evidence. No calibration has used
genuine option market data.
