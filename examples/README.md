# Examples

The legacy configuration validator is configured by
[`../config/default.rvre`](../config/default.rvre). Build and run
`rvre_experiment`, optionally passing a path to another `.rvre` configuration.
The executable validates configuration only and does not price options.

The The `rvre_rough_bergomi_example` executable simulates asset, variance, and
Volterra paths using a fixed seed. It writes CSV rows to stdout by default or
to a requested output file:

```sh
./build/debug/rvre_rough_bergomi_example rough_bergomi_paths.csv 32 8
```

Diagnostics, including empirical driver correlation, forward-variance means,
and discounted-forward statistics, are written to stderr. This is a synthetic
demonstration, not a market-data experiment.

The `rvre_price`, `rvre_calibrate`, and `rvre_hedge` executables are documented
in the corresponding files under `../docs/`. Their default fixtures are
synthetic.
