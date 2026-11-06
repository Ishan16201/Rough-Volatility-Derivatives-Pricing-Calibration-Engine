# Foundation quick start

RVRE requires CMake 3.25 or newer and a C++20 compiler. The core library has no
Python runtime dependency. CMake fetches Catch2 3.8.1 at a pinned commit when
tests are enabled; a network connection is needed on the first configure unless
the dependency has already been populated.

```sh
cmake --preset debug
cmake --build --preset debug
ctest --preset debug
./build/debug/rvre_experiment
```

The `release` preset builds optimized code. The `sanitizers` preset enables
AddressSanitizer and UndefinedBehaviorSanitizer on supported Clang/GNU
toolchains. All presets enable tests.

## Initial configuration format

The `.rvre` file is a strict UTF-8-compatible text file with one `key = value`
per line; `#` starts a comment. Duplicate, unknown, missing, malformed, and
out-of-domain values are errors. Rates and volatility are annualized, maturity
is in years, and the asset price is in the caller's currency units.

The initial schema contains `experiment_name`, `initial_asset`, `volatility`,
`risk_free_rate`, `dividend_yield`, `maturity_years`, `seed`, `paths`, and
`time_steps`. `rvre_experiment` validates this legacy configuration only; its
placeholder result has no price estimate and must not be interpreted as a
simulation. Use the dedicated pricing, calibration and hedging executables for
those workflows.
