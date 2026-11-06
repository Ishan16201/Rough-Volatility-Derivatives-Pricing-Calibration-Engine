# Project instructions

- Preserve existing research documents, experiments, and user changes. Keep
  implementation changes scoped to the requested work.
- Use C++20, separate public interfaces in `include/rvre/` from implementation
  in `src/`, and keep the pricing library free of a Python runtime dependency.
- Represent financial quantities with explicit unit-bearing types. Document
  model assumptions, parameter domains, time conventions, and discretization
  choices; never imply that a placeholder or approximation is an exact model.
- Make randomness explicit and reproducible. Record seeds and numerical
  settings; do not use hidden global random state or report fabricated results.
- Surface invalid inputs and numerical failures with actionable errors; do not
  silently clamp, default, or return success-shaped results.
- Add or update focused automated tests for behavioral changes. Run configure,
  build, tests, and executable smoke checks when changing project foundations;
  report actual commands and outcomes.
- Do not claim empirical calibration or performance improvement without
  documented genuine market data and reproducible evidence.
