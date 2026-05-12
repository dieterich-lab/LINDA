# LINDA 0.2.0

## New

- Added `constraints_mode` to `runLINDA()` with two values:
  - `"hard"` uses AS constraints and forces significantly AS-affected reactions to zero.
  - `"soft"` uses AS-derived reaction penalties in the objective function while keeping reactions feasible.
- Changed the default `splice_effect_sign` from `"both"` to `"negative"`.
- Restored `splice_effect_sign` as a user-facing argument in the hard-constrained workflow.

## Fixes

- Fixed documentation inconsistencies around `pValThresh`, `splice_effect_sign`, and the hard/soft AS modes.
- Fixed AS score integration when `as.input = NULL` in the soft workflow.
- Made AS input handling robust to extra columns by explicitly selecting `id`, `effect`, and `significance`.
- Allowed `top = "all"` in input validation, matching existing `bin_measurements()` behavior.
- Improved validation messages for `solverPath`, `threads`, and `save_res`.
