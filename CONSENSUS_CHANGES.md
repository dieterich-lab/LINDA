# Consensus LINDA package changes

This package merges the previous hard- and soft-constrained LINDA variants into one source package.

## Main user-facing changes

- `runLINDA()` now has `constraints_mode`, with allowed values:
  - `"hard"` - default; uses the original hard AS constraints.
  - `"soft"` - uses AS-derived objective-function penalties and does not add AS exclusion constraints.
- `splice_effect_sign` is now exposed in both modes and defaults to `"negative"`.
- `pValThresh` defaults to `0.05` and is used only in hard mode to decide which AS-affected reactions are forced to zero.

## Preserved hard-constraint formulation

The AS hard constraints are still produced by `write_as_constraints()` in the same form:

```text
reaction_variable = 0
```

for reactions satisfying:

```text
min_fdr <= pValThresh and min_score < 0
```

## Files updated

- `R/runLINDA.R`
- `R/checkInputs.R`
- `R/prepare_bn.R`
- `R/integrate_scores_in_bn.R`
- `R/integrate_scores_in_bn_hard.R`
- `R/integrate_scores_in_bn_soft.R`
- `R/write_all_constraints.R`
- `man/runLINDA.Rd`
- `README.md`
- `NEWS.md`
- `vignettes/LINDA.Rmd`
- `doc/LINDA.Rmd`

## Additional fixes

- Fixed soft-mode handling when `as.input = NULL`.
- Made AS input parsing robust to extra columns by selecting `id`, `effect`, and `significance` explicitly.
- Repaired input validation so `top = "all"` is accepted.
- Added validation for `constraints_mode`, `pValThresh`, and `save_res`.
- Removed transient RStudio/session files from the delivered source package.
