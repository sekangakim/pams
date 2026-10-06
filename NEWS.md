# pams 0.2.0

## User interface

- `BootSmacof()` now returns an object of class `pams_fit` while retaining all
  previously available named-list components.
- Added `plot.pams_fit()`, `summary.pams_fit()`, `print.pams_fit()`, and
  `print.summary.pams_fit()` methods.
- Corrected the documented and implemented default transformation to
  `type = "ordinal"`.

## Reliability and documentation

- Added systematic validation for data, profile count, directions, confidence
  level, bootstrap count, participant indices, labels, and choice arguments.
- Added an explicit warning when fewer than 1,000 bootstrap samples are used.
- Made BCa calculations robust to boundary bias estimates and zero jackknife
  acceleration denominators.
- Applied sign alignment to both bootstrap and jackknife configurations.
- Clarified that person weights are unstandardized no-intercept OLS
  coefficients and that `corDim` values are partial correlations.
- Added explicit documentation of sign indeterminacy and the absence of
  general rotational or dimension-permutation alignment.
- Added a worked vignette, automated tests, and a multi-platform GitHub Actions
  R CMD check workflow.

## Metadata

- Updated the maintainer address to `sekangandroid@gmail.com`.
- Added optional vignette and tidy-tabular dependencies under `Suggests`.
