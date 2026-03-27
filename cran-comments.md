## Resubmission

This is a resubmission addressing the review comments from Konstanze Lauseker:

* Added DOI references to the DESCRIPTION field in the required format
  authors (year) <doi:...>
* Replaced \dontrun{} with \donttest{} in examples; also added a small
  toy example that runs automatically in < 5 sec
* File writing in examples and demo now uses tempdir() instead of the
  working directory
* Added par() save/restore (op <- par(...) / par(op)) throughout
  demo/PAMS_analysis.R; added setwd() save/restore via on.exit()

## R CMD check results

0 errors | 0 warnings | 0 notes

Tested on macOS Tahoe 26.3.1, R 4.5.2 (aarch64-apple-darwin20)
