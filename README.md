# pams: Profile Analysis via Multidimensional Scaling

<!-- badges: start -->
[![CRAN status](https://www.r-pkg.org/badges/version/pams)](https://CRAN.R-project.org/package=pams)
[![R-CMD-check](https://github.com/sekangakim/pams/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/sekangakim/pams/actions/workflows/R-CMD-check.yaml)
<!-- badges: end -->

## Overview

**PAMS** implements Profile Analysis via Multidimensional Scaling for the
identification of population-level core response profiles from cross-sectional
and longitudinal person-score data.

In a typical social-science dataset each row is a *person profile* — a vector
of scores across J related subscales. PAMS decomposes each profile into two
components:

- **Level**: the person's mean across all subscales, capturing overall
  elevation.
- **Pattern**: the ipsatized subscores (deviations from the person mean),
  capturing the shape of the profile — its peaks and valleys across
  subscales.

PAMS then uses nonmetric multidimensional scaling (via the SMACOF algorithm)
on the J × J inter-variable proximity matrix to identify a small number of
**core profiles** — the central response patterns in the population. Because
the proximity matrix is J × J rather than I × I, PAMS scales to any sample
size, unlike cluster profile analysis.

The key inferential contribution of PAMS in R is the estimation of
**bootstrap standard errors and BCa confidence intervals** for every core
profile coordinate, enabling researchers to test which subscale peaks and
valleys are statistically significant. This capability is not available in
other profile analysis methods such as cluster profile analysis (CPA) or
latent profile analysis (LPA). Both CPA and LPA also recover only level
information, whereas PAMS recovers both level and pattern.

PAMS applies equally to **cross-sectional** and **longitudinal** data. For
longitudinal data the input variables are the subscale scores stacked across
time points, and the resulting core profiles are *trajectory profiles* that
show how response patterns evolve over time.

---

## Installation

```r
# Install the development version from GitHub:
# install.packages("devtools")
devtools::install_github("sekangakim/pams")
```

---

## Quick Start

```r
library(pams)
library(smacof)

# ---- Cross-sectional example ------------------------------------------------

# Load your data (persons x subscales)
cross_data <- read.csv("Cross-sectional.csv", header = FALSE)
colnames(cross_data) <- c("OV1","NS2","VA3","LP4","PP5","SR6","VS7","GI8",
                           "CF9","NR10","NP11","NW12","VA13","PR14","AS15",
                           "ON16","PC17","MW18")

# Step 1: inspect stress across dimensionalities to choose nprofile
smacofSym(dist(t(cross_data)), ndim = 2, type = "ordinal")$stress  # 0.069
smacofSym(dist(t(cross_data)), ndim = 3, type = "ordinal")$stress  # 0.047
smacofSym(dist(t(cross_data)), ndim = 4, type = "ordinal")$stress  # 0.027
# Three-dimensional solution chosen (stress <= 0.05, Kruskal 1964)

# Step 2: run PAMS with 2,000 bootstrap samples
set.seed(1)
result <- BootSmacof(
  testdata    = cross_data,
  participant = 1:10,          # persons selected for individual assessment
  mds         = "smacof",
  type        = "ordinal",
  distance    = "euclid",
  nprofile    = 3,
  direction   = c(-1, 1, 1),  # flip dimensions to aid interpretation
  cl          = 0.95,
  nBoot       = 2000,
  testname    = colnames(cross_data)
)

# Step 3: inspect results
result$MDS$stress          # stress of original fit
result$WeightmeanR2        # mean R^2 across all persons

# Core profile coordinates with BCa CIs
round(result$MDSsummary[[1]], 3)  # Core Profile 1
round(result$MDSsummary[[2]], 3)  # Core Profile 2
round(result$MDSsummary[[3]], 3)  # Core Profile 3

# Weights and core-profile correlations for selected persons
round(result$Weight[1:10, ], 2)
```

---

## Key Function

### `BootSmacof()`

The core function of the package. It fits a nonmetric MDS solution to the
J × J inter-variable distance matrix of the input data, bootstraps the
solution to produce empirical sampling distributions of core profile
coordinates, and computes BCa confidence intervals. It also estimates
person-level weights, R-squared values, and correlations with core profiles
for all participants, with optional bootstrap CIs for a selected subset.

| Argument | Description |
|---|---|
| `testdata` | Data frame or matrix of persons × subscales |
| `participant` | Integer vector of persons for individual bootstrap assessment |
| `mds` | MDS algorithm: `"smacof"` (recommended) or `"classical"` |
| `type` | Transformation type: `"ordinal"`, `"interval"`, `"ratio"`, `"mspline"` |
| `distance` | Distance measure: `"euclid"` or `"sqeuclid"` |
| `scale` | Logical; standardise variables before analysis |
| `nprofile` | Number of core profiles (dimensions) to extract |
| `direction` | Integer vector of 1 or −1 to flip dimension signs |
| `cl` | Confidence level for BCa intervals (default `0.95`) |
| `nBoot` | Number of bootstrap samples (minimum `1000` recommended) |
| `testname` | Character vector of subscale names |
| `file` | Optional file stem for saving results as CSV |

**Returns** a named list with components:

| Component | Description |
|---|---|
| `MDS` | Original MDS fit object from smacof |
| `MDSsummary` | List of K data frames: coordinate, SE, BCa CI for each core profile |
| `MDSprofile` | List of K bootstrap coordinate matrices |
| `stresssummary` | Bootstrap summary (mean, SE, BCa CI) of smacof stress |
| `stressprofile` | Vector of 2,000 bootstrap stress values |
| `MDSR2` | R² of Di regressed on other dimensions (collinearity check) |
| `Weight` | Person weights, level, R², and core-profile correlations for all persons |
| `WeightmeanR2` | Mean R² across all persons |
| `WeightB` | Bootstrap CIs for weights of selected participants |
| `PcorrB` | Bootstrap CIs for partial correlations of selected participants |

---

## Worked Examples

Full worked examples for both cross-sectional and longitudinal data —
including dimensionality selection, direction checking, BCa CI plots,
individual person assessment, and residual assumption checks — are provided
in the package vignette:

```r
vignette("PAMS_analysis", package = "pams")
```

### Cross-sectional example
Data: Woodcock-Johnson IV Cognitive Ability Battery (n = 1,650, ages 18–35,
18 cognitive subscales). Three core profiles are extracted and labelled by
their CHC factor structure peaks and valleys:

1. **High Cognitive Processing Speed with Low Long-Term Retrieval**
2. **High Working Memory with Low Cognitive Processing Speed**
3. **High Fluid Reasoning with Low Auditory Processing**

### Longitudinal example
Data: Eating Disorder Inventory-2 (n = 1,261 female anorexic patients,
11 subscales at Pre and Post treatment). Three core trajectory profiles
describe differential treatment response patterns across symptom indicators,
with BCa CIs used to identify which symptom improvements are statistically
significant after treatment.

---

## Comparison with Related Methods

| Feature | PAMS | Cluster Profile Analysis | Latent Profile Analysis |
|---|---|---|---|
| Recovers level information | ✓ | ✓ | ✓ |
| Recovers pattern information | ✓ | ✗ | ✗ |
| Scales to any sample size | ✓ | ✗ | ✓ |
| Bootstrap CIs for coordinates | ✓ | ✗ | ✗ |
| Longitudinal extension | ✓ | Limited | Limited |

---

## Citation

If you use PAMS in your research, please cite:

> Kim, S.-K., & Kim, D. (2024). Utility of profile analysis via
> multidimensional scaling in R for the study of person response profiles
> in cross-sectional and longitudinal data. *The Quantitative Methods for
> Psychology*, *20*(3), 230–247.
> https://doi.org/10.20982/tqmp.20.3.p230

For the theoretical foundation of PAMS please also cite:

> Davison, M. L. (1996). *Multidimensional scaling interest and aptitude
> profiles: Idiographic dimensions, nomothetic factors*. Presidential
> address to Division 5, American Psychological Association, Toronto.

For the SMACOF algorithm used in estimation:

> de Leeuw, J., & Mair, P. (2009). Multidimensional scaling using
> majorization: SMACOF in R. *Journal of Statistical Software*, *31*(3),
> 1–30. https://doi.org/10.18637/jss.v031.i03

---

## Related Package

[**SEPA**](https://github.com/sekangakim/sepa) (Subprofile Extraction via
Pattern Analysis) is a companion package by the same authors. While PAMS
identifies population-level core profiles of response patterns, SEPA
locates individual profiles in a two-dimensional SVD biplot space and
quantifies each person's alignment with cognitive ability domains via
direction cosines. The two packages address complementary scientific
questions and can be used together.

---

## License

MIT © Se-Kang Kim & Donghoh Kim
