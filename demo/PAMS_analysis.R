# =============================================================================
# Reproducible PAMS analysis
# =============================================================================

library(pams)

# The built-in USArrests data are used only to demonstrate the software
# workflow. States are persons and the four variables are related measures.
example_data <- as.data.frame(USArrests[1:12, ])
example_names <- colnames(example_data)

# Step 1: preliminary MDS for dimensionality and direction inspection.
# Standardize because the variables have different measurement units.
standardized <- scale(example_data)
proximity <- dist(t(standardized))
preliminary <- smacof::smacofSym(
    proximity,
    ndim = 2,
    type = "ordinal"
)
preliminary$stress
preliminary$conf

# Inspect the coordinate plots before selecting direction. A value of -1
# reverses an axis; it does not change distances, stress, fit, or whether a
# confidence interval excludes zero.
op <- par(mfrow = c(1, 2))
for (j in 1:2) {
    plot(
        preliminary$conf[, j],
        type = "b",
        main = paste("Preliminary dimension", j),
        xlab = "Variable",
        ylab = "Coordinate",
        xaxt = "n"
    )
    axis(1, at = seq_along(example_names), labels = example_names)
    abline(h = 0, lty = 3)
}
par(op)

# Step 2: bootstrap PAMS. The small bootstrap count keeps the demo fast.
# Use at least 1,000 samples for substantive inference.
set.seed(2026)
fit <- suppressWarnings(BootSmacof(
    testdata = example_data,
    participant = 1:3,
    mds = "smacof",
    type = "ordinal",
    distance = "euclid",
    scale = TRUE,
    nprofile = 2,
    direction = c(1, 1),
    cl = 0.95,
    nBoot = 10,
    testname = example_names
))

# Step 3: standard fitted-model methods.
summary(fit)
plot(fit, profiles = 1:2, interval = "BCa")

# Existing named components remain directly accessible.
round(fit$MDSsummary[[1]], 3)
round(fit$Weight[1:5, ], 3)

# Step 4: optional tidyverse-style tabular handling.
if (requireNamespace("tibble", quietly = TRUE)) {
    tibble::as_tibble(fit$Weight, rownames = "participant")
}
