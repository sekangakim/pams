# =============================================================================
# PAMS_analysis.R
# Profile Analysis via Multidimensional Scaling (PAMS) in R
# Authors: Se-Kang Kim & Donghoh Kim
#
# This script demonstrates PAMS for both cross-sectional and longitudinal data.
# It requires BootSmacof.R (the PAMS core function) to be in the same directory.
#
# Sections
#   0. Setup
#   ---- CROSS-SECTIONAL PAMS -----------------------------------------------
#   1. Cross-Sectional Data: Pre-analysis and direction check
#   2. Cross-Sectional Data: Bootstrap PAMS (Parts 1-4)
#   3. Cross-Sectional Data: Assumption checks (normality & homoscedasticity)
#   ---- LONGITUDINAL PAMS --------------------------------------------------
#   4. Longitudinal Data: Pre-analysis and direction check
#   5. Longitudinal Data: Bootstrap PAMS (Parts 1-4)
#   6. Longitudinal Data: Assumption checks (normality & homoscedasticity)
# =============================================================================


# =============================================================================
# SECTION 0 : Setup
# =============================================================================

library(smacof)
library(lmtest)   # Breusch-Pagan test (Sections 3 & 6)

# Set to the folder containing BootSmacof.R and your data files
setwd("/Users/sekangkim/Documents/J_Sub/pams/")

# Load the PAMS core function
source("BootSmacof.R")


# =============================================================================
# SECTION 1 : Cross-Sectional Data - Pre-analysis and direction check
# =============================================================================

#### Load data
cross_data <- read.csv(file = "Cross-sectional.csv", header = FALSE)
colnames(cross_data) <- c("OV1","NS2","VA3","LP4","PP5","SR6","VS7","GI8",
                           "CF9","NR10","NP11","NW12","VA13","PR14","AS15",
                           "ON16","PC17","MW18")
cross_testname <- colnames(cross_data)

#### Choose number of core profiles
# Run 2-, 3-, and 4-dimensional solutions and inspect stress values.
# Kruskal (1964a) recommends stress <= 0.05.
cross_nprofile <- 3

#### Initial MDS fit (no bootstrapping) to inspect profiles and choose directions
cross_preout <- smacofSym(dist(t(cross_data)), cross_nprofile, type = "ordinal")
cross_preout$stress   # check stress value

# Plot raw dimensions to determine which (if any) need to be flipped
par(mfrow = c(1, cross_nprofile), mar = c(2.9, 4.1, 2, 1))
for (i in 1:cross_nprofile)
    plot(scale(cross_preout$conf[, i]),
         main  = paste("Dim", i, "Core Profile"),
         xlab  = "", ylab = "Coordinate", type = "b")

# --- Optional: compare smacof solution with classical MDS (Torgerson/ProfileR) ---
# install.packages("profileR")
# library(profileR)
# result0 <- pams(cross_data, cross_nprofile)
# result0$dimensional.configuration
# for (i in 1:cross_nprofile)
#     plot(scale(cross_preout$conf[, i]),
#          scale(result0$dimensional.configuration[, i]),
#          xlab = "smacof", ylab = "ProfileR")

# Classical MDS (cmdscale) is identical to ProfileR
# preout2 <- cmdscale(dist(t(cross_data)), k = cross_nprofile)
# for (i in 1:cross_nprofile)
#     plot(scale(preout2[, i]),
#          -scale(result0$dimensional.configuration[, i]),
#          xlab = "Classical", ylab = "ProfileR",
#          main = paste("Corr =",
#                       round(cor(scale(preout2[, i]),
#                                 -scale(result0$dimensional.configuration[, i])), 3)))
# ---------------------------------------------------------------------------------

#### Set flip directions after inspecting the plots above
# -1 flips a dimension; 1 leaves it unchanged.
cross_direction <- c(-1, 1, 1)

for (i in 1:cross_nprofile)
    cross_preout$conf[, i] <- cross_preout$conf[, i] * cross_direction[i]

# Replot with corrected directions
par(mfrow = c(1, cross_nprofile), mar = c(2.9, 4.1, 2, 1))
for (i in 1:cross_nprofile)
    plot(scale(cross_preout$conf[, i]),
         main  = paste("Dim", i, "Core Profile"),
         xlab  = "", ylab = "Coordinate", type = "b")


# =============================================================================
# SECTION 2 : Cross-Sectional Data - Bootstrap PAMS (Parts 1-4)
# =============================================================================

set.seed(1)   # for reproducibility

cross_participant <- 1:10   # participants selected for individual assessment

cross_empprofile <- BootSmacof(
    testdata    = cross_data,
    participant = cross_participant,
    mds         = "smacof",
    type        = "ordinal",
    distance    = "euclid",
    scale       = FALSE,
    nprofile    = cross_nprofile,
    direction   = cross_direction,
    cl          = 0.95,
    nBoot       = 2000,
    testname    = cross_testname,
    file        = "cross_profile"   # saves cross_profileMDS.csv, cross_profileWeight.csv
)

# Save/reload to avoid re-running the bootstrap
save(cross_empprofile, file = "cross_empprofile.RData")
# load("cross_empprofile.RData")   # uncomment to skip bootstrap next session

# ---- Part 1 : MDS result and bootstrap stress -------------------------------

## Collinearity among core profiles: R^2 of Di = a*Dj + b*Dk + error
cross_empprofile$MDSR2

## Stress value of MDS fit to original data
cross_empprofile$MDS$stress

## 2-D scatter of core profile coordinates
pairs(cross_empprofile$MDS$conf)

## Correlation between the first core profile and the subscale mean profile
## (r ~ 1 indicates Dim 1 captures the g-factor structure)
cor(apply(cross_data, 2, mean), cross_empprofile$MDS$conf[, 1])

## Bootstrap stress summary and histogram
round(cross_empprofile$stresssummary, 4)

par(mfrow = c(1, 1), mar = c(2.9, 4.1, 2, 1))
hist(cross_empprofile$stressprofile, main = "Cross-Sectional: Bootstrap Stress")
abline(v = cross_empprofile$stresssummary$Ori,      lwd = 2)
abline(v = c(cross_empprofile$stresssummary$BCaLower,
             cross_empprofile$stresssummary$BCaUpper), col = "blue", lwd = 2, lty = 2)

# ---- Part 2 : Core profile coordinates with 95% BCa CIs --------------------

## Summary table for all core profiles (coordinates, SE, BCa CIs)
for (i in 1:cross_nprofile) {
    cat("\n--- Cross-Sectional Core Profile", i, "---\n")
    print(round(cross_empprofile$MDSsummary[[i]], 3))
}

## BCa CI plot for a chosen core profile (change i as needed)
i <- 3
par(mfrow = c(1, 1), mar = c(4.1, 4.1, 2, 1))
plot(cross_empprofile$MDSsummary[[i]]$Ori,
     type = "b", xlab = "", ylab = "Coordinate",
     main = paste0("Cross-Sectional Core Profile ", i, ": 95% BCa CI"),
     ylim = range(cross_empprofile$MDSsummary[[i]]$BCaLower,
                  cross_empprofile$MDSsummary[[i]]$BCaUpper))
lines(cross_empprofile$MDSsummary[[i]]$BCaLower, lty = 4, col = "blue")
lines(cross_empprofile$MDSsummary[[i]]$BCaUpper, lty = 4, col = "blue")
abline(h = 0, lty = 5, col = "red")

# ---- Part 3 : Person R^2 values and correlations with core profiles ---------

## Mean R^2 across all participants
round(cross_empprofile$WeightmeanR2, 2)

## Weights, levels, R^2, and core-profile correlations for first 10 participants
round(cross_empprofile$Weight[1:10, ], 2)

# ---- Part 4 : Individual person assessment ----------------------------------

## Compare a specific person's profile against each core profile
p <- 9
par(mfrow = c(1, cross_nprofile), mar = c(2.9, 4.1, 2, 1))
for (i in 1:cross_nprofile) {
    plot(scale(as.numeric(cross_data[p, ])),
         main = paste("Core Profile", i),
         xlab = "", ylab = "Coordinate",
         col = "blue", type = "b")
    lines(scale(cross_empprofile$MDS$conf[, i]), lty = 2)
    abline(h = 0, lty = 3)
}

## Bootstrap CIs for weights and partial correlations of person p
personnumber <- paste0("#", p)
round(cross_empprofile$WeightB[rownames(cross_empprofile$WeightB) == personnumber, ], 2)
round(cross_empprofile$PcorrB[ rownames(cross_empprofile$PcorrB)  == personnumber, ], 2)


# =============================================================================
# SECTION 3 : Cross-Sectional Data - Assumption Checks
#             (Univariate normality and homoscedasticity of PAMS residuals)
# =============================================================================

# PAMS requires normality and homogeneous variance for the error terms.
# These are evaluated via residuals from the person-level regression:
#   ipsatized person profile ~ -1 + core profile 1 + core profile 2 + core profile 3
#
# Shapiro-Wilks tests H0: residuals are normal.
# Breusch-Pagan tests H0: residual variance is homogeneous.
# For both: p-value > alpha --> H0 is not rejected (assumption holds).
#
# Note: with only J = 18 observations per person regression, individual
# test power is modest. The overall pattern across all participants matters
# more than any single test result.

cross_profileOri <- cross_empprofile$MDS$conf   # use directly from BootSmacof output
cross_nsubject   <- nrow(cross_data)
cross_testres    <- NULL

for (i in 1:cross_nsubject) {
    y  <- as.numeric(cross_data[i, ])
    my <- mean(y)
    y  <- y - my
    out <- lm(y ~ -1 + D1 + D2 + D3,
              data = data.frame(y,
                                D1 = cross_profileOri[, 1],
                                D2 = cross_profileOri[, 2],
                                D3 = cross_profileOri[, 3]))
    sw <- shapiro.test(out$residuals)$p.value
    bp <- bptest(y ~ -1 + D1 + D2 + D3,
                 data = data.frame(y,
                                   D1 = cross_profileOri[, 1],
                                   D2 = cross_profileOri[, 2],
                                   D3 = cross_profileOri[, 3]))$p.value
    cross_testres <- rbind(cross_testres, c(sw, bp))
}
colnames(cross_testres) <- c("Shapiro-Wilks p", "Breusch-Pagan p")
round(cross_testres, 3)

# Proportion of participants where H0 is rejected at alpha = 0.05
alpha <- 0.05
cat("Cross-sectional: proportion of significant Shapiro-Wilks tests:",
    round(mean(cross_testres[, 1] < alpha), 3), "\n")
cat("Cross-sectional: proportion of significant Breusch-Pagan tests:",
    round(mean(cross_testres[, 2] < alpha), 3), "\n")

## Diagnostic plots for a specific person (change i as needed)
i <- 100
y  <- as.numeric(cross_data[i, ])
my <- mean(y)
y  <- y - my
out <- lm(y ~ -1 + D1 + D2 + D3,
          data = data.frame(y,
                            D1 = cross_profileOri[, 1],
                            D2 = cross_profileOri[, 2],
                            D3 = cross_profileOri[, 3]))

par(mfrow = c(1, 2))
qqnorm(out$residuals, main = paste("Cross-Sectional: Q-Q Plot, Person", i))
qqline(out$residuals)
plot(out$fitted.values, out$residuals,
     xlab = "Predicted", ylab = "Residual",
     main = paste("Cross-Sectional: Residual Plot, Person", i))
abline(h = 0)


# =============================================================================
# SECTION 4 : Longitudinal Data - Pre-analysis and direction check
# =============================================================================

#### Load data
# Columns E1_1-E1_11 = Pre (Time 1) subscales; E3_1-E3_11 = Post (Time 2) subscales
long_data <- read.csv(file = "Longitudinal.csv", header = TRUE)
colnames(long_data) <- c("E1_1","E1_2","E1_3","E1_4","E1_5","E1_6","E1_7",
                          "E1_8","E1_9","E1_10","E1_11",
                          "E3_1","E3_2","E3_3","E3_4","E3_5","E3_6","E3_7",
                          "E3_8","E3_9","E3_10","E3_11")
long_testname <- colnames(long_data)

#### Choose number of core profiles
long_nprofile <- 3

#### Initial MDS fit to inspect profiles and choose directions
long_preout <- smacofSym(dist(t(long_data)), long_nprofile, type = "ordinal")
long_preout$stress

par(mfrow = c(1, long_nprofile), mar = c(2.9, 4.1, 2, 1))
for (i in 1:long_nprofile)
    plot(scale(long_preout$conf[, i]),
         main  = paste("Dim", i, "Scaled Profile"),
         xlab  = "", ylab = "Coordinate", type = "b")

#### Set flip directions after inspecting the plots above
long_direction <- c(1, -1, -1)

for (i in 1:long_nprofile)
    long_preout$conf[, i] <- long_preout$conf[, i] * long_direction[i]

par(mfrow = c(1, long_nprofile), mar = c(2.9, 4.1, 2, 1))
for (i in 1:long_nprofile)
    plot(scale(long_preout$conf[, i]),
         main  = paste("Dim", i, "Scaled Profile"),
         xlab  = "", ylab = "Coordinate", type = "b")


# =============================================================================
# SECTION 5 : Longitudinal Data - Bootstrap PAMS (Parts 1-4)
# =============================================================================

set.seed(1)

long_participant <- 1:10

long_empprofile <- BootSmacof(
    testdata    = long_data,
    participant = long_participant,
    mds         = "smacof",
    type        = "ordinal",
    distance    = "euclid",
    scale       = FALSE,
    nprofile    = long_nprofile,
    direction   = long_direction,
    cl          = 0.95,
    nBoot       = 2000,
    testname    = long_testname,
    file        = "long_profile"
)

save(long_empprofile, file = "long_empprofile.RData")
# load("long_empprofile.RData")

# ---- Part 1 : MDS result and bootstrap stress -------------------------------

long_empprofile$MDSR2
long_empprofile$MDS$stress
pairs(long_empprofile$MDS$conf)

## Correlation between first core profile and subscale mean profile
cor(apply(long_data, 2, mean), long_empprofile$MDS$conf[, 1])

round(long_empprofile$stresssummary, 4)

par(mfrow = c(1, 1), mar = c(2.9, 4.1, 2, 1))
hist(long_empprofile$stressprofile, main = "Longitudinal: Bootstrap Stress")
abline(v = long_empprofile$stresssummary$Ori,      lwd = 2)
abline(v = c(long_empprofile$stresssummary$BCaLower,
             long_empprofile$stresssummary$BCaUpper), col = "blue", lwd = 2, lty = 2)

# ---- Part 2 : Pre and Post core profile coordinates with 95% BCa CIs -------

for (i in 1:long_nprofile) {
    cat("\n--- Longitudinal Core Profile", i, "---\n")
    print(round(long_empprofile$MDSsummary[[i]], 3))
}

## BCa CI plot for a chosen core profile (change i as needed)
## Longitudinal plots use paired Pre (1:11) and Post (12:22) x-axis positions.
long_colors <- c("blue", "red", "darkgreen")

for (i in 1:long_nprofile) {
    par(mfrow = c(1, 1), mar = c(2.1, 4.1, 2, 1), oma = c(2, 0, 0, 0))
    plot(long_empprofile$MDSsummary[[i]]$Ori,
         type = "b", lty = 2, col = long_colors[i], lwd = 2,
         xlab = "", ylab = "Coordinate",
         main = paste0("Longitudinal Core Profile ", i, ": 95% BCa CI"),
         ylim = range(long_empprofile$MDSsummary[[i]]$BCaLower,
                      long_empprofile$MDSsummary[[i]]$BCaUpper))
    lines(long_empprofile$MDSsummary[[i]]$BCaLower, lty = 4, col = long_colors[i])
    lines(long_empprofile$MDSsummary[[i]]$BCaUpper, lty = 4, col = long_colors[i])
    abline(h = 0, lty = 5, col = long_colors[i])
    # Add Pre/Post separator line
    abline(v = 11.5, lty = 3, col = "grey50")
    mtext("Pre", side = 1, at = 6,   line = 0.5, outer = FALSE, cex = 0.85)
    mtext("Post", side = 1, at = 17, line = 0.5, outer = FALSE, cex = 0.85)
}

# ---- Part 3 : Person R^2 values and correlations with core profiles ---------

round(long_empprofile$WeightmeanR2, 2)
round(long_empprofile$Weight[1:10, ], 2)

# ---- Part 4 : Individual person assessment ----------------------------------

p <- 9
par(mfrow = c(1, long_nprofile), mar = c(2.9, 4.1, 2, 1))
for (i in 1:long_nprofile) {
    plot(scale(as.numeric(long_data[p, ])),
         main = paste("Long. Core Profile", i),
         xlab = "", ylab = "Coordinate",
         col = "blue", type = "b")
    lines(scale(long_empprofile$MDS$conf[, i]), lty = 2)
    abline(h = 0, lty = 3)
    abline(v = 11.5, lty = 3, col = "grey50")
}

personnumber <- paste0("#", p)
round(long_empprofile$WeightB[rownames(long_empprofile$WeightB) == personnumber, ], 2)
round(long_empprofile$PcorrB[ rownames(long_empprofile$PcorrB)  == personnumber, ], 2)


# =============================================================================
# SECTION 6 : Longitudinal Data - Assumption Checks
#             (Normality and homoscedasticity of PAMS residuals)
# =============================================================================

# PAMS requires normality and homogeneous variance for the error terms.
# These are evaluated identically to Section 3, using residuals from the
# person-level regression:
#   ipsatized person profile ~ -1 + core profile 1 + core profile 2 + core profile 3
#
# Shapiro-Wilks tests H0: residuals are normal.
# Breusch-Pagan tests H0: residual variance is homogeneous.
# For both: p-value > alpha --> H0 is not rejected (assumption holds).
#
# Note: with J = 22 observations per person regression (11 Pre + 11 Post),
# individual test power is modest. The overall pattern across all
# participants matters more than any single result.

long_profileOri <- long_empprofile$MDS$conf
long_nsubject   <- nrow(long_data)
long_testres    <- NULL

for (i in 1:long_nsubject) {
    y  <- as.numeric(long_data[i, ])
    my <- mean(y)
    y  <- y - my
    df_reg <- data.frame(y,
                         D1 = long_profileOri[, 1],
                         D2 = long_profileOri[, 2],
                         D3 = long_profileOri[, 3])
    out <- lm(y ~ -1 + D1 + D2 + D3, data = df_reg)
    sw  <- shapiro.test(out$residuals)$p.value
    bp  <- bptest(y ~ -1 + D1 + D2 + D3, data = df_reg)$p.value
    long_testres <- rbind(long_testres, c(sw, bp))
}
colnames(long_testres) <- c("Shapiro-Wilks p", "Breusch-Pagan p")

alpha <- 0.05
cat("Longitudinal: proportion of significant Shapiro-Wilks tests:",
    round(mean(long_testres[, 1] < alpha), 3), "\n")
cat("Longitudinal: proportion of significant Breusch-Pagan tests:",
    round(mean(long_testres[, 2] < alpha), 3), "\n")

## Diagnostic plots for a specific person (change i as needed)
i <- 9
y  <- as.numeric(long_data[i, ])
my <- mean(y)
y  <- y - my
df_reg <- data.frame(y,
                     D1 = long_profileOri[, 1],
                     D2 = long_profileOri[, 2],
                     D3 = long_profileOri[, 3])
out <- lm(y ~ -1 + D1 + D2 + D3, data = df_reg)

par(mfrow = c(1, 2))
qqnorm(out$residuals, main = paste("Longitudinal: Q-Q Plot, Person", i))
qqline(out$residuals)
plot(out$fitted.values, out$residuals,
     xlab = "Predicted", ylab = "Residual",
     main = paste("Longitudinal: Residual Plot, Person", i))
abline(h = 0)

# =============================================================================
# End of PAMS_analysis.R
# =============================================================================
