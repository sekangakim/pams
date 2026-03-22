#' Profile Analysis via Multidimensional Scaling
#'
#' @description
#' Identifies population-level core response profiles from cross-sectional or
#' longitudinal person-score data using nonmetric multidimensional scaling
#' (SMACOF algorithm). Each person profile is decomposed into a level
#' component (person mean) and a pattern component (ipsatized subscores).
#' \code{BootSmacof} fits a nonmetric MDS solution to the \eqn{J \times J}
#' inter-variable distance matrix, bootstraps the solution to generate
#' empirical sampling distributions of core profile coordinates, and computes
#' bias-corrected and accelerated (BCa) confidence intervals for each
#' coordinate. Person-level weights, R-squared values, and correlations with
#' core profiles are estimated for all participants, with optional bootstrap
#' confidence intervals for a selected subset.
#'
#' @param testdata A data frame or matrix of persons (rows) by subscales
#'   (columns). Subscores are assumed to be related and continuous. For
#'   longitudinal data, columns should be ordered as all subscales at Time 1
#'   followed by all subscales at Time 2, and so on.
#' @param participant An integer vector of row indices identifying persons for
#'   whom individual bootstrap confidence intervals on weights and partial
#'   correlations are computed. If \code{NULL} (the default), individual
#'   bootstrapping is skipped and only population-level results are returned.
#' @param mds Character string specifying the MDS algorithm. Either
#'   \code{"smacof"} (default, recommended; uses the majorization algorithm
#'   of de Leeuw & Mair, 2009) or \code{"classical"} (Torgerson's classical
#'   metric MDS via \code{\link[stats]{cmdscale}}).
#' @param type Character string specifying the optimal scaling transformation
#'   passed to \code{\link[smacof]{smacofSym}}. One of \code{"ordinal"}
#'   (default, nonmetric; recommended for most social-science data),
#'   \code{"interval"}, \code{"ratio"}, or \code{"mspline"}. Ignored when
#'   \code{mds = "classical"}.
#' @param distance Character string specifying the distance measure used to
#'   compute the \eqn{J \times J} inter-variable proximity matrix. Either
#'   \code{"euclid"} (default, Euclidean distance) or \code{"sqeuclid"}
#'   (squared Euclidean distance). Note that squaring amplifies large
#'   distances and compresses small ones; \code{"euclid"} is recommended
#'   unless faster convergence is specifically required.
#' @param scale Logical. If \code{TRUE}, columns of \code{testdata} are
#'   standardised (zero mean, unit variance) before analysis. Set to
#'   \code{TRUE} when subscales have different measurement units.
#'   Default is \code{FALSE}.
#' @param nprofile A positive integer specifying the number of core profiles
#'   (MDS dimensions) to extract. Choose by inspecting stress values across
#'   2-, 3-, and 4-dimensional solutions; Kruskal's (1964) criterion of
#'   stress \eqn{\leq 0.05} is recommended. Must be less than the number of
#'   subscales (columns) in \code{testdata}.
#' @param direction An integer vector of length \code{nprofile}, with each
#'   element either \code{1} or \code{-1}. Multiplying a dimension by
#'   \code{-1} flips its sign to aid substantive interpretation (e.g., so
#'   that the first core profile aligns with the subscale mean profile).
#'   Default is \code{rep(1, nprofile)} (no flipping).
#' @param cl Numeric confidence level for BCa intervals. Default is
#'   \code{0.95}. Common alternatives are \code{0.99} and \code{0.90}.
#' @param nBoot A positive integer specifying the number of bootstrap
#'   samples. A minimum of \code{1000} is recommended for stable confidence
#'   interval estimation (Efron & Tibshirani, 1993); \code{2000} is the
#'   default.
#' @param testname An optional character vector of length equal to the number
#'   of columns in \code{testdata}, giving subscale names used as row labels
#'   in summary output and plots. If \code{NULL}, labels \code{"T1"},
#'   \code{"T2"}, \ldots are generated automatically.
#' @param file An optional character string giving a file path stem. If
#'   supplied, three CSV files are written: \code{<file>MDS.csv} (stress
#'   summary and core profile coordinates with BCa CIs),
#'   \code{<file>Weight.csv} (person weights, levels, R-squared values, and
#'   core-profile correlations), and \code{<file>WeightB.csv} (bootstrap
#'   summaries for selected participants). If \code{NULL} (the default), no
#'   files are written.
#'
#' @return A named list with the following components:
#'   \describe{
#'     \item{\code{MDS}}{The MDS fit object for the original data. When
#'       \code{mds = "smacof"} this is the full object returned by
#'       \code{\link[smacof]{smacofSym}}, including \code{$conf} (core
#'       profile coordinate matrix, \eqn{J \times K}) and \code{$stress}.
#'       When \code{mds = "classical"} this is a minimal list with
#'       \code{$conf} only.}
#'     \item{\code{MDSsummary}}{A list of \code{nprofile} data frames, one
#'       per core profile. Each data frame has rows corresponding to
#'       subscales and columns: \code{Ori} (original coordinate),
#'       \code{Mean} (bootstrap mean), \code{SE} (bootstrap standard error),
#'       \code{Lower} and \code{Upper} (percentile CI bounds),
#'       \code{BCaLower} and \code{BCaUpper} (BCa CI bounds). Coordinates
#'       whose BCa CI does not include zero are statistically significant.}
#'     \item{\code{MDSprofile}}{A list of \code{nprofile} matrices, each of
#'       dimension \code{nBoot} \eqn{\times} \eqn{J}, containing the full
#'       bootstrap distribution of core profile coordinates.}
#'     \item{\code{stresssummary}}{A one-row data frame with bootstrap
#'       summary statistics for the smacof stress value: \code{Ori},
#'       \code{Mean}, \code{SE}, \code{Lower}, \code{Upper},
#'       \code{BCaLower}, \code{BCaUpper}. \code{NULL} when
#'       \code{mds = "classical"}.}
#'     \item{\code{stressprofile}}{A numeric vector of length \code{nBoot}
#'       containing bootstrap stress values. \code{NULL} when
#'       \code{mds = "classical"}.}
#'     \item{\code{MDSR2}}{A numeric vector of length \code{nprofile}
#'       containing the R-squared values from regressing each core profile
#'       dimension on the remaining dimensions. Low values confirm that the
#'       core profiles are not collinear.}
#'     \item{\code{Weight}}{A matrix of dimension \eqn{I \times (2K + 2)}
#'       containing, for every person: raw weights (\code{w1}, \ldots,
#'       \code{wK}), level estimate, R-squared value, and correlations with
#'       each core profile (\code{corDim1}, \ldots, \code{corDimK}).
#'       Row names are \code{"#1"}, \code{"#2"}, \ldots}
#'     \item{\code{WeightmeanR2}}{The mean R-squared value across all
#'       \eqn{I} persons, summarising how well the \code{nprofile} core
#'       profiles account for pattern variance in the sample.}
#'     \item{\code{WeightB}}{A matrix of bootstrap summary statistics
#'       (original estimate, mean, SE, lower and upper CI bounds) for the
#'       weights of each person in \code{participant}. \code{NULL} if
#'       \code{participant} is \code{NULL}.}
#'     \item{\code{PcorrB}}{A matrix of bootstrap summary statistics for
#'       the partial correlations of each person in \code{participant}.
#'       \code{NULL} if \code{participant} is \code{NULL}.}
#'     \item{\code{nprofile}}{The number of core profiles extracted.}
#'     \item{\code{nBoot}}{The number of bootstrap samples used.}
#'     \item{\code{scale}}{Logical; whether columns were standardised.}
#'     \item{\code{testname}}{Character vector of subscale names used.}
#'   }
#'
#' @references
#' Davison, M. L. (1996). \emph{Multidimensional scaling interest and
#' aptitude profiles: Idiographic dimensions, nomothetic factors}.
#' Presidential address to Division 5, American Psychological Association,
#' Toronto.
#'
#' de Leeuw, J., & Mair, P. (2009). Multidimensional scaling using
#' majorization: SMACOF in R. \emph{Journal of Statistical Software},
#' \emph{31}(3), 1--30. \doi{10.18637/jss.v031.i03}
#'
#' Efron, B., & Tibshirani, R. J. (1993). \emph{An introduction to the
#' bootstrap}. Chapman & Hall.
#'
#' Kim, S.-K., & Kim, D. (2024). Utility of profile analysis via
#' multidimensional scaling in R for the study of person response profiles
#' in cross-sectional and longitudinal data. \emph{The Quantitative Methods
#' for Psychology}, \emph{20}(3), 230--247.
#' \doi{10.20982/tqmp.20.3.p230}
#'
#' Kruskal, J. B. (1964). Multidimensional scaling by optimizing goodness
#' of fit to a nonmetric hypothesis. \emph{Psychometrika}, \emph{29},
#' 1--27. \doi{10.1007/BF02289565}
#'
#' @seealso \code{\link[smacof]{smacofSym}} for the underlying MDS
#'   algorithm.
#'
#' @examples
#' \dontrun{
#' # Cross-sectional example using bundled WJ-IV cognitive ability data.
#' library(pams)
#' library(smacof)
#'
#' cross_data <- read.csv(
#'   system.file("extdata", "Cross-sectional.csv", package = "pams"),
#'   header = FALSE
#' )
#' colnames(cross_data) <- c(
#'   "OV1","NS2","VA3","LP4","PP5","SR6","VS7","GI8",
#'   "CF9","NR10","NP11","NW12","VA13","PR14","AS15",
#'   "ON16","PC17","MW18"
#' )
#'
#' # Inspect stress across dimensionalities to choose nprofile
#' smacofSym(dist(t(cross_data)), ndim = 2, type = "ordinal")$stress
#' smacofSym(dist(t(cross_data)), ndim = 3, type = "ordinal")$stress
#' smacofSym(dist(t(cross_data)), ndim = 4, type = "ordinal")$stress
#'
#' # Run PAMS with 3 core profiles and 2,000 bootstrap samples
#' set.seed(1)
#' result <- BootSmacof(
#'   testdata    = cross_data,
#'   participant = 1:10,
#'   mds         = "smacof",
#'   type        = "ordinal",
#'   distance    = "euclid",
#'   nprofile    = 3,
#'   direction   = c(-1, 1, 1),
#'   cl          = 0.95,
#'   nBoot       = 2000,
#'   testname    = colnames(cross_data)
#' )
#'
#' result$MDS$stress       # should be <= 0.05
#' result$WeightmeanR2     # mean R^2 across all persons
#'
#' # Core profile coordinates with BCa CIs
#' round(result$MDSsummary[[1]], 3)
#' round(result$MDSsummary[[2]], 3)
#' round(result$MDSsummary[[3]], 3)
#'
#' # Weights and correlations for first 10 persons
#' round(result$Weight[1:10, ], 2)
#' }
#'
#' @importFrom smacof smacofSym
#' @importFrom stats cmdscale dist cor lm quantile qnorm pnorm sd
#'   update coefficients as.formula
#' @importFrom utils write.table
#'
#' @export
BootSmacof <- function(testdata, participant = NULL,
                       mds      = c("smacof", "classical"),
                       type     = c("ratio", "interval", "ordinal", "mspline"),
                       distance = c("euclid", "sqeuclid"),
                       scale    = FALSE,
                       nprofile  = 3,
                       direction = rep(1, nprofile),
                       cl        = 0.95,
                       nBoot     = 2000,
                       testname  = NULL,
                       file      = NULL)
{
    ntest    <- ncol(testdata)
    nsubject <- nrow(testdata)
    lalpha   <- (1 - cl) / 2
    ualpha   <- 1 - lalpha

    if (is.null(testname))
        testname <- paste0("T", 1:ntest)
    if (ntest < nprofile)
        stop("The number of profiles must be less than the number of test.")
    if (length(direction) != nprofile)
        stop("The number of profile direction must be the same to the number of profile.")
    if (!(any(distance == c("euclid", "sqeuclid"))))
        stop("Distance for smacof must be one of euclid, sqeuclid.")

    if (scale) testdata <- scale(testdata)

    if (distance == "euclid")
        distance0 <- dist(t(testdata))
    else if (distance == "sqeuclid")
        distance0 <- dist(t(testdata))^2

    # ------------------------------------------------------------------
    # Fit MDS to original sample
    # ------------------------------------------------------------------

    if (mds == "smacof") {
        MDS <- smacofSym(distance0, ndim = nprofile, type = type)
        for (i in 1:nprofile)
            MDS$conf[, i] <- MDS$conf[, i] * direction[i]
        profileOri <- MDS$conf
        stressOri  <- MDS$stress
    }

    if (mds == "classical") {
        MDS        <- NULL
        MDS$conf   <- cmdscale(distance0, k = nprofile)
        MDS$stress <- stressOri <- NULL
        colnames(MDS$conf) <- paste0("D", 1:nprofile)
        for (i in 1:nprofile)
            MDS$conf[, i] <- MDS$conf[, i] * direction[i]
        profileOri <- MDS$conf
    }

    # ------------------------------------------------------------------
    # Bootstrap: empirical distribution of profiles
    # ------------------------------------------------------------------

    profileBoot <- vector("list", nprofile)
    stressBoot  <- NULL

    for (i in 1:nBoot) {
        testdataBoot <- testdata[sample(1:nsubject, replace = TRUE), ]

        if (distance == "euclid")
            distanceB <- dist(t(testdataBoot))
        else if (distance == "sqeuclid")
            distanceB <- dist(t(testdataBoot))^2

        if (mds == "smacof") {
            tmp0       <- smacofSym(distanceB, ndim = nprofile, type = type)
            tmp1       <- tmp0$conf
            stressBoot <- c(stressBoot, tmp0$stress)
        } else if (mds == "classical") {
            tmp1 <- cmdscale(distanceB, k = nprofile)
        }

        # Procrustes sign-alignment to original profile
        tmp1 <- tmp1 %*% diag(sign(diag(cor(tmp1, profileOri))))

        for (j in 1:nprofile)
            profileBoot[[j]] <- rbind(profileBoot[[j]], tmp1[, j])
    }

    # ------------------------------------------------------------------
    # BCa acceleration constants via jackknife
    # (Fix 1: inner loop variable renamed from 'i' to 'k' to avoid
    #  shadowing the outer subject-index loop variable.)
    #
    # BCa implementation based on "bootBCa" R package Version 1.0.
    # Author: S original, from StatLib, by Rob Tibshirani. R port by
    # Friedrich Leisch. Enhancements by David Flater.
    # License: BSD_3_clause + file LICENSE
    # ------------------------------------------------------------------

    profileu <- uu <- vector("list", nprofile)
    stressu  <- suu <- NULL

    for (i in 1:nsubject) {

        if (distance == "euclid")
            distance0 <- dist(t(testdata[-i, ]))
        else if (distance == "sqeuclid")
            distance0 <- dist(t(testdata[-i, ]))^2

        if (mds == "smacof") {
            tmp <- smacofSym(distance0, ndim = nprofile, type = type)
            for (k in 1:nprofile)
                tmp$conf[, k] <- tmp$conf[, k] * direction[k]
            stressu <- c(stressu, tmp$stress)
        }

        if (mds == "classical") {
            tmp      <- NULL
            tmp$conf <- cmdscale(distance0, k = nprofile)
            for (k in 1:nprofile)
                tmp$conf[, k] <- tmp$conf[, k] * direction[k]
        }

        for (j in 1:nprofile)
            profileu[[j]] <- rbind(profileu[[j]], tmp$conf[, j])
    }

    acc <- NULL
    for (j in 1:nprofile) {
        uu[[j]] <- -sweep(profileu[[j]], 2, apply(profileu[[j]], 2, mean))
        acc     <- cbind(acc,
                         apply(uu[[j]], 2,
                               function(x) sum(x*x*x) / (6 * (sum(x*x))^1.5)))
    }

    # ------------------------------------------------------------------
    # Compute BCa confidence intervals for profile coordinates
    # ------------------------------------------------------------------

    zalpha <- qnorm(c(lalpha, ualpha))
    BCACI  <- vector("list", nprofile)

    for (i in 1:nprofile) {
        tmp <- NULL
        for (j in 1:ntest) {
            z0   <- qnorm(sum(profileBoot[[i]][, j] < profileOri[, i][j]) / nBoot)
            tt   <- pnorm(z0 + (z0 + zalpha) / (1 - acc[j, i] * (z0 + zalpha)))
            tmp1 <- quantile(profileBoot[[i]][, j], probs = tt)
            tmp  <- rbind(tmp, tmp1)
        }
        BCACI[[i]] <- tmp
    }

    # ------------------------------------------------------------------
    # Summary statistics for profile coordinates
    # (Fix 2: column names standardised to 'BCaLower' / 'BCaUpper'.)
    # ------------------------------------------------------------------

    profile <- vector("list", nprofile)
    for (i in 1:nprofile)
        profile[[i]] <- data.frame(
            Ori      = profileOri[, i],
            Mean     = apply(profileBoot[[i]], 2, mean),
            SE       = apply(profileBoot[[i]], 2, sd),
            Lower    = apply(profileBoot[[i]], 2, quantile, lalpha),
            Upper    = apply(profileBoot[[i]], 2, quantile, ualpha),
            BCaLower = BCACI[[i]][, 1],
            BCaUpper = BCACI[[i]][, 2],
            row.names = testname
        )

    # Fix 3: condition corrected from 'nprofile > 2' to 'nprofile > 1'
    # so that the second profile is included when nprofile == 2.
    tmp2 <- profile[[1]]
    if (nprofile > 1)
        for (i in 2:nprofile)
            tmp2 <- cbind(tmp2, profile[[i]])

    # ------------------------------------------------------------------
    # Summary statistics for smacof stress (including BCa CI)
    # ------------------------------------------------------------------

    stresssummary <- NULL
    if (mds == "smacof") {
        suu  <- -(stressu - mean(stressu))
        sacc <- sum(suu*suu*suu) / (6 * (sum(suu*suu))^1.5)
        z0   <- qnorm(sum(stressBoot < stressOri) / nBoot)
        tt   <- pnorm(z0 + (z0 + zalpha) / (1 - sacc * (z0 + zalpha)))
        BCACIstress <- quantile(stressBoot, probs = tt)

        stresssummary <- data.frame(
            Ori      = stressOri,
            Mean     = mean(stressBoot),
            SE       = sd(stressBoot),
            Lower    = quantile(stressBoot, lalpha),
            Upper    = quantile(stressBoot, ualpha),
            BCaLower = BCACIstress[1],
            BCaUpper = BCACIstress[2],
            row.names = NULL
        )
    }

    # ------------------------------------------------------------------
    # Write MDS results to CSV (if file path supplied)
    # ------------------------------------------------------------------

    if (!is.null(file)) {
        cat("Summary Statistics for Stress",
            file = paste0(file, "MDS.csv"), sep = "", "\n")
        cat(c("Ori", "Mean", "SE", "Lower", "Upper", "BCaLower", "BCaUpper"),
            file = paste0(file, "MDS.csv"), sep = ",", "\n", append = TRUE)

        if (mds == "classical")
            cat(stressOri,
                file = paste0(file, "MDS.csv"), sep = ",", "\n\n", append = TRUE)
        if (mds == "smacof")
            cat(as.numeric(stresssummary),
                file = paste0(file, "MDS.csv"), sep = ",", "\n\n", append = TRUE)

        cat("Summary Statistics for Profile",
            file = paste0(file, "MDS.csv"), sep = "", "\n", append = TRUE)
        cat(c("Name",
              paste0(rep(c("Ori", "Mean", "SE", "Lower", "Upper",
                           "BCaLower", "BCaUpper"),
                         nprofile),
                     rep(1:nprofile, each = 7))),
            file = paste0(file, "MDS.csv"), sep = ",", "\n", append = TRUE)
        write.table(tmp2,
                    file = paste0(file, "MDS.csv"),
                    sep = ",", append = TRUE, col.names = FALSE)
    }

    # ------------------------------------------------------------------
    # Person weights and partial correlations (regression step)
    # ------------------------------------------------------------------

    formulaall  <- as.formula(paste("y ~ -1 +",
                                    paste(paste0("D", 1:nprofile), collapse = "+")))
    formulaeach <- NULL
    for (k in 1:nprofile)
        formulaeach <- c(formulaeach,
                         paste(paste0("D", k, " ~ -1+"),
                               paste(paste0("D", setdiff(1:nprofile, k),
                                            collapse = "+"))))

    R2            <- rep(0, nprofile)
    outDresiduals <- NULL
    for (j in 1:nprofile) {
        outD          <- lm(formulaeach[[j]], data.frame(profileOri))
        outDresiduals <- cbind(outDresiduals, outD$residuals)
        R2[j]         <- summary(outD)$r.squared
    }

    result <- NULL
    pcorr  <- rep(0, nprofile)
    for (i in 1:nsubject) {
        y          <- as.numeric(testdata[i, ])
        my         <- mean(y)
        y          <- y - my
        out        <- lm(formulaall, data.frame(y, profileOri))
        coeffmulti <- out$coefficients
        R2multi    <- summary(out)$r.squared

        for (j in 1:nprofile) {
            outy     <- update(out, as.formula(paste0(". ~ . -D", j)))
            pcorr[j] <- cor(outy$residuals, outDresiduals[, j])
        }
        result <- rbind(result, c(coeffmulti, my, R2multi, pcorr))
    }
    meanR2 <- mean(result[, nprofile + 2])

    colnames(result) <- c(paste0("w", 1:nprofile),
                          "level", "R^2",
                          paste0("corDim", 1:nprofile))
    rownames(result) <- paste0("#", 1:nsubject)

    if (!is.null(file)) {
        cat(c("Overall R^2", meanR2),
            file = paste0(file, "Weight.csv"), sep = ",", "\n")
        cat(c("id", paste0("w", 1:nprofile), "level", "R^2",
              paste0("corDim", 1:nprofile)),
            file = paste0(file, "Weight.csv"), sep = ",", "\n", append = TRUE)
        write.table(result,
                    file = paste0(file, "Weight.csv"),
                    sep = ",", append = TRUE, col.names = FALSE)
    }

    # ------------------------------------------------------------------
    # Bootstrap CIs for person weights and partial correlations
    # (for specified participants only)
    # ------------------------------------------------------------------

    formulaallB  <- as.formula(paste("y ~ -1 +",
                                     paste(paste0("X", 1:nprofile), collapse = "+")))
    formulaeachB <- NULL
    for (k in 1:nprofile)
        formulaeachB <- c(formulaeachB,
                          paste(paste0("X", k, " ~ -1+"),
                                paste(paste0("X", setdiff(1:nprofile, k),
                                             collapse = "+"))))

    resultB <- resultBP <- NULL

    if (!is.null(participant)) {
        sumstat <- rep(0, 5 * nprofile)
        index1  <- seq(1, 5 * nprofile, by = 5)
        index2  <- seq(2, 5 * nprofile, by = 5)
        index3  <- seq(3, 5 * nprofile, by = 5)
        index4  <- seq(4, 5 * nprofile, by = 5)
        index5  <- seq(5, 5 * nprofile, by = 5)

        xBoot <- vector("list", nBoot)
        for (k in 1:nBoot)
            for (j in 1:nprofile)
                xBoot[[k]] <- cbind(xBoot[[k]], profileBoot[[j]][k, ])

        outDresiduals <- vector("list", nBoot)
        for (k in 1:nBoot)
            for (j in 1:nprofile) {
                outD               <- lm(formulaeachB[[j]],
                                         data.frame(y, xBoot[[k]]))
                outDresiduals[[k]] <- cbind(outDresiduals[[k]], outD$residuals)
            }

        for (i in 1:length(participant)) {
            y  <- as.numeric(testdata[participant[i], ])
            my <- mean(y)
            y  <- y - my

            tmpweight <- tmppcorr <- NULL
            pcorr     <- rep(0, nprofile)

            for (k in 1:nBoot) {
                tmplm     <- lm(formulaallB, data.frame(y, xBoot[[k]]))
                tmpweight <- rbind(tmpweight, coefficients(tmplm))

                for (j in 1:nprofile) {
                    outy     <- update(tmplm, as.formula(paste0(". ~ . -X", j)))
                    pcorr[j] <- cor(outy$residuals, outDresiduals[[k]][, j])
                }
                tmppcorr <- rbind(tmppcorr, pcorr)
            }

            sumstat[index1] <- result[participant[i], 1:nprofile]
            sumstat[index2] <- apply(tmpweight, 2, mean)
            sumstat[index3] <- apply(tmpweight, 2, sd)
            sumstat[index4] <- apply(tmpweight, 2, quantile, lalpha)
            sumstat[index5] <- apply(tmpweight, 2, quantile, ualpha)
            resultB <- rbind(resultB,
                             c(sumstat, result[participant[i], -(1:nprofile)]))

            sumstat[index1] <- result[participant[i], -(1:(nprofile + 2))]
            sumstat[index2] <- apply(tmppcorr, 2, mean)
            sumstat[index3] <- apply(tmppcorr, 2, sd)
            sumstat[index4] <- apply(tmppcorr, 2, quantile, lalpha)
            sumstat[index5] <- apply(tmppcorr, 2, quantile, ualpha)
            resultBP <- rbind(resultBP, sumstat)
        }

        colnames(resultB) <- c(paste0(rep(c("w", "m", "se", "L", "U"),
                                          nprofile),
                                      rep(1:nprofile, each = 5)),
                               "level", "R^2",
                               paste0("corDim", 1:nprofile))
        rownames(resultB) <- paste0("#", participant)

        if (!is.null(file)) {
            cat(c("id",
                  paste0(rep(c("w", "m", "se", "L", "U"), nprofile),
                         rep(1:nprofile, each = 5)),
                  "level", "R^2",
                  paste0("corDim", 1:nprofile)),
                file = paste0(file, "WeightB.csv"), sep = ",", "\n")
            write.table(resultB,
                        file = paste0(file, "WeightB.csv"),
                        sep = ",", append = TRUE, col.names = FALSE)
        }

        colnames(resultBP) <- paste0(rep(c("corDim", "m", "se", "L", "U"),
                                         nprofile),
                                     rep(1:nprofile, each = 5))
        rownames(resultBP) <- paste0("#", participant)

        if (!is.null(file)) {
            cat(c("id",
                  paste0(rep(c("corDim", "m", "se", "L", "U"), nprofile),
                         rep(1:nprofile, each = 5))),
                file = paste0(file, "PcorrB.csv"), sep = ",", "\n")
            write.table(resultBP,
                        file = paste0(file, "PcorrB.csv"),
                        sep = ",", append = TRUE, col.names = FALSE)
        }
    }

    # ------------------------------------------------------------------
    # Return all results as a named list
    # ------------------------------------------------------------------

    list(MDS           = MDS,
         MDSsummary    = profile,
         MDSprofile    = profileBoot,
         stresssummary = stresssummary,
         stressprofile = stressBoot,
         MDSR2         = R2,
         Weight        = result,
         WeightmeanR2  = meanR2,
         WeightB       = resultB,
         PcorrB        = resultBP,
         nprofile      = nprofile,
         nBoot         = nBoot,
         scale         = scale,
         testname      = testname)
}
