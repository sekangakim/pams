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
#' coordinate. Person-level weights, R-squared values, and partial correlations with
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
#'   Inspect the preliminary \code{smacofSym()} coordinate plots and choose
#'   signs so that prespecified anchor variables appear on the desired side of
#'   each axis. Reversing a sign does not alter distances, stress, fit, or
#'   whether an interval excludes zero. Default is \code{rep(1, nprofile)}.
#' @param cl Numeric confidence level for BCa intervals. Default is
#'   \code{0.95}. Common alternatives are \code{0.99} and \code{0.90}.
#' @param nBoot A positive integer specifying the number of bootstrap
#'   samples. A minimum of \code{1000} is recommended for stable confidence
#'   interval estimation (Efron & Tibshirani, 1993); values below 1000 issue a
#'   warning and values below 10 are rejected. \code{2000} is the default.
#' @param testname An optional character vector of length equal to the number
#'   of columns in \code{testdata}, giving subscale names used as row labels
#'   in summary output and plots. If \code{NULL}, labels \code{"T1"},
#'   \code{"T2"}, \ldots are generated automatically.
#' @param file An optional character string giving a file path stem. If
#'   supplied, two CSV files are always written: \code{<file>MDS.csv} (stress
#'   summary and core profile coordinates with BCa CIs),
#'   \code{<file>Weight.csv} (person weights, levels, R-squared values, and
#'   core-profile partial correlations). When \code{participant} is not \code{NULL},
#'   \code{<file>WeightB.csv} and \code{<file>PcorrB.csv} are also written.
#'   If \code{NULL} (the default), no files are written.
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
#'       whose pointwise BCa CI does not include zero are statistically
#'       significant at the stated coordinate-wise level; these intervals are
#'       not adjusted for multiplicity.}
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
#'       \code{wK}), level estimate, R-squared value, and partial correlations
#'       with each core profile (\code{corDim1}, \ldots, \code{corDimK}).
#'       The raw weights are unstandardized no-intercept OLS coefficients from
#'       regressing the person's ipsatized pattern on the retained coordinate
#'       vectors. Each \code{corDim} value is the correlation between the
#'       residualized person pattern and residualized focal core profile after
#'       controlling for the other profiles. Row names are \code{"#1"},
#'       \code{"#2"}, \ldots}
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
#'   The returned list has class \code{"pams_fit"}, while retaining direct
#'   access to every component listed above.
#'
#' @section Sign indeterminacy and alignment:
#' MDS coordinate signs are arbitrary. For each bootstrap and jackknife
#' sample, \code{BootSmacof()} aligns each dimension to its original-sample
#' counterpart by its correlation sign and then reports the orientation
#' selected through \code{direction}. This is sign alignment only: the method
#' does not perform general rotational alignment or dimension-permutation
#' matching. Coordinate-wise inference should therefore be interpreted
#' cautiously when dimensions are weak or nearly interchangeable. A sign
#' reversal changes only the reported orientation; it does not change
#' interpoint distances, stress, model fit, or interval exclusion of zero.
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
#' # Small toy example (runs automatically)
#' set.seed(42)
#' toy_data <- as.data.frame(matrix(rnorm(8 * 4, mean = 10, sd = 2),
#'                                  nrow = 8, ncol = 4))
#' colnames(toy_data) <- paste0("S", 1:4)
#'
#' result <- suppressWarnings(BootSmacof(
#'   testdata    = toy_data,
#'   participant = NULL,
#'   mds         = "smacof",
#'   type        = "ordinal",
#'   distance    = "euclid",
#'   nprofile    = 2,
#'   direction   = c(1, 1),
#'   cl          = 0.95,
#'   nBoot       = 10,
#'   testname    = colnames(toy_data)
#' ))
#' result$MDS$stress
#' round(result$WeightmeanR2, 2)
#' round(result$MDSsummary[[1]], 3)
#'
#' @export
BootSmacof <- function(testdata, participant = NULL,
                       mds      = c("smacof", "classical"),
                       type     = c("ordinal", "interval", "ratio", "mspline"),
                       distance = c("euclid", "sqeuclid"),
                       scale    = FALSE,
                       nprofile  = 3,
                       direction = rep(1, nprofile),
                       cl        = 0.95,
                       nBoot     = 2000,
                       testname  = NULL,
                       file      = NULL)
{
    call <- match.call()
    mds <- match.arg(mds)
    type <- match.arg(type)
    distance <- match.arg(distance)

    if (!is.matrix(testdata) && !is.data.frame(testdata)) {
        stop("`testdata` must be a numeric matrix or data frame.", call. = FALSE)
    }
    if (is.data.frame(testdata) &&
        !all(vapply(testdata, is.numeric, logical(1L)))) {
        stop("Every column of `testdata` must be numeric.", call. = FALSE)
    }
    testdata <- as.matrix(testdata)
    if (!is.numeric(testdata)) {
        stop("`testdata` must contain only numeric values.", call. = FALSE)
    }
    storage.mode(testdata) <- "double"
    if (length(dim(testdata)) != 2L || nrow(testdata) < 3L || ncol(testdata) < 2L) {
        stop("`testdata` must contain at least three persons and two variables.",
             call. = FALSE)
    }
    if (anyNA(testdata) || any(!is.finite(testdata))) {
        stop("`testdata` must not contain missing or infinite values.", call. = FALSE)
    }

    nsubject <- nrow(testdata)
    ntest <- ncol(testdata)
    nprofile <- .pams_integer(nprofile, "nprofile")
    nBoot <- .pams_integer(nBoot, "nBoot", minimum = 10L)

    if (nprofile >= ntest) {
        stop("`nprofile` must be smaller than the number of variables in `testdata`.",
             call. = FALSE)
    }
    if (nBoot < 1000L) {
        warning("Fewer than 1,000 bootstrap samples may yield unstable confidence intervals.",
                call. = FALSE)
    }
    if (length(scale) != 1L || is.na(scale) || !is.logical(scale)) {
        stop("`scale` must be either TRUE or FALSE.", call. = FALSE)
    }
    if (length(cl) != 1L || !is.numeric(cl) || !is.finite(cl) || cl <= 0 || cl >= 1) {
        stop("`cl` must be a single finite number strictly between 0 and 1.",
             call. = FALSE)
    }
    if (!is.numeric(direction) || length(direction) != nprofile ||
        anyNA(direction) || any(!direction %in% c(-1, 1))) {
        stop("`direction` must contain exactly one value (1 or -1) per profile.",
             call. = FALSE)
    }
    direction <- as.integer(direction)

    if (is.null(participant) || length(participant) == 0L) {
        participant <- NULL
    } else {
        if (!is.numeric(participant) || anyNA(participant) ||
            any(!is.finite(participant)) || any(participant != floor(participant))) {
            stop("`participant` must be NULL or an integer vector of row indices.",
                 call. = FALSE)
        }
        participant <- as.integer(participant)
        if (any(participant < 1L | participant > nsubject)) {
            stop("Every `participant` index must identify a row of `testdata`.",
                 call. = FALSE)
        }
        if (anyDuplicated(participant)) {
            stop("`participant` indices must be unique.", call. = FALSE)
        }
    }

    if (is.null(testname)) {
        supplied_names <- colnames(testdata)
        testname <- if (!is.null(supplied_names) &&
                       all(!is.na(supplied_names)) &&
                       all(nzchar(supplied_names))) {
            supplied_names
        } else {
            paste0("T", seq_len(ntest))
        }
    }
    if (!is.character(testname) || length(testname) != ntest ||
        anyNA(testname) || any(!nzchar(testname)) || anyDuplicated(testname)) {
        stop("`testname` must contain one unique, non-empty label per variable.",
             call. = FALSE)
    }
    if (!is.null(file) &&
        (!is.character(file) || length(file) != 1L || is.na(file) || !nzchar(file))) {
        stop("`file` must be NULL or a single non-empty character string.",
             call. = FALSE)
    }

    variable_sd <- apply(testdata, 2L, stats::sd)
    if (any(variable_sd == 0)) {
        stop("Every variable in `testdata` must have non-zero variance.", call. = FALSE)
    }
    if (any(apply(testdata, 1L, stats::sd) == 0)) {
        warning("At least one person has no within-profile variation; its R-squared and correlations are undefined.",
                call. = FALSE)
    }
    if (scale) testdata <- base::scale(testdata)

    lalpha <- (1 - cl) / 2
    ualpha <- 1 - lalpha
    distance0 <- .pams_distance(testdata, distance)

    if (identical(mds, "smacof")) {
        MDS <- smacof::smacofSym(distance0, ndim = nprofile, type = type)
        MDS$conf <- sweep(MDS$conf, 2L, direction, `*`)
        stressOri <- MDS$stress
    } else {
        MDS <- list(
            conf = stats::cmdscale(distance0, k = nprofile),
            stress = NULL
        )
        MDS$conf <- sweep(MDS$conf, 2L, direction, `*`)
        stressOri <- NULL
    }
    profileOri <- MDS$conf
    rownames(profileOri) <- testname
    colnames(profileOri) <- paste0("D", seq_len(nprofile))
    MDS$conf <- profileOri

    profileBoot <- lapply(
        seq_len(nprofile),
        function(...) matrix(NA_real_, nrow = nBoot, ncol = ntest)
    )
    stressBoot <- if (identical(mds, "smacof")) numeric(nBoot) else NULL

    for (b in seq_len(nBoot)) {
        boot_data <- testdata[sample.int(nsubject, nsubject, replace = TRUE), , drop = FALSE]
        boot_distance <- .pams_distance(boot_data, distance)
        if (identical(mds, "smacof")) {
            boot_fit <- smacof::smacofSym(
                boot_distance,
                ndim = nprofile,
                type = type
            )
            configuration <- boot_fit$conf
            stressBoot[b] <- boot_fit$stress
        } else {
            configuration <- stats::cmdscale(boot_distance, k = nprofile)
        }
        configuration <- .align_pams_signs(configuration, profileOri)
        for (j in seq_len(nprofile)) profileBoot[[j]][b, ] <- configuration[, j]
    }

    profileJack <- lapply(
        seq_len(nprofile),
        function(...) matrix(NA_real_, nrow = nsubject, ncol = ntest)
    )
    stressJack <- if (identical(mds, "smacof")) numeric(nsubject) else NULL

    for (i in seq_len(nsubject)) {
        jack_distance <- .pams_distance(testdata[-i, , drop = FALSE], distance)
        if (identical(mds, "smacof")) {
            jack_fit <- smacof::smacofSym(
                jack_distance,
                ndim = nprofile,
                type = type
            )
            configuration <- jack_fit$conf
            stressJack[i] <- jack_fit$stress
        } else {
            configuration <- stats::cmdscale(jack_distance, k = nprofile)
        }
        configuration <- .align_pams_signs(configuration, profileOri)
        for (j in seq_len(nprofile)) profileJack[[j]][i, ] <- configuration[, j]
    }

    profile <- vector("list", nprofile)
    for (j in seq_len(nprofile)) {
        values <- t(vapply(
            seq_len(ntest),
            function(k) .pams_summary_row(
                original = profileOri[k, j],
                bootstrap = profileBoot[[j]][, k],
                jackknife = profileJack[[j]][, k],
                lower = lalpha,
                upper = ualpha
            ),
            numeric(7L)
        ))
        profile[[j]] <- as.data.frame(values)
        rownames(profile[[j]]) <- testname
    }
    names(profile) <- paste0("Profile", seq_len(nprofile))
    names(profileBoot) <- names(profile)

    stresssummary <- NULL
    if (identical(mds, "smacof")) {
        stresssummary <- as.data.frame(t(.pams_summary_row(
            original = stressOri,
            bootstrap = stressBoot,
            jackknife = stressJack,
            lower = lalpha,
            upper = ualpha
        )))
        rownames(stresssummary) <- NULL
    }

    R2 <- .pams_profile_r2(profileOri)
    names(R2) <- paste0("D", seq_len(nprofile))

    result <- matrix(
        NA_real_,
        nrow = nsubject,
        ncol = 2L * nprofile + 2L,
        dimnames = list(
            paste0("#", seq_len(nsubject)),
            c(paste0("w", seq_len(nprofile)), "level", "R^2",
              paste0("corDim", seq_len(nprofile)))
        )
    )
    for (i in seq_len(nsubject)) {
        observed <- as.numeric(testdata[i, ])
        level <- mean(observed)
        pattern <- observed - level
        person_fit <- .pams_fit_no_intercept(pattern, profileOri)
        partial_correlations <- .pams_partial_correlations(pattern, profileOri)
        result[i, ] <- c(
            person_fit$coefficients,
            level,
            person_fit$r_squared,
            partial_correlations
        )
    }
    meanR2 <- if (all(is.na(result[, nprofile + 2L]))) {
        NA_real_
    } else {
        mean(result[, nprofile + 2L], na.rm = TRUE)
    }

    resultB <- resultBP <- NULL
    if (!is.null(participant)) {
        xBoot <- lapply(seq_len(nBoot), function(b) {
            do.call(cbind, lapply(profileBoot, function(x) x[b, ]))
        })
        weight_names <- unlist(lapply(
            seq_len(nprofile),
            function(j) c(paste0("w", j), paste0(c("m", "se", "L", "U"), j))
        ))
        pcorr_names <- unlist(lapply(
            seq_len(nprofile),
            function(j) c(paste0("corDim", j), paste0(c("m", "se", "L", "U"), j))
        ))
        resultB <- matrix(
            NA_real_,
            nrow = length(participant),
            ncol = 5L * nprofile + 2L + nprofile,
            dimnames = list(
                paste0("#", participant),
                c(weight_names, "level", "R^2", paste0("corDim", seq_len(nprofile)))
            )
        )
        resultBP <- matrix(
            NA_real_,
            nrow = length(participant),
            ncol = 5L * nprofile,
            dimnames = list(paste0("#", participant), pcorr_names)
        )

        for (i in seq_along(participant)) {
            observed <- as.numeric(testdata[participant[i], ])
            level <- mean(observed)
            pattern <- observed - level
            boot_weights <- matrix(NA_real_, nrow = nBoot, ncol = nprofile)
            boot_pcorr <- matrix(NA_real_, nrow = nBoot, ncol = nprofile)
            for (b in seq_len(nBoot)) {
                boot_weights[b, ] <- .pams_fit_no_intercept(pattern, xBoot[[b]])$coefficients
                boot_pcorr[b, ] <- .pams_partial_correlations(pattern, xBoot[[b]])
            }

            weight_summary <- unlist(lapply(seq_len(nprofile), function(j) {
                c(
                    result[participant[i], j],
                    mean(boot_weights[, j]),
                    stats::sd(boot_weights[, j]),
                    as.numeric(stats::quantile(boot_weights[, j], lalpha, names = FALSE)),
                    as.numeric(stats::quantile(boot_weights[, j], ualpha, names = FALSE))
                )
            }))
            pcorr_summary <- unlist(lapply(seq_len(nprofile), function(j) {
                c(
                    result[participant[i], nprofile + 2L + j],
                    mean(boot_pcorr[, j], na.rm = TRUE),
                    stats::sd(boot_pcorr[, j], na.rm = TRUE),
                    as.numeric(stats::quantile(
                        boot_pcorr[, j], lalpha, na.rm = TRUE, names = FALSE
                    )),
                    as.numeric(stats::quantile(
                        boot_pcorr[, j], ualpha, na.rm = TRUE, names = FALSE
                    ))
                )
            }))
            resultB[i, ] <- c(
                weight_summary,
                result[participant[i], nprofile + 1L],
                result[participant[i], nprofile + 2L],
                result[participant[i], (nprofile + 3L):(2L * nprofile + 2L)]
            )
            resultBP[i, ] <- pcorr_summary
        }
    }

    if (!is.null(file)) {
        profile_table <- do.call(cbind, profile)
        profile_names <- unlist(lapply(
            seq_len(nprofile),
            function(j) paste0(names(profile[[j]]), j)
        ))
        colnames(profile_table) <- profile_names
        mds_file <- paste0(file, "MDS.csv")
        cat("Summary Statistics for Stress\n", file = mds_file)
        if (is.null(stresssummary)) {
            cat("Not available for classical MDS\n\n", file = mds_file, append = TRUE)
        } else {
            utils::write.table(
                stresssummary,
                file = mds_file,
                sep = ",",
                row.names = FALSE,
                col.names = TRUE,
                append = TRUE
            )
            cat("\n", file = mds_file, append = TRUE)
        }
        cat("Summary Statistics for Profile\n", file = mds_file, append = TRUE)
        utils::write.table(
            profile_table,
            file = mds_file,
            sep = ",",
            row.names = TRUE,
            col.names = NA,
            append = TRUE
        )

        utils::write.table(
            result,
            file = paste0(file, "Weight.csv"),
            sep = ",",
            row.names = TRUE,
            col.names = NA
        )
        if (!is.null(resultB)) {
            utils::write.table(
                resultB,
                file = paste0(file, "WeightB.csv"),
                sep = ",",
                row.names = TRUE,
                col.names = NA
            )
            utils::write.table(
                resultBP,
                file = paste0(file, "PcorrB.csv"),
                sep = ",",
                row.names = TRUE,
                col.names = NA
            )
        }
    }

    output <- list(
        MDS = MDS,
        MDSsummary = profile,
        MDSprofile = profileBoot,
        stresssummary = stresssummary,
        stressprofile = stressBoot,
        MDSR2 = R2,
        Weight = result,
        WeightmeanR2 = meanR2,
        WeightB = resultB,
        PcorrB = resultBP,
        nprofile = nprofile,
        nBoot = nBoot,
        scale = scale,
        testname = testname,
        call = call,
        mds = mds,
        type = if (identical(mds, "smacof")) type else NA_character_,
        distance = distance,
        nsubject = nsubject,
        ntest = ntest,
        cl = cl,
        direction = direction
    )
    class(output) <- c("pams_fit", "list")
    output
}
