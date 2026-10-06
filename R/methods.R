#' Summarize a fitted PAMS model
#'
#' @param object A fitted object returned by [BootSmacof()].
#' @param digits Number of digits retained when the summary is printed.
#' @param ... Additional arguments passed to or from methods.
#'
#' @return `summary.pams_fit()` returns an object of class
#'   `summary.pams_fit` containing the analysis settings, stress, fit indices,
#'   and the number of coordinates whose pointwise BCa interval excludes zero.
#'
#' @examples
#' \donttest{
#' set.seed(42)
#' x <- matrix(rnorm(32), nrow = 8, ncol = 4)
#' fit <- suppressWarnings(BootSmacof(x, nprofile = 2, nBoot = 10))
#' summary(fit)
#' }
#'
#' @export
summary.pams_fit <- function(object, digits = 3L, ...) {
    digits <- .pams_integer(digits, "digits", minimum = 1L)
    significant <- vapply(object$MDSsummary, function(x) {
        sum(x$BCaLower > 0 | x$BCaUpper < 0, na.rm = TRUE)
    }, integer(1L))
    result <- list(
        call = object$call,
        mds = object$mds,
        type = object$type,
        distance = object$distance,
        nsubject = object$nsubject,
        ntest = object$ntest,
        nprofile = object$nprofile,
        nBoot = object$nBoot,
        cl = object$cl,
        stress = if (is.null(object$MDS$stress)) NA_real_ else object$MDS$stress,
        stresssummary = object$stresssummary,
        WeightmeanR2 = object$WeightmeanR2,
        MDSR2 = object$MDSR2,
        significant_coordinates = significant,
        digits = digits
    )
    class(result) <- "summary.pams_fit"
    result
}

#' @rdname summary.pams_fit
#' @param x An object returned by `summary.pams_fit()` or [BootSmacof()].
#' @export
print.summary.pams_fit <- function(x, ...) {
    digits <- x$digits
    cat("Profile Analysis via Multidimensional Scaling\n\n")
    cat("Call:\n")
    print(x$call)
    cat("\nModel:", x$mds)
    if (!is.na(x$type)) cat("(", x$type, ")")
    cat("with", x$distance, "distances\n")
    cat("Persons:", x$nsubject, " | Variables:", x$ntest,
        " | Core profiles:", x$nprofile, "\n")
    cat("Bootstrap samples:", x$nBoot, " | Confidence level:", x$cl, "\n")
    if (!is.na(x$stress)) cat("Original-sample stress:", round(x$stress, digits), "\n")
    cat("Mean person-level R-squared:", round(x$WeightmeanR2, digits), "\n")
    cat("Core-profile collinearity R-squared:\n")
    print(round(x$MDSR2, digits))
    cat("Coordinates with pointwise BCa intervals excluding zero:\n")
    print(x$significant_coordinates)
    invisible(x)
}

#' @rdname summary.pams_fit
#' @export
print.pams_fit <- function(x, ...) {
    print(summary(x), ...)
    invisible(x)
}

#' Plot fitted PAMS core profiles
#'
#' Displays selected core-profile coordinates and, optionally, their
#' pointwise percentile or BCa confidence limits.
#'
#' @param x A fitted object returned by [BootSmacof()].
#' @param profiles Integer vector identifying the profiles to display.
#' @param interval Confidence limits to draw: `"BCa"` (default),
#'   `"percentile"`, or `"none"`.
#' @param layout Optional integer vector of length two giving the rows and
#'   columns of the plotting layout.
#' @param labels Logical; show variable names on the horizontal axis.
#' @param label.cex Character expansion used for variable labels.
#' @param type,pch,col Standard graphical controls for the estimated profile.
#' @param ci.col,ci.lty Color and line type for confidence limits.
#' @param zero.line Logical; draw a horizontal reference line at zero.
#' @param zero.col Color of the zero reference line.
#' @param xlab,ylab Axis labels.
#' @param main Optional plot title or character vector of titles.
#' @param ylim Optional common vertical-axis limits. By default each panel
#'   includes its estimate and selected confidence limits.
#' @param xaxt Axis annotation setting passed to [graphics::plot()].
#' @param ... Additional graphical arguments passed to [graphics::plot()].
#'
#' @return The fitted object, invisibly.
#'
#' @examples
#' \donttest{
#' set.seed(42)
#' x <- matrix(rnorm(32), nrow = 8, ncol = 4)
#' fit <- suppressWarnings(BootSmacof(x, nprofile = 2, nBoot = 10))
#' plot(fit, profiles = 1:2)
#' }
#'
#' @export
plot.pams_fit <- function(x,
                          profiles = seq_len(x$nprofile),
                          interval = c("BCa", "percentile", "none"),
                          layout = NULL,
                          labels = TRUE,
                          label.cex = 0.75,
                          type = "b",
                          pch = 19,
                          col = "black",
                          ci.col = "steelblue",
                          ci.lty = 2,
                          zero.line = TRUE,
                          zero.col = "grey60",
                          xlab = "Variable",
                          ylab = "Coordinate",
                          main = NULL,
                          ylim = NULL,
                          xaxt = if (labels) "n" else "s",
                          ...) {
    interval <- match.arg(interval)
    if (!is.numeric(profiles) || !length(profiles) || anyNA(profiles) ||
        any(profiles != floor(profiles)) ||
        any(profiles < 1L | profiles > x$nprofile) || anyDuplicated(profiles)) {
        stop("`profiles` must contain unique profile indices between 1 and `x$nprofile`.",
             call. = FALSE)
    }
    profiles <- as.integer(profiles)
    if (is.null(layout)) {
        layout <- if (length(profiles) <= 3L) {
            c(1L, length(profiles))
        } else {
            rows <- ceiling(sqrt(length(profiles)))
            c(rows, ceiling(length(profiles) / rows))
        }
    }
    if (!is.numeric(layout) || length(layout) != 2L ||
        anyNA(layout) || any(layout < 1) || any(layout != floor(layout))) {
        stop("`layout` must be NULL or two positive integers.", call. = FALSE)
    }
    if (!is.logical(labels) || length(labels) != 1L || is.na(labels) ||
        !is.logical(zero.line) || length(zero.line) != 1L || is.na(zero.line)) {
        stop("`labels` and `zero.line` must each be TRUE or FALSE.", call. = FALSE)
    }

    oldpar <- graphics::par(no.readonly = TRUE)
    on.exit(graphics::par(oldpar), add = TRUE)
    graphics::par(mfrow = as.integer(layout))

    profile_col <- rep_len(col, length(profiles))
    confidence_col <- rep_len(ci.col, length(profiles))
    titles <- if (is.null(main)) {
        paste("Core profile", profiles)
    } else {
        rep_len(main, length(profiles))
    }

    for (index in seq_along(profiles)) {
        profile_number <- profiles[index]
        values <- x$MDSsummary[[profile_number]]
        lower_name <- if (identical(interval, "BCa")) "BCaLower" else "Lower"
        upper_name <- if (identical(interval, "BCa")) "BCaUpper" else "Upper"
        limits <- if (identical(interval, "none")) {
            values$Ori
        } else {
            c(values$Ori, values[[lower_name]], values[[upper_name]])
        }
        panel_ylim <- if (is.null(ylim)) range(limits, finite = TRUE) else ylim
        horizontal <- seq_len(nrow(values))
        graphics::plot(
            horizontal,
            values$Ori,
            type = type,
            pch = pch,
            col = profile_col[index],
            xlab = xlab,
            ylab = ylab,
            main = titles[index],
            ylim = panel_ylim,
            xaxt = xaxt,
            ...
        )
        if (labels) {
            graphics::axis(
                1,
                at = horizontal,
                labels = x$testname,
                las = 2,
                cex.axis = label.cex
            )
        }
        if (!identical(interval, "none")) {
            graphics::lines(
                horizontal,
                values[[lower_name]],
                col = confidence_col[index],
                lty = ci.lty
            )
            graphics::lines(
                horizontal,
                values[[upper_name]],
                col = confidence_col[index],
                lty = ci.lty
            )
        }
        if (zero.line) graphics::abline(h = 0, col = zero.col, lty = 3)
    }
    invisible(x)
}
