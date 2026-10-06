.pams_integer <- function(x, name, minimum = 1L) {
    if (length(x) != 1L || is.na(x) || !is.numeric(x) || !is.finite(x) ||
        x != floor(x) || x < minimum) {
        stop(sprintf("`%s` must be a single integer greater than or equal to %d.",
                     name, minimum),
             call. = FALSE)
    }
    as.integer(x)
}

.pams_distance <- function(x, distance) {
    out <- stats::dist(t(x))
    if (identical(distance, "sqeuclid")) out <- out^2
    out
}

.align_pams_signs <- function(configuration, reference) {
    correlations <- diag(stats::cor(configuration, reference))
    signs <- ifelse(is.na(correlations) | correlations == 0, 1, sign(correlations))
    sweep(configuration, 2L, signs, `*`)
}

.pams_fit_no_intercept <- function(y, x) {
    fit <- stats::lm.fit(x = x, y = y)
    coefficients <- stats::coef(fit)
    coefficients[is.na(coefficients)] <- 0
    residuals <- as.numeric(y - x %*% coefficients)
    denominator <- sum(y^2)
    r_squared <- if (denominator > 0) 1 - sum(residuals^2) / denominator else NA_real_
    list(coefficients = coefficients, residuals = residuals, r_squared = r_squared)
}

.pams_residualize <- function(y, x) {
    if (is.null(x) || ncol(x) == 0L) return(as.numeric(y))
    .pams_fit_no_intercept(y, x)$residuals
}

.pams_safe_cor <- function(x, y) {
    if (length(x) < 2L || stats::sd(x) == 0 || stats::sd(y) == 0) return(NA_real_)
    stats::cor(x, y)
}

.pams_partial_correlations <- function(y, configuration) {
    nprofile <- ncol(configuration)
    result <- numeric(nprofile)
    for (j in seq_len(nprofile)) {
        others <- setdiff(seq_len(nprofile), j)
        other_configuration <- if (length(others)) {
            configuration[, others, drop = FALSE]
        } else {
            NULL
        }
        y_residual <- .pams_residualize(y, other_configuration)
        profile_residual <- .pams_residualize(configuration[, j], other_configuration)
        result[j] <- .pams_safe_cor(y_residual, profile_residual)
    }
    result
}

.pams_profile_r2 <- function(configuration) {
    nprofile <- ncol(configuration)
    result <- numeric(nprofile)
    if (nprofile == 1L) return(result)
    for (j in seq_len(nprofile)) {
        others <- setdiff(seq_len(nprofile), j)
        result[j] <- .pams_fit_no_intercept(
            configuration[, j],
            configuration[, others, drop = FALSE]
        )$r_squared
    }
    result
}

.pams_acceleration <- function(jackknife) {
    influence <- -(jackknife - mean(jackknife))
    denominator <- 6 * sum(influence^2)^1.5
    if (!is.finite(denominator) || denominator <= .Machine$double.eps) return(0)
    sum(influence^3) / denominator
}

.pams_bca_interval <- function(bootstrap, original, jackknife, probabilities) {
    n_boot <- length(bootstrap)
    proportion <- mean(bootstrap < original)
    boundary <- 1 / (2 * n_boot)
    proportion <- min(max(proportion, boundary), 1 - boundary)
    z0 <- stats::qnorm(proportion)
    acceleration <- .pams_acceleration(jackknife)
    z_alpha <- stats::qnorm(probabilities)
    adjusted <- stats::pnorm(
        z0 + (z0 + z_alpha) / (1 - acceleration * (z0 + z_alpha))
    )
    adjusted <- pmin(pmax(adjusted, 0), 1)
    as.numeric(stats::quantile(bootstrap, probs = adjusted, names = FALSE))
}

.pams_summary_row <- function(original, bootstrap, jackknife, lower, upper) {
    bca <- .pams_bca_interval(
        bootstrap = bootstrap,
        original = original,
        jackknife = jackknife,
        probabilities = c(lower, upper)
    )
    c(
        Ori = original,
        Mean = mean(bootstrap),
        SE = stats::sd(bootstrap),
        Lower = as.numeric(stats::quantile(bootstrap, lower, names = FALSE)),
        Upper = as.numeric(stats::quantile(bootstrap, upper, names = FALSE)),
        BCaLower = bca[1L],
        BCaUpper = bca[2L]
    )
}
