test_that("summary produces concise fitted-model information", {
    out <- summary(pams_test_fit)
    expect_s3_class(out, "summary.pams_fit")
    expect_equal(out$nsubject, 12)
    expect_equal(out$nprofile, 2)
    expect_length(out$significant_coordinates, 2)
    expect_output(print(out), "Profile Analysis via Multidimensional Scaling")
    expect_output(print(pams_test_fit), "Mean person-level R-squared")
})

test_that("plot supports selected profiles and confidence choices", {
    path <- tempfile(fileext = ".pdf")
    grDevices::pdf(path)
    on.exit(grDevices::dev.off(), add = TRUE)
    expect_invisible(plot(pams_test_fit, profiles = 1, interval = "BCa"))
    expect_invisible(plot(pams_test_fit, profiles = 2, interval = "none"))
    expect_error(plot(pams_test_fit, profiles = 3), "between 1")
})

test_that("Weight can be converted directly to a tibble", {
    skip_if_not_installed("tibble")
    out <- tibble::as_tibble(pams_test_fit$Weight, rownames = "participant")
    expect_s3_class(out, "tbl_df")
    expect_equal(nrow(out), 12)
    expect_true(all(c("participant", "w1", "corDim1") %in% names(out)))
})
