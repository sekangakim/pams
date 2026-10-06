test_that("invalid data and scalar arguments fail informatively", {
    expect_error(BootSmacof("not data"), "numeric matrix or data frame")

    bad_type <- as.data.frame(pams_test_data)
    bad_type[[1]] <- letters[seq_len(nrow(bad_type))]
    expect_error(BootSmacof(bad_type), "Every column")

    missing_data <- pams_test_data
    missing_data[1, 1] <- NA_real_
    expect_error(BootSmacof(missing_data), "missing or infinite")

    expect_error(
        BootSmacof(pams_test_data, nprofile = ncol(pams_test_data)),
        "smaller than"
    )
    expect_error(
        BootSmacof(pams_test_data, nprofile = 2, direction = c(1, 0)),
        "1 or -1"
    )
    expect_error(BootSmacof(pams_test_data, nBoot = 9), "greater than or equal")
    expect_error(BootSmacof(pams_test_data, cl = 1), "strictly between")
    expect_error(
        BootSmacof(pams_test_data, participant = c(1, 1)),
        "unique"
    )
    expect_error(BootSmacof(pams_test_data, mds = "unknown"), "arg")
})
