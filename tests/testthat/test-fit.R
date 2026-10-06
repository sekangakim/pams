test_that("BootSmacof returns a backward-compatible classed fit", {
    expect_s3_class(pams_test_fit, "pams_fit")
    expect_named(
        pams_test_fit,
        c(
            "MDS", "MDSsummary", "MDSprofile", "stresssummary",
            "stressprofile", "MDSR2", "Weight", "WeightmeanR2",
            "WeightB", "PcorrB", "nprofile", "nBoot", "scale",
            "testname", "call", "mds", "type", "distance", "nsubject",
            "ntest", "cl", "direction"
        ),
        ignore.order = FALSE
    )
    expect_equal(dim(pams_test_fit$Weight), c(12, 6))
    expect_equal(dim(pams_test_fit$MDSsummary[[1]]), c(5, 7))
    expect_equal(dim(pams_test_fit$WeightB), c(2, 14))
    expect_equal(dim(pams_test_fit$PcorrB), c(2, 10))
    expect_true(all(is.finite(pams_test_fit$Weight[, 1:4])))
})

test_that("the original-sample fit agrees with smacofSym", {
    direct <- smacof::smacofSym(
        stats::dist(t(pams_test_data)),
        ndim = 2,
        type = "ratio"
    )
    expect_equal(pams_test_fit$MDS$stress, direct$stress, tolerance = 1e-10)
    expect_equal(pams_test_fit$MDS$conf, direct$conf, tolerance = 1e-8)
})

test_that("direction reverses only the requested original-sample axis", {
    set.seed(1042)
    reversed <- suppressWarnings(BootSmacof(
        pams_test_data,
        mds = "smacof",
        type = "ratio",
        nprofile = 2,
        direction = c(-1, 1),
        nBoot = 10
    ))
    expect_equal(
        reversed$MDS$conf[, 1],
        -pams_test_fit$MDS$conf[, 1],
        tolerance = 1e-8
    )
    expect_equal(
        reversed$MDS$conf[, 2],
        pams_test_fit$MDS$conf[, 2],
        tolerance = 1e-8
    )
    expect_equal(reversed$MDS$stress, pams_test_fit$MDS$stress)
})

test_that("the classical-MDS branch has documented stress outputs", {
    set.seed(1042)
    classical <- suppressWarnings(BootSmacof(
        pams_test_data,
        mds = "classical",
        nprofile = 2,
        nBoot = 10
    ))
    expect_s3_class(classical, "pams_fit")
    expect_null(classical$MDS$stress)
    expect_null(classical$stresssummary)
    expect_null(classical$stressprofile)
    expect_equal(dim(classical$MDS$conf), c(5, 2))
})
