set.seed(1042)
pams_test_data <- matrix(rnorm(12 * 5), nrow = 12, ncol = 5)
colnames(pams_test_data) <- paste0("V", 1:5)

set.seed(1042)
pams_test_fit <- suppressWarnings(BootSmacof(
    pams_test_data,
    mds = "smacof",
    type = "ratio",
    nprofile = 2,
    nBoot = 10,
    participant = 1:2
))
