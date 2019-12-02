context("MSDial functions")
library(testthat)

test_that("get_samples_info get accurate number of samples", {
    test_dir <- system.file("extdata", "dirtest", package = "mscleanr")
    expect_equal(nrow(get_samples_info(test_dir)), 8)
    expect_equal(nrow(get_samples_info(test_dir, source="pos")), 8)
    expect_equal(nrow(get_samples_info(test_dir, source="neg")), 8)
})
