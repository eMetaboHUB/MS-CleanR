context("MSDial functions")
library(testthat)

test_that("get_samples_info get accurate number of samples", {
    choose_project_directory(system.file("extdata", "dirtest", package = "mscleanr"))
    expect_equal(nrow(get_samples_info()), 8)
    expect_equal(nrow(get_samples_info(source="pos")), 8)
    expect_equal(nrow(get_samples_info(source="neg")), 8)
})
