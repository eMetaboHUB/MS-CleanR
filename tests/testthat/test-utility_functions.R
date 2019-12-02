context("Utility functions")
library(testthat)

test_that("extract_concatenated_data_from_column works fine", {
    test_data <- data.frame(test = c(
        "Classyfire_subclass=Sesquiterpenoids,Class=1b,Internal_id=1037912,Links=dnp:JHF25",
        "Classyfire_subclass=,Class=1b,Internal_id=255710,Links=unpd:UNPD163703",
        "Classyfire_subclass=,Class=1b,Internal_id=1056487,Links=dnp:LLL10",
        "Classyfire_subclass=Hydroxycinnamic acids and derivatives,Class=1b,Internal_id=979382,Links=dnp:QJJ32",
        "Classyfire_subclass=1-benzopyrans,Class=1b,Internal_id=255097,Links=dnp:BFM61|unpd:UNPD19697|cas:17623-62-0",
        "Classyfire_subclass=,Class=1b,Internal_id=967765,Links=dnp:MJT17")
        )
    test_other_name <- data.frame(other_name = test_data$test)
    test_data_na <- data.frame(test = c(
        "Classyfire_subclass=,Class=,Internal_id=,Links=",
        "Classyfire_subclass=,Class=,Internal_id=,Links=",
        "Classyfire_subclass=,Class=,Internal_id=,Links=",
        "Classyfire_subclass=,Class=,Internal_id=,Links=",
        "Classyfire_subclass=,Class=,Internal_id=,Links=",
        "Classyfire_subclass=,Class=,Internal_id=,Links=")
    )
    for (test_df in list(test_data, test_other_name, test_data_na)) {
        test_extract <- extract_concatenated_data_from_column(test_df)
        expect_equal(ncol(test_extract), 4)
        expect_equal(nrow(test_extract), 6)
        expect_true("Classyfire_subclass" %in% names(test_extract))
        expect_true("Class" %in% names(test_extract))
        expect_true("Internal_id" %in% names(test_extract))
        expect_true("Links" %in% names(test_extract))
        expect_false("test" %in% names(test_extract))
    }

    test_extract <- extract_concatenated_data_from_column(test_data)
    expect_equal(nrow(test_extract[is.na(test_extract$Classyfire_subclass),]), 3)

    test_extract <- extract_concatenated_data_from_column(test_data_na)
    expect_equal(nrow(test_extract[is.na(test_extract$Classyfire_subclass),]), 6)

    expect_equal(extract_concatenated_data_from_column(data.frame(test = c("truc"))), NULL)
    expect_equal(extract_concatenated_data_from_column(data.frame(test = c(""))), NULL)
})
