

testthat::test_that("errors", {
  testthat::expect_error(
    StructuralDecompose(Data = c(1,2,3), frequency = '3', break_algorithm = 'strucchange', break_level = 0.05, median_level = 0.5, mean_level = 0.1, level_length = 0.5, conf_level = 0.5),
    "Value needs to be numeric"
  )

  testthat::expect_error(
    StructuralDecompose(Data = c(1,2,3), frequency = '2', break_algorithm = 'strucchange', break_level = 0.05, median_level = 0.5, mean_level = 0.1, level_length = 0.5, conf_level = 0.5),
    "Value needs to be numeric"
  )
})

# Data,
# frequency = 12,
# break_algorithm = 'strucchange',
# smoothening_algorithm = 'lowess',
# break_level = 0.05,
# median_level = 0.5,
# mean_level = 0.5,
# level_length = 0.5,
# conf_level = 0.5
