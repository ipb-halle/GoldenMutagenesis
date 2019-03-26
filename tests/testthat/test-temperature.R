context("Temperature")

test_that("Temperature calculation works correct (nnb)", {
  expect_identical(round(calculate_tm_nnb("GTATGTGTGTATATATATGT", salt_concentration = 50, primer_concentration = 200, offset = 0), digits = 5), 40.50067)
  expect_identical(round(calculate_tm_nnb("GTATGTGTGTATATATATGT", salt_concentration = 50, primer_concentration = 50, offset = 0), digits = 5), 38.26582)
})

test_that("Temperature calculation works correct", {
  expect_identical(round(calculate_tm("GTATGTGTGTATATATATGT", salt_concentration = 50, primer_concentration = 200, offset = 0), digits = 4), 48.1529)
  expect_identical(round(calculate_tm("GTATGTGTGTATATATATGT", salt_concentration = 50, primer_concentration = 50, offset = 0), digits = 4), 48.1529)
})