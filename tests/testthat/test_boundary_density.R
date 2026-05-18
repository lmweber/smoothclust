# run example from boundary_density() function documentation
example(boundary_density, echo = FALSE)


test_that("results object has correct structure", {
  expect_is(res, "list")
  expect_length(res, 2)
})
