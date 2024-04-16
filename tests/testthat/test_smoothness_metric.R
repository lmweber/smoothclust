# run example from smoothness_metric() function documentation
example(smoothness_metric, echo = FALSE)


test_that("results object has correct structure", {
  expect_is(res, "list")
  expect_length(res, 2)
})
