# run example from smoothness_metric() function documentation
example(smoothness_metric, echo = FALSE)


test_that("results object has correct structure", {
  expect_is(res, "list")
  expect_length(res, 2)
})

test_that("metric output values for example dataset are correct", {
  expect_equal(sum(res$n_discordant), 2046)
})
