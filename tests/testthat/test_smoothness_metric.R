# run example from smoothness_metric() function documentation
example(smoothness_metric, echo = FALSE)


test_that("results object has correct structure", {
  expect_is(res, "list")
  expect_length(res, 2)
})

test_that("value of metric for example dataset is correct", {
  expect_equal(signif(res$mean_discordant, 4), 0.5622)
})
