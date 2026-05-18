# run example from boundary_density() function documentation
example(boundary_density, echo = FALSE)


test_that("results object has correct structure", {
  expect_is(res, "list")
  expect_named(res, c("n_discordant", "n_neighbors", "local_boundary_density", 
                      "mean_discordant", "boundary_density"))
})
