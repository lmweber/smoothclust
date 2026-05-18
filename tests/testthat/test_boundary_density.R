test_that("results object has correct structure", {
  spatial_coords <- cbind(c(0, 1, 3, 6), 0)
  labels <- c("A", "B", "A", "A")
  
  res <- boundary_density(spatial_coords, labels, k = 2)
  
  expect_is(res, "list")
  expect_named(res, c("n_discordant", "n_neighbors", "local_boundary_density", 
                      "mean_discordant", "boundary_density"))
})

test_that("boundary_density returns analytic composition adjustment", {
  spatial_coords <- cbind(c(0, 1, 3, 6), 0)
  labels <- c(1, 2, 1, 1)
  
  res <- boundary_density(spatial_coords, labels, k = 2, adjust = "analytic")
  
  expect_equal(res$boundary_density, 0.625)
  expect_equal(res$expected_boundary_density, 0.5)
  expect_equal(res$relative_boundary_density, 1.25)
  expect_equal(res$excess_boundary_density, 0.125)
})

test_that("boundary_density returns permutation composition adjustment", {
  spatial_coords <- cbind(c(0, 1, 3, 6), 0)
  labels <- c(1, 2, 1, 1)
  
  res1 <- boundary_density(spatial_coords, labels, k = 2, 
                           adjust = "permutation", n_permutations = 20, 
                           seed = 123)
  res2 <- boundary_density(spatial_coords, labels, k = 2, 
                           adjust = "permutation", n_permutations = 20, 
                           seed = 123)
  
  expect_equal(res1$permutation_mean, res2$permutation_mean)
  expect_equal(res1$permutation_sd, res2$permutation_sd)
  expect_true(is.finite(res1$permutation_mean))
  expect_true(is.finite(res1$permutation_sd))
  expect_true(is.finite(res1$permutation_relative_boundary_density))
  expect_true(is.finite(res1$permutation_z))
  expect_equal(res1$n_permutations, 20)
})

test_that("boundary_density handles one-label composition references", {
  spatial_coords <- cbind(c(0, 1, 3, 6), 0)
  labels <- c(1, 1, 1, 1)
  
  analytic <- boundary_density(spatial_coords, labels, k = 2, 
                               adjust = "analytic")
  permutation <- boundary_density(spatial_coords, labels, k = 2, 
                                  adjust = "permutation", 
                                  n_permutations = 20, seed = 123)
  
  expect_equal(analytic$boundary_density, 0)
  expect_equal(analytic$expected_boundary_density, 0)
  expect_true(is.na(analytic$relative_boundary_density))
  expect_equal(analytic$excess_boundary_density, 0)
  expect_equal(permutation$permutation_mean, 0)
  expect_equal(permutation$permutation_sd, 0)
  expect_true(is.na(permutation$permutation_relative_boundary_density))
  expect_true(is.na(permutation$permutation_z))
})
