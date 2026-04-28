test_that("smoothclust includes self in uniform neighborhoods", {
  vals <- matrix(1:4, nrow = 1)
  spatial_coords <- cbind(c(0, 1, 3, 6), 0)
  
  out <- smoothclust(vals, spatial_coords = spatial_coords, 
                     method = "uniform", bandwidth = 1.1 / 6)
  
  expect_equal(as.numeric(out), c(1.5, 1.5, 3, 4))
})

test_that("smoothclust includes self in kernel neighborhoods", {
  vals <- matrix(1:3, nrow = 1)
  spatial_coords <- cbind(c(0, 100, 101), 0)
  
  out <- smoothclust(vals, spatial_coords = spatial_coords, 
                     method = "kernel", bandwidth = 1 / 101, truncate = 0.14)
  
  w <- exp(-1)
  expect_equal(as.numeric(out), c(1, (2 + 3 * w) / (1 + w), 
                                  (3 + 2 * w) / (1 + w)))
})

test_that("smoothclust includes self in knn neighborhoods", {
  vals <- matrix(1:4, nrow = 1)
  spatial_coords <- cbind(c(0, 1, 3, 6), 0)
  
  out <- smoothclust(vals, spatial_coords = spatial_coords, method = "knn", k = 1)
  
  expect_equal(as.numeric(out), c(1.5, 1.5, 2.5, 3.5))
})

test_that("smoothness_metric does not skip the nearest non-self neighbor", {
  spatial_coords <- cbind(c(0, 1, 3, 6), 0)
  labels <- c(1, 2, 1, 1)
  
  res <- smoothness_metric(spatial_coords, labels, k = 1)
  
  expect_equal(res$n_discordant, c(1, 1, 1, 0))
  expect_equal(res$mean_discordant, 0.75)
})
