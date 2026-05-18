test_that("smoothclust validates parameter values", {
  input <- matrix(1:12, nrow = 3)
  spatial_coords <- cbind(0:3, 0)
  
  expect_error(smoothclust(input, spatial_coords = spatial_coords, bandwidth = 0))
  expect_error(smoothclust(input, spatial_coords = spatial_coords, bandwidth = c(0.1, 0.2)))
  expect_error(smoothclust(input, spatial_coords = spatial_coords, k = 1.5))
  expect_error(smoothclust(input, spatial_coords = spatial_coords, truncate = 1))
  expect_error(smoothclust(input, spatial_coords = spatial_coords, truncate = NA_real_))
  expect_error(smoothclust(input, spatial_coords = spatial_coords, n_threads = 0))
  expect_error(smoothclust(input, spatial_coords = spatial_coords, n_threads = 1.5))
})

test_that("smoothclust validates spatial coordinate compatibility", {
  input <- matrix(1:12, nrow = 3)
  spatial_coords <- cbind(0:3, 0)
  
  expect_error(smoothclust(input, spatial_coords = spatial_coords[-1, ]))
  
  spatial_coords_bad <- spatial_coords
  spatial_coords_bad[1, 1] <- Inf
  expect_error(smoothclust(input, spatial_coords = spatial_coords_bad))
})

test_that("smoothclust validates k when using knn smoothing", {
  input <- matrix(1:12, nrow = 3)
  spatial_coords <- cbind(0:3, 0)
  
  expect_error(smoothclust(input, spatial_coords = spatial_coords, method = "knn", k = 4))
})

test_that("boundary_density validates inputs", {
  spatial_coords <- cbind(0:3, 0)
  labels <- c(1, 1, 2, 2)
  
  expect_error(boundary_density(spatial_coords[-1, ], labels, k = 1))
  
  spatial_coords_bad <- spatial_coords
  spatial_coords_bad[1, 1] <- NA_real_
  expect_error(boundary_density(spatial_coords_bad, labels, k = 1))
  
  expect_error(boundary_density(spatial_coords, labels, k = 0))
  expect_error(boundary_density(spatial_coords, labels, k = 1.5))
  expect_error(boundary_density(spatial_coords, labels, k = 4))
  expect_error(boundary_density(spatial_coords, labels, n_threads = 0))
  expect_error(boundary_density(spatial_coords, labels, n_threads = Inf))
})
