# create random data
set.seed(123)
input <- matrix(rpois(2000, 10), ncol = 100)
spatial_coords <- matrix(runif(200), ncol = 2)

# run smoothclust with non-SpatialExperiment inputs
set.seed(123)
out <- smoothclust(input, spatial_coords = spatial_coords)


test_that("smoothclust runs with non-SPE inputs", {
  expect_s4_class(out, "dgTMatrix")
  expect_s4_class(out, "Matrix")
  expect_equal(dim(out), c(20, 100))
})
