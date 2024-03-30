# run example from smoothclust() function documentation
example(smoothclust, echo = FALSE)


test_that("example object has correct class", {
  expect_s4_class(spe, "SpatialExperiment")
})

test_that("example object has correct dimensions", {
  expect_equal(dim(spe), c(33538, 3639))
})

test_that("example object has correct assays", {
  expect_equal(assayNames(spe), c("counts", "counts_smooth"))
})

test_that("first few output values in example object are correct", {
  expect_equal(
    unname(signif(as.numeric(assays(spe)[["counts_smooth"]][1:6, 1], 6))), 
    c(0, 0, 0, 0.142857, 0, 0))
})
