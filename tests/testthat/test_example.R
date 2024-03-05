# run example from smoothclust() function documentation
example(smoothclust, echo = FALSE)


test_that("exapmle object has correct class", {
  expect_s4_class(spe, "SpatialExperiment")
})

test_that("example object has correct dimensions", {
  expect_equal(dim(spe), c(33538, 3639))
})

test_that("example object has correct assays", {
  expect_equal(assayNames(spe), c("counts", "counts_unsmoothed"))
})

test_that("example object has correct first few values", {
  expect_equal(unname(signif(counts(spe)[1:6, 1], 7)), 
               c(0, 0, 0, 0.1428571, 0, 0))
})
