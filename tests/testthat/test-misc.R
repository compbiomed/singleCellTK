context("misc functions")

test_that("summarizeSCE", {
  ta <- summarizeSCE(mouseBrainSubsetSCE, sample = NULL)
  expect_is(ta, "data.frame")
})
