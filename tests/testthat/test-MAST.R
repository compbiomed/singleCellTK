context("MAST testing")

test_that("MAST should fail if too few samples pass filter", {
  expect_error(MAST(mouseBrainSubsetSCE, condition = "level1class",
                    freqExpressed = "1", useThresh = TRUE,
                    useAssay = "logcounts"),
               "Not enough genes pass frequency expressed filter of 1")
})

test_that("MAST should fail if data has NAs", {
  tmpdata <- mouseBrainSubsetSCE
  colData(tmpdata)$testcondition <- c(rep(NA, ncol(tmpdata) - 1), "a")
  expect_error(MAST(tmpdata, condition = "testcondition",
                    freqExpressed = "0.2", useThresh = TRUE,
                    useAssay = "logcounts"),
               "Annotation data has NA values. Filter them to continue.")
})
