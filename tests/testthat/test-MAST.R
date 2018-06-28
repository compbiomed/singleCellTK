context("MAST testing")

data(maits, package="MAST")
maits_sce <- createSCE(assayFile = t(maits$expressionmat),
                       annotFile = maits$cdat,
                       featureFile = maits$fdat,
                       assayName = "logtpm",
                       inputDataFrames = TRUE,
                       createLogCounts = FALSE)
rm(maits)

test_that("MAST should fail if too few samples pass filter", {
  expect_error(MAST(maits_sce, condition = "condition", freqExpressed = "1", useThresh = TRUE, useAssay = "logtpm"),
               "Not enough genes pass frequency expressed filter of 1")
})

