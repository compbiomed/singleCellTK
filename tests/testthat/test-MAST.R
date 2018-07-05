context("MAST testing")

data(maits, package="MAST")
maits_sce <- createSCE(assayFile = t(maits$expressionmat),
                       annotFile = maits$cdat,
                       featureFile = maits$fdat,
                       assayName = "logtpm",
                       inputDataFrames = TRUE,
                       createLogCounts = FALSE)

library(scRNAseq)

data(allen, package = "scRNAseq")

tempsce <- as(allen, "SingleCellExperiment")

original <- as(tempsce, "SCtkExperiment")


assay(original, "counts") <- assay(original, "tophat_counts")
assay(original, "logcounts") <- log2(assay(original, "counts")+1)
assay(original, "tophat_counts") <- NULL

test_that("MAST should fail if too few samples pass filter", {
  expect_error(MAST(maits_sce, condition = "condition", freqExpressed = "1", useThresh = TRUE, useAssay = "logtpm"),
               "Not enough genes pass frequency expressed filter of 1")
})
test_that("MAST should fail if data has NAs", {
  expect_error(MAST(original, condition = "Primary.Type", interest.level = "L4 Arf5", freqExpressed = "0.2", useThresh = TRUE, useAssay = "logcounts" ),
               "Data has NAs, use Filter samples by annotation to fix")
})

