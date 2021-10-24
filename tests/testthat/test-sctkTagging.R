# sctkTagging functions
library(singleCellTK)
context("Testing sctkTagging functions")
data(scExample, package = "singleCellTK")

test_that(desc = "Testing tagging mechanism", {
  sce <- expSetDataTag(
    inSCE = sce, 
    assayType = "raw",
    assays = assayNames(sce))
  testthat::expect_true(!is.null(metadata(sce)$assayType))
  testthat::expect_equal(metadata(sce)$assayType$assayTag[1], "raw")
  testthat::expect_equal(metadata(sce)$assayType$assayName[1], "counts")
  
  tags <- expTaggedData(
    inSCE = sce,
    tags = "raw",
    redDims = TRUE,
    recommended = "raw",
    showTags = FALSE
  )
  testthat::expect_length(tags, 1)
  testthat::expect_equal(tags, "counts")
  testthat::expect_true(inherits(tags, "character"))
  
  tags <- expTaggedData(
    inSCE = sce,
    tags = "raw",
    redDims = TRUE,
    recommended = "raw",
    showTags = TRUE
  )
  testthat::expect_length(tags, 1)
  testthat::expect_equal(names(tags), "raw (recommended)")
  testthat::expect_true(inherits(tags, "list"))
  
  counts <- expData(sce, "counts")
  testthat::expect_true(inherits(expData(sce, "counts"), "dgCMatrix"))
  
  expData(sce, "counts") <- counts
  testthat::expect_true(inherits(expData(sce, "counts"), "dgCMatrix"))
  
  assayNames <- expDataNames(sce)
  testthat::expect_length(assayNames, 1)
  testthat::expect_equal(assayNames, "counts")
  
  sce <- expDeleteDataTag(sce, "counts")
  testthat::expect_true(!is.null(metadata(sce)$assayType))
  testthat::expect_equal(nrow(metadata(sce)$assayType), 0)
})
