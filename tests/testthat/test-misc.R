context("misc functions")

data(scExample, package = "singleCellTK")
sce <- subsetSCECols(sce, colData = "type != 'EmptyDroplet'")

test_that("summarizeSCE", {
  ta <- summarizeSCE(sce, sample = NULL)
  expect_is(ta, "data.frame")
})

test_that(desc = "Testing sampleSummaryStats", {
  data(scExample)
  stats <- sampleSummaryStats(sce, simple = FALSE)
  expect_true("matrix" %in% class(stats))
})


test_that(desc = "Testing subsetSCECols", {
  data(scExample)
  sce <- subsetSCECols(sce, colData = "type != 'EmptyDroplet'")
  expect_true(sum(colData(sce)$type == "EmptyDroplet") == 0)
})

test_that(desc = "Testing subsetSCERows", {
  data(scExample)
  rowData(sce)$isMito <- ifelse(grepl("^MT-", rowData(sce)$feature_name),"yes", "no")
  sce <- subsetSCERows(sce, rowData = "isMito == 'yes'")
  expect_true(length(altExps(sce)) == 1)
})

