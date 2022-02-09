context("misc functions")

### this code runs outside of test_that({})
# data(scExample, package = "singleCellTK")
# sce <- subsetSCECols(sce, colData = "type != 'EmptyDroplet'")

test_that("summarizeSCE", {
  data(scExample, package = "singleCellTK")
  sce <- subsetSCECols(sce, colData = "type != 'EmptyDroplet'")
  ta <- summarizeSCE(sce, sample = NULL)
  expect_is(ta, "data.frame")
})

test_that(desc = "Testing sampleSummaryStats", {
  data(scExample)
  sce <- sampleSummaryStats(sce, simple = FALSE)
  expect_true("sctk" %in% names(metadata(sce)))
  expect_is(metadata(sce)$sctk$sample_summary$qc_table,
            "matrix")
  expect_is(getSampleSummaryStatsTable(sce, statsName = "qc_table"),
            "matrix")
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


test_that(desc = "Testing runVAM", {
  data(scExample)
  sce <- subsetSCECols(sce, colData = "type != 'EmptyDroplet'")
  sce <- scaterlogNormCounts(sce, assayName = "logcounts")
  sce <- importGeneSetsFromMSigDB(inSCE = sce,
                                  categoryIDs = "H",
                                  species = "Homo sapiens",
                                  mapping = "gene_symbol",
                                  by = "feature_name")

  sce <- runVAM(inSCE = sce, geneSetCollectionName = "H", useAssay = "logcounts")
  expect_true(validObject(reducedDim(sce)))
})


test_that(desc = "Testing runGSVA", {
  data(scExample)
  sce <- subsetSCECols(sce, colData = "type != 'EmptyDroplet'")
  sce <- scaterlogNormCounts(sce, assayName = "logcounts")
  sce <- importGeneSetsFromMSigDB(inSCE = sce,
                                  categoryIDs = "H",
                                  species = "Homo sapiens",
                                  mapping = "gene_symbol",
                                  by = "feature_name")

  sce <- runGSVA(inSCE = sce, geneSetCollectionName = "H", useAssay = "logcounts")
  expect_true(validObject(reducedDim(sce)))
})

