# Enrichment Analysis
library(singleCellTK)
context("Testing enrichment analysis functions")
data(scExample, package = "singleCellTK")
sce <- subsetSCECols(sce, colData = "type != 'EmptyDroplet'")
internetConnection <- suppressWarnings(Biobase::testBioCConnection())

test_that(desc = "Testing wrong getter method usage", {
    testthat::expect_error(getEnrichRResult(sce, "analysis1"),
                           "EnrichR analysis not performed yet")
    getEnrichRResult(sce, "analysis1") <- "placeholder"
    testthat::expect_error(getEnrichRResult(sce, "analysis2"), 
                           '"analysis2" not found in EnrichR analysis names')
})

test_that(desc = "Testing wrong function usage", {
    testthat::expect_error(runEnrichR("character"), 
                           "inSCE has to inherit")
    testthat::expect_error(runEnrichR(sce, "hello", "analysis2"), 
                           "Not all features found in `rownames[(]inSCE[)]`.")
    testthat::expect_error(runEnrichR(sce, rownames(sce)[1:5], "analysis2", 
                                      by="ensembl"), 
                           "`by` not found in rowData[(]inSCE[)].")
    testthat::expect_error(runEnrichR(sce, rownames(sce)[1:5], "analysis2", 
                                      by="feature_name"), 
                           "in `rowData[(]inSCE[)][$]feature_name`.")
    testthat::expect_error(runEnrichR(sce, rownames(sce)[1:5], "analysis2", 
                                      featureName = "sample"), 
                           "featureName not found in")
    testthat::expect_error(runEnrichR(sce, rownames(sce)[1:5], "analysis2", 
                                      featureName = c("gene1", "gene2")), 
                           "Invalid featureName specification.")
    testthat::expect_error(runEnrichR(sce, rownames(sce)[1:5], "analysis2", 
                                      featureName = c("gene1", "gene2")), 
                           "Invalid featureName specification.")
    
    testthat::skip_if_not(internetConnection)
    testthat::expect_error(runEnrichR(sce, rownames(sce)[1:5], "analysis2", 
                                      featureName = "feature_name",
                                      db = "db1"), 
                           "database db1 do not exist.")
})

test_that(desc = "Testing correct function usage", {
    testthat::skip_if_not(internetConnection)
    sce <- runEnrichR(sce, rownames(sce)[1:5], "analysis2", 
                      featureName = "feature_name",
                      db = "HDSigDB_Human_2021")
    testthat::expect_true("analysis2" %in% 
                              names(metadata(sce)$sctk$runEnrichR))
    res <- getEnrichRResult(sce, "analysis2")
    testthat::expect_equal(ncol(res$result), 11)
    testthat::expect_gt(nrow(res$result), 200)
})
