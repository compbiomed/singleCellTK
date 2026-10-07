library(singleCellTK)

test_that(desc = "Testing runClusterSummaryMetrics.R", {
  data("scExample")
  
  B2M <- runClusterSummaryMetrics(sce, useAssay="counts", featureNames=c("B2M"), displayName="feature_name", groupNames="type")
  percExpr <- c(1, 0, 1)
  aveExpr <- c(94, 2, 54)
  
  testthat::expect_true(identical(floor(B2M$percExpr[1:3]), percExpr) & 
                          
                          identical(floor(B2M$avgExpr[1:3]), aveExpr))
  
  testthat::expect_error(runClusterSummaryMetrics(sce, useAssay="counts", featureNames=c("B2M"), displayName="feature_name", groupNames="howdy"),
                         "Specified variable 'howdy' not found in colData(inSCE)", fixed=TRUE)
  
  testthat::expect_warning(runClusterSummaryMetrics(sce, useAssay="counts", featureNames=c("B2M", "applesauce"), displayName="feature_name", groupNames="type"),
                         "Specified genes 'applesauce' not found in rowData(inSCE)$feature_name", fixed=TRUE)
  
  testthat::expect_error(runClusterSummaryMetrics(sce, useAssay="counts", featureNames=c("applesauce", "cowboy"), displayName="feature_name", groupNames="type"),
                         "All genes in 'applesauce, cowboy' not found in rowData(inSCE)$feature_name", fixed=TRUE)
})
test_that(desc = "runClusterSummaryMetrics matches base R and does not warn", {
  data("scExample")
  genes <- c("B2M", "MALAT1", "ACTB")
  expect_no_warning(
    res <- runClusterSummaryMetrics(sce, useAssay = "counts",
                                    featureNames = genes,
                                    displayName = "feature_name",
                                    groupNames = "type")
  )
  ix <- match(genes, rowData(sce)$feature_name)
  mat <- as.matrix(counts(sce)[ix, ])
  groups <- sort(unique(sce$type))
  expMean <- vapply(groups, function(g) rowMeans(mat[, sce$type == g]),
                    numeric(length(genes)))
  expDet <- vapply(groups, function(g) rowMeans(mat[, sce$type == g] > 0),
                   numeric(length(genes)))
  expect_equal(unname(res$avgExpr), unname(expMean))
  expect_equal(unname(res$percExpr), unname(expDet))
  expect_equal(colnames(res$avgExpr), groups)
  expect_equal(colnames(res$percExpr), groups)
  expect_equal(res$featureNames, genes)
})
