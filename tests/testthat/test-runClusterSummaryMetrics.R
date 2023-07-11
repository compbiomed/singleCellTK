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