library(singleCellTK)

test_that(desc = "Testing runClusterSummaryMetrics.R", {
  data("scExample")
  
  B2M <- runClusterSummaryMetrics(sce, useAssay="counts", feature=c("B2M"), displayName="feature_name", clusters="type")
  percExpr <- c(1, 0, 1)
  aveExpr <- c(94, 2, 54)
  
  testthat::expect_true(identical(floor(B2M$percExpr), percExpr) & 
                          
                          identical(floor(B2M$avgExpr), aveExpr))
  
  testthat::expect_error(runClusterSummaryMetrics(sce, useAssay="counts", feature=c("B2M"), displayName="feature_name", clusters="howdy"),
                         "Specified variable 'howdy' not found in colData(inSCE)", fixed=TRUE)
  
  testthat::expect_warning(runClusterSummaryMetrics(sce, useAssay="counts", feature=c("B2M", "applesauce"), displayName="feature_name", clusters="type"),
                         "Specified genes 'applesauce' not found in rowData(inSCE)$feature_name", fixed=TRUE)
  
  testthat::expect_error(runClusterSummaryMetrics(sce, useAssay="counts", feature=c("applesauce", "cowboy"), displayName="feature_name", clusters="type"),
                         "All genes in 'applesauce, cowboy' not found in rowData(inSCE)$feature_name", fixed=TRUE)
})