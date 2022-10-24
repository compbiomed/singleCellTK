library(singleCellTK)

test_that(desc = "Testing runClusterSummaryMetrics.R", {
  data("pbmc3k_2.7.1_sce")
  topMarkers <- getFindMarkerTopTable(pbmc3k_2.7.1_sce, topN = 1, log2fcThreshold = 0, fdrThreshold = 0.05,
                                      minClustExprPerc = 0.5, maxCtrlExprPerc = 0.5,
                                      minMeanExpr = 0)
  IL7R <- dplyr::select(dplyr::filter(topMarkers, Gene=="IL7R"), clusterExprPerc, clusterAveExpr)
  rownames(IL7R) <- NULL
  test <- runClusterSummaryMetrics(pbmc3k_2.7.1_sce, useAssay="logcounts", "IL7R", clusters="cluster")
  test <- dplyr::select(dplyr::filter(test, cluster==1), clusterExprPerc, clusterAveExpr)
  rownames(test) <- NULL
  test$clusterExprPerc <- as.numeric(test$clusterExprPerc)
  
  testthat::expect_true(identical(IL7R$clusterExprPerc, test$clusterExprPerc) & 
                          identical(IL7R$clusterAveExpr, test$clusterAveExpr))
  
  testthat::expect_error(runClusterSummaryMetrics(pbmc3k_2.7.1_sce, useAssay="logcounts", "IL7R", clusters="howdy"),
                         "Specified variable 'howdy' not found in colData(inSCE)", fixed=TRUE)
  
  testthat::expect_error(runClusterSummaryMetrics(pbmc3k_2.7.1_sce, useAssay="logcounts", 
                                                  gene=c("IL7R", "CD3E", "applesauce", "cowboy"), clusters="SingleR_hpca_fine_pruned.labels"),
                         "Specified genes 'applesauce, cowboy' not found in rowData(inSCE)$Symbol", fixed=TRUE)
})