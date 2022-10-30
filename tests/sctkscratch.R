library(singleCellTK)
sce <- runNormalization(sce, useAssay = "counts", outAssayName = "logcounts",
                        normalizationMethod = "logNormCounts")

topMarkers <- getFindMarkerTopTable(sce, topN = 1, log2fcThreshold = 0, fdrThreshold = 0.05,
                                 minClustExprPerc = 0.5, maxCtrlExprPerc = 0.5,
                                 minMeanExpr = 0)
IL7R <- dplyr::select(dplyr::filter(topMarkers, Gene=="IL7R"), clusterExprPerc, clusterAveExpr)
rownames(IL7R) <- NULL

test <- runClusterSummaryMetrics(pbmc3k_2.7.1_sce, useAssay="logcounts", "IL7R", clusters="cluster")
test <- dplyr::select(dplyr::filter(test, cluster==1), clusterExprPerc, clusterAveExpr)
rownames(test) <- NULL
test$clusterExprPerc <- as.numeric(test$clusterExprPerc)
identical(IL7R$clusterAveExpr, test$clusterAveExpr)
data("pbmc3k_2.7.1_sce")
data("scExample")

B2M <- runClusterSummaryMetrics(sce, useAssay="counts", c("B2M"), displayName="feature_name", clusters="type")
percExpr <- c(1.0000000, 0.7282051, 1.0000000)
aveExpr <- c(94.311111, 2.517949, 54.473333)


y# See number of cells after filtering
dimnames(counts(sce))[[1]] <- rowData(sce)$feature_name
test <- runClusterSummaryMetrics(sce, useAssay="counts", c("cfvgbhn", "B2M"), displayName="feature_name", clusters="type")

test <- scuttle::aggregateAcrossCells(sce, ids=SingleCellExperiment::colData(sce)[,"type"], 
                              statistics="prop.detected", use.assay.type="counts", 
                              subset.row="B2M")
test <- assay(sce, use.assay.type="counts")
test2 <- assay(pbmc3k_2.7.1_sce, use.assay.type="logcounts")
test <- scuttle:::.summarize_assay(assay(sce, I="counts"), ids=SingleCellExperiment::colData(sce)[,"type"], statistics="prop.detected", subset.row="B2M")
                              
plotBubble(pbmc3k_2.7.1_sce, gene="CD3E", displayName="Symbol_TENx", title="cluster test")
plotBubble(pbmc3k_2.7.1_sce, gene=c("applesauce", "cowboy"), displayName="Symbol_TENx", clusters="SingleR_hpca_fine_pruned.labels", title="cell type test")
