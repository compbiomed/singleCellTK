# Internal helpers built on scrapper. They replace functions deprecated in
# scuttle 1.22 and scran 1.40 while keeping singleCellTK's outputs unchanged.

# Mean expression and proportion of cells with expression > 0 for each group
# of cells. Returns a list of two features x groups matrices, with the
# feature names as rownames and the sorted group labels as colnames, like
# scuttle::aggregateAcrossCells() with statistics "mean" and "prop.detected".
.aggregateMeanDetected <- function(mat, ids) {
  agg <- scrapper::aggregateAcrossCells(mat, factors = list(ids))
  groupNames <- as.character(agg$combinations[[1]])
  counts <- agg$counts
  avg <- sweep(agg$sums, 2, counts, "/")
  det <- sweep(agg$detected, 2, counts, "/")
  dimnames(avg) <- list(rownames(mat), groupNames)
  dimnames(det) <- list(rownames(mat), groupNames)
  list(mean = avg, prop.detected = det)
}

# Per-cluster mean of a reducedDim, with one row per cluster in factor
# order (levels for a factor, otherwise sorted values), like reducedDim() of
# scuttle::aggregateAcrossCells().
.clusterCentroids <- function(inSCE, clusters, useReducedDim) {
  emb <- SingleCellExperiment::reducedDim(inSCE, useReducedDim)
  ids <- factor(clusters)
  rowsum(emb, ids) / as.vector(table(ids))
}

# Shared-nearest-neighbor graph of cells, as scran::buildSNNGraph() built it
# before it was deprecated. For a features x cells assay, cells are first
# projected onto the top nComp principal components of the centered data
# (skipped when nComp is NA or not smaller than the number of features).
# For a cells x dimensions embedding (transposed = TRUE), it is used as is.
.snnGraph <- function(mat, k, weightType, nComp = NA, transposed = FALSE,
                      BPPARAM = BiocParallel::SerialParam()) {
  if (!transposed) {
    mat <- t(mat)
    if (!is.na(nComp) && nComp < ncol(mat)) {
      svd <- BiocSingular::runSVD(mat, k = nComp, nu = nComp, nv = 0,
                                  center = TRUE,
                                  BSPARAM = BiocSingular::bsparam(),
                                  BPPARAM = BPPARAM)
      mat <- sweep(svd$u, 2, svd$d, "*")
    }
  }
  bluster::makeSNNGraph(mat, k = k, type = weightType, BPPARAM = BPPARAM)
}
