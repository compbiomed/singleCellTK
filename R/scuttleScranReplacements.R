# Internal replacements for functions deprecated in scuttle 1.22 and scran
# 1.40, built on scrapper, bluster, and base R. They keep singleCellTK's outputs
# unchanged.

# Mean expression and proportion of cells with expression > 0 for each group
# of cells. Returns a list of two features x groups matrices, with the
# feature names as rownames and the sorted group labels as colnames, like
# scuttle::aggregateAcrossCells() with statistics "mean" and "prop.detected".
# Cells with a missing group are dropped, as scuttle did.
.aggregateMeanDetected <- function(mat, ids) {
  keep <- !is.na(ids)
  agg <- scrapper::aggregateAcrossCells(mat[, keep, drop = FALSE],
                                        factors = list(ids[keep]))
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
# scuttle::aggregateAcrossCells(). Cells without a cluster are dropped.
.clusterCentroids <- function(inSCE, clusters, useReducedDim) {
  emb <- SingleCellExperiment::reducedDim(inSCE, useReducedDim)
  keep <- !is.na(clusters)
  ids <- factor(clusters[keep])
  rowsum(emb[keep, , drop = FALSE], ids) / as.vector(table(ids))
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

# Two-sided Wilcoxon rank-sum test of each row of mat between the cells in
# ix1 and ix2 (logical or integer indices), with the normal approximation, tie correction, and continuity
# correction, as in stats::wilcox.test(exact = FALSE, correct = TRUE) and
# scran::pairwiseWilcox(), which is deprecated without a replacement. Rows
# without variation get p = 1. Returns a DataFrame with p.value and FDR
# (Benjamini-Hochberg), one row per row of mat. Rows are processed in chunks
# so that only chunkSize rows are densified at a time.
.wilcoxTest <- function(mat, ix1, ix2, chunkSize = 1000) {
  if (is.logical(ix1)) ix1 <- which(ix1)
  if (is.logical(ix2)) ix2 <- which(ix2)
  n1 <- length(ix1)
  n2 <- length(ix2)
  n <- n1 + n2
  pValue <- numeric(nrow(mat))
  for (start in seq(1, nrow(mat), by = chunkSize)) {
    rows <- seq(start, min(start + chunkSize - 1, nrow(mat)))
    x <- as.matrix(mat[rows, c(ix1, ix2), drop = FALSE])
    ranks <- matrixStats::rowRanks(x, ties.method = "average")
    tieSize <- matrixStats::rowRanks(x, ties.method = "max") -
      matrixStats::rowRanks(x, ties.method = "min") + 1
    stat <- rowSums(ranks[, seq_len(n1), drop = FALSE]) - n1 * (n1 + 1) / 2
    centered <- stat - n1 * n2 / 2
    tieTerm <- rowSums(tieSize^2 - 1)
    sigma <- sqrt((n1 * n2 / 12) * ((n + 1) - tieTerm / (n * (n - 1))))
    z <- (centered - sign(centered) * 0.5) / sigma
    p <- 2 * pmin(stats::pnorm(z), stats::pnorm(z, lower.tail = FALSE))
    p[sigma == 0] <- 1
    pValue[rows] <- p
  }
  S4Vectors::DataFrame(p.value = pValue,
                       FDR = stats::p.adjust(pValue, method = "BH"),
                       row.names = rownames(mat))
}
