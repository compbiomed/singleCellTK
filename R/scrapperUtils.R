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
