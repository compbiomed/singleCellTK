
library(data.table)
library(Matrix)

matrixLocation <- "../../20200117_optimus_test/bb4a2a5e-ff34-41b6-97d2-0c0c0c534530/call-MergeCountFiles/sparse_counts.npz"
colIndexLocation <- "../../20200117_optimus_test/bb4a2a5e-ff34-41b6-97d2-0c0c0c534530/call-MergeCountFiles/sparse_counts_col_index.npy"
rowIndexLocation <- "../../20200117_optimus_test/bb4a2a5e-ff34-41b6-97d2-0c0c0c534530/call-MergeCountFiles/sparse_counts_row_index.npy"
CellMetricsLocation <- "../../20200117_optimus_test/bb4a2a5e-ff34-41b6-97d2-0c0c0c534530/call-MergeCellMetrics/merged-cell-metrics.csv.gz"
GeneMetricsLocation <- "../../20200117_optimus_test/bb4a2a5e-ff34-41b6-97d2-0c0c0c534530/call-MergeGeneMetrics/merged-gene-metrics.csv.gz"
EmptyDropsLocation <- "../../20200117_optimus_test/bb4a2a5e-ff34-41b6-97d2-0c0c0c534530/call-RunEmptyDrops/empty_drops_result.csv"

sparse <- reticulate::import("scipy.sparse")
np <- reticulate::import("numpy")

mat <- sparse$load_npz(matrixLocation)
colIndex <- np$load(colIndexLocation)
rowIndex <- np$load(rowIndexLocation)


# Only keep top 1000 cells
rorder <- order(rowSums(mat), decreasing = T)[1:1000]
matrorder <- mat[rorder, ]
corder <- order(colSums(matrorder), decreasing = T)[1:20]

matExample <- mat[rorder, ][, corder]
colIndexeExample <- colIndex[corder]
rowIndexeExample <- rowIndex[rorder]

sparse$save_npz("./Optimus/bb4a2a5e-ff34-41b6-97d2-0c0c0c534530/call-MergeCountFiles/sparse_counts.npz", matExample)
np$save("./Optimus/bb4a2a5e-ff34-41b6-97d2-0c0c0c534530/call-MergeCountFiles/sparse_counts_col_index.npy", colIndexeExample)
np$save("./Optimus/bb4a2a5e-ff34-41b6-97d2-0c0c0c534530/call-MergeCountFiles/sparse_counts_row_index.npy", rowIndexeExample)

cellMetrics <- data.table::fread(CellMetricsLocation)
cm <- merge(cellMetrics,
  data.table::data.table(V1 = rowIndexeExample),
  by.x = "V1",
  by.y = "V1",
  all.x = FALSE,
  all.y = TRUE,
  sort = FALSE)

fwrite(cm, file = "./Optimus/bb4a2a5e-ff34-41b6-97d2-0c0c0c534530/call-MergeCellMetrics/merged-cell-metrics.csv.gz")

rowMetrics <- data.table::fread(GeneMetricsLocation)
rm <- merge(rowMetrics,
  data.table::data.table(V1 = colIndexeExample),
  by.x = "V1",
  by.y = "V1",
  all.x = FALSE,
  all.y = TRUE,
  sort = FALSE)

fwrite(rm, file = "./Optimus/bb4a2a5e-ff34-41b6-97d2-0c0c0c534530/call-MergeGeneMetrics/merged-gene-metrics.csv.gz")

emptyDrops <- data.table::fread(EmptyDropsLocation)
ed <- merge(emptyDrops,
  data.table::data.table(CellId = rowIndexeExample),
  by.x = "CellId",
  by.y = "CellId",
  all.x = FALSE,
  all.y = TRUE,
  sort = FALSE)

fwrite(ed, file = "./Optimus/bb4a2a5e-ff34-41b6-97d2-0c0c0c534530/call-RunEmptyDrops/empty_drops_result.csv")




