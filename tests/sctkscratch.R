library(singleCellTK)
data("pbmc3k_2.7.1_sce")

plotBubble(pbmc3k_2.7.1_sce, gene=c("IL7R", "CD3E"), title="cluster test", clusters="SingleR_hpca_fine_pruned.labels")
plotBubble(pbmc3k_2.7.1_sce, gene=c("IL7R", "CD3E", "applesauce", "cowboy"), clusters="SingleR_hpca_fine_pruned.labels", title="cell type test")
