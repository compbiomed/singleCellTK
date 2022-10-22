library(singleCellTK)
load("C:/Users/creli/PycharmProjects/BU/Campbell/singleCellTK/data/pbmc3k_2.7.1_sce.rda")

plotBubble(pbmc3k_2.7.1_sce, gene="CD3E", title="cluster test")
plotBubble(pbmc3k_2.7.1_sce, gene=c("IL7R", "CD3E", "applesauce", "cowboy"), clusters="SingleR_hpca_fine_pruned.labels", title="cell type test")
