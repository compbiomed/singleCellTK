#' @export
differentialAbundanceFET <- function(inSCE, variable, phenotype, control, case) {
    #if (method == "celda"){
    #    cluster <- celdaClusters(inSCE)
    #}else if(method == "seurat"){
    #    cluster <- inSCE$Seurat_selected_res
    #}

    #variable <- inSCE$sample
    cluster <- colData(inSCE)[, variable]
    variable <- colData(inSCE)[, phenotype]

    control.lab <- paste(control, collapse=",")
    case.lab <- paste(case, collapse=",")
    label <- rep(NA, length(variable))
    label[variable %in% case] <- case.lab
    label[variable %in% control] <- control.lab
    label <- factor(label, levels=c(control.lab, case.lab))

    res <- matrix(NA, nrow=length(unique(cluster)), ncol=10)
    cluster.label <- sort(unique(cluster))
    for(i in seq_along(cluster.label)) {
        cluster.factor <- factor(ifelse(cluster %in% cluster.label[i], cluster.label[i], "Other"), levels=c(i, "Other"))

        ta <- table(cluster.factor, label)
        ta.p <- prop.table(ta, 2)
        fet <- fisher.test(ta)
        res[i,] <- c(ta, ta.p, fet$estimate, fet$p.value)
    }
    colnames(res) <- c(paste0("Number of cells in cluster and in ", control.lab),
                       paste0("Number of cells NOT in cluster and in ", control.lab),
                       paste0("Number of cells in cluster and in ", case.lab),
                       paste0("Number of cells NOT in cluster and in ", case.lab),
                       paste0("Fraction of cells in cluster and in ", control.lab),
                       paste0("Fraction of cells NOT in cluster and in ", control.lab),
                       paste0("Fraction of cells in cluster and in ", case.lab),
                       paste0("Fraction of cells NOT in cluster and in ", case.lab),
                       "Odds_Ratio", "Pvalue")
    res <- data.frame(Cluster=cluster.label, res, FDR=p.adjust(res[,"Pvalue"], 'fdr'), check.names=FALSE)
    return(res)
}

#' @export
plotClusterAbundance <- function(inSCE, variable, phenotype) {
    #if (method == "celda"){
    #    cluster <- celdaClusters(inSCE)
    #    color_palette <- celda::distinctColors(length(unique(sce$Celda)))
    #}else if(method == "seurat"){
    #    cluster <- inSCE$Seurat_selected_res
    #    palette_data <- ggplot2::ggplot_build(UMAP_coloredby_cluster)$data[[1]]
    #    color_palette <- unique((palette_data[order(palette_data$group), ])$colour)
    #}

    cluster <- colData(inSCE)[, variable]
    color_palette <- celda::distinctColors(length(unique(cluster)))

    #label <- inSCE$sample
    label <- colData(inSCE)[, phenotype]

    cluster.color <- color_palette
    df <- data.frame(Cluster=as.factor(cluster), Sample=as.factor(label))

    g1 <- ggplot(df, aes(Cluster)) + geom_bar(aes(fill=Sample)) +
        theme_classic() +
        scale_y_continuous(expand = c(0, 0))

    g2 <- ggplot(df, aes(Cluster)) + geom_bar(aes(fill=Sample), position="fill") +
        theme_classic() +
        scale_y_continuous(expand = c(0, 0)) +
        ylab("Fraction")

    g3 <- ggplot(df, aes(Sample)) + geom_bar(aes(fill=Cluster)) +
        theme_classic() +
        scale_y_continuous(expand = c(0, 0)) +
        theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
        scale_fill_manual(values=cluster.color)

    g4 <- ggplot(df, aes(Sample)) + geom_bar(aes(fill=Cluster), position="fill") +
        theme_classic() +
        scale_y_continuous(expand = c(0, 0)) +
        theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
        scale_fill_manual(values=cluster.color) +
        ylab("Fraction")
    g <- list(g1, g2, g3, g4)
    return(g)
}




