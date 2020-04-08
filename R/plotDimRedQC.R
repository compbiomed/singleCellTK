.plotDimRedQC <- function(inSCE, algo,
                     reducedDimName, runDimensionReduction,
                     useAssay, size){
    
    if(!(reducedDimName %in% names(SingleCellExperiment::reducedDims(inSCE)))){
        if (runDimensionReduction){
            if(reducedDimName == "TSNE"){
                inSCE <- getTSNE(inSCE, useAssay = useAssay,
                             reducedDimName = reducedDimName)
            }else if(reducedDimName == "UMAP"){
                inSCE <- getUMAP(inSCE, useAssay = useAssay,
                                 reducedDimName = reducedDimName)
            }
        } else {
            stop(reducedDimName,
                 " dimension not found. Run dimension reduction via getTSNE() or getUMAP(). 
                 Alternatively, set runDimensionReduction to TRUE.")
        }
    }
    ix <- which(names(colData(inSCE)) == algo)
    label <- colData(inSCE)[,ix]
    
    dimRedDf <- data.frame(SingleCellExperiment::reducedDim(inSCE,
                                                          reducedDimName),label)
    if (ncol(dimRedDf) > 3){
        warning("More than two dimensions. Using the first two.")
    }
    
    if(reducedDimName == "TSNE"){
        colnames(dimRedDf)[1] <- "TSNE1"
        xlab <- "TSNE1"
        colnames(dimRedDf)[2] <- "TSNE2"
        ylab <- "TSNE2"
    }else if(reducedDimName == "UMAP"){
        colnames(dimRedDf)[1] <- "UMAP1"
        xlab <- "UMAP1"
        colnames(dimRedDf)[2] <- "UMAP2"
        ylab <- "UMAP2"
    }
    
    colnames(dimRedDf)[3] <- "Label"
    xdim <- colnames(dimRedDf)[1]
    ydim <- colnames(dimRedDf)[2]
    
    
    g <- ggplot2::ggplot(dimRedDf, ggplot2::aes_string(x = xlab, y = ylab)) +
        ggplot2::geom_point(stat = "identity",
                            size = size,
                            ggplot2::aes_string(color = "Label")) +
        ggplot2::theme(
            panel.grid.major = ggplot2::element_blank(),
            panel.grid.minor = ggplot2::element_blank(),
            panel.background = ggplot2::element_blank(),
            axis.line = ggplot2::element_line(color = "black"))# +
    
    return(g)
}

#' @title Dimension reduction plot of QC outputs.
#' @description Visualizes QC output values via a scatterplot.
#'  Users may either define specific outputs to plot, or can specify
#'  all outputs from a specific QC algorithm. Currently supports TSNE and UMAP.
#'  \code{fileName} is defined.
#' @param inSCE A \link[SingleCellExperiment]{SingleCellExperiment} object.
#'  Needs \code{counts} in assays slot.
#' @param outputs Names of QC values in the colData that will be plotted.
#' @param log.counts Determines whether counts should be log transformed 
#'  prior to dimension reduction. Default TRUE. 
#' @param reducedDimName Determines which dimension reduction coordinates 
#'  to plot. 
#' @param runDimensionReduction Determines whether dimension reduction needs
#'  to be run prior to plotting. Default FALSE. If TRUE, will run dimension
#'  reduction algorithm specified in "reducedDimName".
#' @param useAssay Indicate which assay to use. Default "counts".
#' @param size Numeric. Size of dots in plot. Default 0.5.
#' @param figRows Number of rows in figure.
#' @param figCols Number of columns in figure.
#' @examples
#' data("mouseBrainSubsetSCE")
#' plotDimRedQC(inSCE = mouseBrainSubsetSCE, outputs = c("tissue","group"), 
#' runDimensionReduction=TRUE)
#' @export
plotDimRedQC <- function(inSCE, outputs, log.counts = TRUE,
                       reducedDimName="TSNE", runDimensionReduction=FALSE,
                       useAssay="counts", size = 0.5, figRows = NULL,
                       figCols = NULL){
    dimred.list <- list()
    
    if(log.counts){
        assays(inSCE)$logcounts <- log(as.matrix(assays(inSCE)$counts)+1)
        useAssay = "logcounts"
    }
    
    for(algo in outputs){
      dimred <- .plotDimRedQC(inSCE = inSCE, 
                            algo = algo,
                            reducedDimName = reducedDimName,
                            runDimensionReduction = runDimensionReduction,
                            useAssay = useAssay,
                            size = size)
      dimred.list[[algo]] = dimred
    }
    v.plot <- gridExtra::grid.arrange(grobs = dimred.list, nrow = figRows, ncol = figCols)
    return(v.plot)
}

