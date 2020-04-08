.violinPlotQC <- function(inSCE,
                         output,
                         sample,
                         useAssay,
                         fileName,
                         fileWidth,
                         fileHeight,
                         fileUnits,
                         figRows,
                         figCols){
    samples <- unique(sample)
    for(i in seq_len(length(samples))){
        sceSampleInd <- which(sample == samples[i])
        sceSample <- inSCE[, sceSampleInd]
        mat <- SummarizedExperiment::assay(sceSample, i = useAssay)
        value <- colData(sceSample)[, output]
        
        df <- data.frame(x = "Sample", value = value)
        
        p <- ggplot2::ggplot(df) +
            ggplot2::aes_string(
                x = "x",
                y = "value") +
            ggplot2::geom_violin(trim = TRUE, scale = "width") +
            ggplot2::geom_jitter(height = 0, size = 0.1) +
            ggplot2::ggtitle(output) +
            ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.5, size = 14, face = "bold"),
                           axis.title = ggplot2::element_blank(),
                           axis.text.x = ggplot2::element_blank(),  ## Since the function won't know what the units of each QC metric will be
                           panel.background = ggplot2::element_blank(),
                           panel.border = ggplot2::element_blank())
        return(p)
    }
}

#' @title Violin plot of QC outputs.
#' @description Visualizes QC output values via a violin plot.
#'  Users may either define specific outputs to plot, or can specify
#'  all outputs from a specific QC algorithm. Will save plot if parameter 
#'  \code{fileName} is defined.
#' @param inSCE A \link[SingleCellExperiment]{SingleCellExperiment} object.
#'  Needs \code{counts} in assays slot.
#' @param outputnames Names of QC values in the colData that will be plotted.
#' @param sample Character vector. Indicates which sample each cell belongs to.
#' @param useAssay Indicate which assay to use. Default "counts".
#' @param fileName Desired file name of the output plot file. (e.g. .pdf, .png) 
#'  If NULL, a plot object will be outputted, with no file output. Default NULL.
#' @param fileWidth Width of output plot file.
#' @param fileHeight Height of output plot file.
#' @param fileUnits Units of output plot file dimensions. Default "cm".
#' @param figRows Number of rows in figure.
#' @param figCols Number of columns in figure.
#' @examples
#' data("mouseBrainSubsetSCE")
#' violinPlotQC(inSCE = mouseBrainSubsetSCE, 
#'  outputnames = "total.mRNA.mol")
#' @export
violinPlotQC <- function(inSCE,
    outputnames,
    sample = NULL,
    useAssay = "counts",
    fileName = NULL,
    fileWidth = 20,
    fileHeight = 20,
    fileUnits = "cm",
    figRows = NULL,
    figCols = NULL){

    
    if(!is.null(sample)) {
        if(length(sample) != ncol(inSCE)) {
        stop("'sample' must be the same length as the number of columns in 'inSCE'")
        }
    } else {
        sample <- rep(1, ncol(inSCE))
    }

    if(!is.null(outputnames)){
        all.output.names <- names(colData(inSCE))
        nomatch <- outputnames[!outputnames %in% all.output.names]
        
        if (length(nomatch) > 0) {
            stop("'", paste(nomatch, collapse=","), "' is not found in ColData.")
        }
        
        output.list <- outputnames
    }else{
        stop("You must either define the desired output.")
    }
    
    list.plot <- list()
    
    for(n in 1:length(output.list)){
        output <- output.list[n]
        p <- .violinPlotQC(inSCE = inSCE,
                           sample = sample,
                           output = output,
                           useAssay = useAssay,
                           fileName = fileName,
                           fileWidth = fileWidth,
                           fileHeight = fileHeight,
                           fileUnits = fileUnits,
                           figRows = figRows,
                           figCols = figCols)
        
        list.plot[[n]] <- p   
    }

    v.plot <- gridExtra::grid.arrange(grobs = list.plot, nrow = figRows, ncol = figCols)
    if(!is.null(fileName)){
        ggplot2::ggsave(file = fileName, plot = v.plot, width = fileWidth, height = fileHeight, units = fileUnits)
    }
    
    return(v.plot)

}

