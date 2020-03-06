#' Summarize SCtkExperiment
#'
#' Creates a table of summary metrics from an input SCtkExperiment.
#'
#' @param inSCE Input SCtkExperiment object. Required
#' @param useAssay Indicate which assay to summarize. Default is "counts"
#' @param expressionCutoff Count number of samples with fewer than
#' expressionCutoff genes. The default is 1700.
#'
#' @return A data.frame object of summary metrics.
#' @export
#' @examples
#' data("mouseBrainSubsetSCE")
#' summarizeTable(mouseBrainSubsetSCE)
#'
summarizeTable <- function(inSCE, useAssay="counts", expressionCutoff=1700){
  return(
    data.frame(
      "Metric" = c(
        "Number of Samples",
        "Number of Genes",
        "Average number of reads per cell",
        "Average number of genes per cell",
        paste0("Samples with <", expressionCutoff, " detected genes"),
        "Genes with no expression across all samples"
      ),
      "Value" = c(
        ncol(inSCE),
        nrow(inSCE),
        as.integer(mean(DelayedArray::colSums(
          SummarizedExperiment::assay(inSCE, useAssay)))),
        as.integer(mean(DelayedArray::colSums(
          SummarizedExperiment::assay(inSCE, useAssay) > 0))),
        sum(DelayedArray::colSums(
          SummarizedExperiment::assay(inSCE, useAssay) != 0) <
            expressionCutoff),
        sum(DelayedArray::rowSums(
          SummarizedExperiment::assay(inSCE, useAssay)) == 0)
      )
    )
  )
}

#' Create a SCtkExperiment object
#'
#' From a file of counts and a file of annotation information, create a
#' SCtkExperiment object.
#'
#' @param assayFile The path to a text file that contains a header row of sample
#' names, and rows of raw counts per gene for those samples.
#' @param annotFile The path to a text file that contains columns of annotation
#' information for each sample in the assayFile. This file should have the same
#' number of rows as there are columns in the assayFile.
#' @param featureFile The path to a text file that contains columns of
#' annotation information for each gene in the count matrix. This file should
#' have the same genes in the same order as assayFile. This is optional.
#' @param assayName The name of the assay that you are uploading. The default
#' is "counts".
#' @param inputDataFrames If TRUE, assayFile and annotFile are read as data
#' frames instead of file paths. The default is FALSE.
#' @param createLogCounts If TRUE, create a log2(counts+1) normalized assay
#' and include it in the object. The default is TRUE
#'
#' @return a SCtkExperiment object
#' @export
#' @examples
#' data("mouseBrainSubsetSCE")
#' counts_mat <- assay(mouseBrainSubsetSCE, "counts")
#' sample_annot <- colData(mouseBrainSubsetSCE)
#' row_annot <- rowData(mouseBrainSubsetSCE)
#' newSCE <- createSCE(assayFile = counts_mat, annotFile = sample_annot,
#'                     featureFile = row_annot, assayName = "counts",
#'                     inputDataFrames = TRUE, createLogCounts = TRUE)
createSCE <- function(assayFile=NULL, annotFile=NULL, featureFile=NULL,
                      assayName="counts", inputDataFrames=FALSE,
                      createLogCounts=TRUE){
  
  if (is.null(assayFile)){
    stop("You must supply a count file.")
  }
  if (inputDataFrames){
    countsin <- assayFile
    annotin <- annotFile
    featurein <- featureFile
  } else{
    countsin <- utils::read.table(assayFile, sep = "\t", header = TRUE,
                                  row.names = 1)
    if (!is.null(annotFile)){
      annotin <- utils::read.table(annotFile, sep = "\t", header = TRUE,
                                   row.names = 1)
    }
    if (!is.null(featureFile)){
      featurein <- utils::read.table(featureFile, sep = "\t", header = TRUE,
                                     row.names = 1)
    }
  }
  if (is.null(annotFile)){
    annotin <- data.frame(row.names = colnames(countsin))
    annotin$Sample <- rownames(annotin)
    annotin <- S4Vectors::DataFrame(annotin)
  }
  if (is.null(featureFile)){
    featurein <- data.frame(Gene = rownames(countsin))
    rownames(featurein) <- featurein$Gene
    featurein <- S4Vectors::DataFrame(featurein)
  }
  if (nrow(annotin) != ncol(countsin)){
    stop("Different number of samples in input matrix and annotations: annot: ",
         nrow(annotin), ", counts: ", ncol(countsin))
  }
  if (nrow(featurein) != nrow(countsin)){
    stop("Different number of samples in input matrix and feature annotation",
         nrow(featurein), ", counts: ", nrow(countsin))
  }
  if (any(rownames(annotin) != colnames(countsin))){
    stop("Sample names in input matrix and annotation do not match!\nExample: ",
         rownames(annotin)[rownames(annotin) != colnames(countsin)][1], " vs. ",
         colnames(countsin)[rownames(annotin) != colnames(countsin)][1])
  }
  if (any(rownames(featurein) != rownames(countsin))){
    stop("Sample names in input matrix and feature annotation do not match!")
  }
  assaylist <- list()
  assaylist[[assayName]] <- as.matrix(countsin)
  newassay <- SCtkExperiment(assays = assaylist,
                             colData = annotin,
                             rowData = featurein)
  if (createLogCounts){
    SummarizedExperiment::assay(newassay, paste0("log", assayName)) <-
      log2(SummarizedExperiment::assay(newassay, assayName) + 1)
  }
  return(newassay)
}

#' Filter Genes and Samples from a Single Cell Object
#'
#' @param inSCE Input SCtkExperiment object. Required
#' @param useAssay Indicate which assay to use for filtering. Default is
#' "counts"
#' @param deletesamples List of samples to delete from the object.
#' @param removeNoExpress Remove genes that have no expression across all
#' samples. The default is true
#' @param removeBottom Fraction of low expression genes to remove from the
#' single cell object. This occurs after removeNoExpress. The default is 0.50.
#' @param minimumDetectGenes Minimum number of genes with at least 1
#' count to include a sample in the single cell object. The default is 1700.
#' @param filterSpike Apply filtering to Spike in controls (indicated by
#' isSpike).
#' The default is TRUE.
#'
#' @return The filtered single cell object.
#' @export
#' @examples
#' data("mouseBrainSubsetSCE")
#' mouseBrainSubsetSCE <- filterSCData(mouseBrainSubsetSCE,
#'                                     deletesamples="X1772063061_G11")
filterSCData <- function(inSCE, useAssay="counts", deletesamples=NULL,
                         removeNoExpress=TRUE, removeBottom=0.5,
                         minimumDetectGenes=1700, filterSpike=TRUE){
  #delete specified samples
  inSCE <- inSCE[, !(colnames(inSCE) %in% deletesamples)]

  if (filterSpike){
    nkeeprows <- ceiling((1 - removeBottom) * as.numeric(nrow(inSCE)))
    tokeeprow <- order(rowSums(SummarizedExperiment::assay(inSCE, useAssay)),
                       decreasing = TRUE)[seq_len(nkeeprows)]
  } else {
    nkeeprows <- ceiling((1 - removeBottom) * as.numeric(nrow(inSCE))) -
      sum(SingleCellExperiment::isSpike(inSCE))
    tokeeprow <- order(rowSums(SummarizedExperiment::assay(inSCE, useAssay)),
                       decreasing = TRUE)
    tokeeprow <- setdiff(tokeeprow,
                         which(SingleCellExperiment::isSpike(inSCE)))
    tokeeprow <- tokeeprow[seq_len(nkeeprows)]
    tokeeprow <- c(tokeeprow, which(SingleCellExperiment::isSpike(inSCE)))
  }
  tokeepcol <- colSums(SummarizedExperiment::assay(inSCE, useAssay) != 0) >=
    minimumDetectGenes
  inSCE <- inSCE[tokeeprow, tokeepcol]

  #remove genes with no expression
  if (removeNoExpress){
    if (filterSpike){
      inSCE <- inSCE[rowSums(SummarizedExperiment::assay(inSCE,
                                                              useAssay)) != 0, ]
    } else {
      inSCE <- inSCE[(rowSums(
        SummarizedExperiment::assay(inSCE, useAssay)) != 0 |
          SingleCellExperiment::isSpike(inSCE)), ]
    }
  }

  return(inSCE)
}

#' Generate a distinct palette for coloring different clusters
#'
#' @param n Integer; Number of colors to generate
#' @param hues Character vector of R colors available from the colors()
#' function. These will be used as the base colors for the clustering scheme.
#' Different saturations and values (i.e. darkness) will be generated for each
#' hue.
#' @param saturation.range Numeric vector of length 2 with values between 0 and
#' 1. Default: c(0.25, 1)
#' @param value.range Numeric vector of length 2 with values between 0 and 1.
#' Default: c(0.5, 1)
#' @return A vector of distinct colors that have been converted to  HEX from
#' HSV.
#' @export
#' @examples
#' distinctColors(10)
distinctColors <- function(n, hues = c("red", "cyan", "orange", "blue",
                                        "yellow", "purple", "green", "magenta"),
                           saturation.range = c(0.7, 1),
                           value.range = c(0.7, 1)) {
  #Adapted from compbiomed/celda, thanks to all celda developers
  if (!(all(hues %in% grDevices::colors()))) {
    stop("Only color names listed in the 'color'",
         " function can be used in 'hues'")
  }

  ## Convert R colors to RGB and then to HSV color format
  hues.hsv <- grDevices::rgb2hsv(grDevices::col2rgb(hues))

  ## Calculate all combination of saturation/value pairs
  ## Note that low saturation with low value (i.e. high darkness) is too dark
  ## for all hues
  ## Likewise, high saturation with high value (i.e. low darkness) is hard to
  ## distinguish
  ## Therefore, saturation and value are set to be anticorrelated
  num.vs <- ceiling(n / length(hues))
  s <- seq(from = saturation.range[1], to = saturation.range[2],
           length = num.vs)
  v <- seq(from = value.range[2], to = value.range[1], length = num.vs)

  ## Create all combination of hues with saturation/value pairs
  new.hsv <- c()
  for (i in seq_len(num.vs)) {
    temp <- rbind(hues.hsv[1, ], s[i], v[i])
    new.hsv <- cbind(new.hsv, temp)
  }

  ## Convert to hex
  col <- grDevices::hsv(new.hsv[1, ], new.hsv[2, ], new.hsv[3, ])

  return(col[seq_len(n)])
}

#test shiny functions
.testFunctions <- function(){
  if (interactive()){
    res <- DT::datatable(matrix(1, 2))
    shinyjs::runExample("basic")
    shinyalert::runExample()
    p <- plotly::plot_ly(data = data.frame(test = c(1, 2, 3)),
                         x = "test", type = "histogram")
    colourpicker::runExample()
    rt <- ape::rtree(10)
    gt <- ggtree::ggtree(rt)
    shinycssloaders::withSpinner(shiny::plotOutput("my_plot"))
    x <- rbind(cbind(stats::rnorm(200, 0, 8), stats::rnorm(200, 0, 8)),
               cbind(stats::rnorm(300, 50, 8), stats::rnorm(300, 50, 8)))
    clarax <- cluster::clara(x, 2, samples = 50)
    circlize::colorRamp2(c(1, 2, 3), c("red", "blue", "black"))
  }
}
# 
# log = function(f){
#   function(...){
#     args <- listToString(list(...))
#     fName <- as.character(match.call()[1])
#     res = f(...)
#     call <- paste0(fName, "(", args, ")")
#     # print("logging a function:")
#     print(noquote(call))
#     # insertUI(paste0("#", "console"), where = "beforeEnd",
#     # ui = tags$p(paste0(call, "\n", collapse = ""))
#     # )
#     return(res)
#   }
# }
# 
# listToString <- function(l) {
#   res = ""
#   for(i in l[1:length(l)]){
#     if (is.character(i)) {
#       res = paste0(res, sQuote(i), ",")
#     } else {
#       res = paste0(res, i, ",")
#     }
#   }
#   return(substr(res, 1, nchar(res)-1))
# }

