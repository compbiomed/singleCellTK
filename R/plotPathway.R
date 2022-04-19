#' List pathway analysis result names
#' @details Pathway analysis results will be stored as matrices in 
#' \code{reducedDims} slot of \code{inSCE}. This function lists the result names
#' stored in \code{metadata} slot when analysis is performed. 
#' @param inSCE Input \linkS4class{SingleCellExperiment} object.
#' @param stopIfNone Whether to stop and raise an error if no results found. If
#' \code{FALSE}, will return an empty character vector. 
#' @return A character vector of valid pathway analysis result names.
#' @export
getPathwayResultNames <- function(inSCE, stopIfNone = FALSE){
    if (!"pathwayAnalysisResultNames" %in% names(S4Vectors::metadata(inSCE))) {
        if (isTRUE(stopIfNone)) {
            stop("No pathway analysis has been performed via singleCellTK.",
                 "Please try `runVAM()` or `runGSVA()`.")
        } else {
            warning("No pathway analysis has been performed via singleCellTK.",
                    "Please try `runVAM()` or `runGSVA()`.")
            return(character())
        }
    } else {
        return(S4Vectors::metadata(inSCE)[["pathwayAnalysisResultNames"]])
    }
}

#' Generate violin plots for pathway analysis results
#' @details \code{runGSVA()} or \code{runVAM()} should be applied in advance of 
#' using this function. Users can group the data by specifying \code{groupby}.
#' @param inSCE Input \linkS4class{SingleCellExperiment} object. With 
#' \code{runGSVA()} or \code{runVAM()} applied in advance. 
#' @param resultName A single character of the name of a score matrix, which 
#' should be found in \code{getPathwayResultNames(inSCE)}. 
#' @param geneset A single character specifying the geneset of interest. Should
#' be found in the geneSetCollection used for performing the analysis.
#' @param groupBy Either a single character specifying a column of 
#' \code{colData(inSCE)} or a vector of equal length as the number of cells. 
#' Default \code{NULL}.
#' @param boxplot Boolean, Whether to add a boxplot. Default \code{FALSE}.
#' @param violin Boolean, Whether to add a violin plot. Default \code{TRUE}.
#' @param dots Boolean, If \code{TRUE}, will plot dots for each violin plot.
#' Default \code{TRUE}.
#' @param summary Adds a summary statistic, as well as a crossbar to the violin 
#' plot. Options are \code{"mean"} or \code{"median"}, and \code{NULL} for not 
#' adding. Default \code{"median"}.
#' @param axisSize Size of x/y-axis ticks. Default \code{10}.
#' @param axisLabelSize Size of x/y-axis labels. Default \code{10}.
#' @param dotSize Size of dots. Default \code{0.5}.
#' @param transparency Transparency of the dots, values will be 0-1. Default 
#' \code{1}.
#' @param defaultTheme Removes grid in plot and sets axis title size to 
#' \code{10} when \code{TRUE}. Default \code{TRUE}.
#' @param gridLine Adds a horizontal grid line if \code{TRUE}. Will still be 
#' drawn even if \code{defaultTheme} is \code{TRUE}. Default \code{FALSE}.
#' @param title Title of plot. Default using \code{geneset}.
#' @param titleSize Size of the title of the plot. Default \code{15}.
#' @return A \code{ggplot} object for the violin plot
#' @export
#' @examples 
#' data("scExample", package = "singleCellTK")
#' sce <- subsetSCECols(sce, colData = "type != 'EmptyDroplet'")
#' sce <- scaterlogNormCounts(sce, assayName = "logcounts")
#' gs1 <- rownames(sce)[seq(10)]
#' gs2 <- rownames(sce)[seq(11,20)]
#' gs <- list("geneset1" = gs1, "geneset2" = gs2)
#' sce <- importGeneSetsFromList(inSCE = sce, geneSetList = gs,
#'                               by = "rownames")
#' sce <- runVAM(inSCE = sce, geneSetCollectionName = "GeneSetCollection",
#'               useAssay = "logcounts")
#' plotPathway(sce, "VAM_GeneSetCollection_CDF", "geneset1")
plotPathway <- function(inSCE, 
                        resultName, 
                        geneset,
                        groupBy = NULL,
                        boxplot = FALSE,
                        violin = TRUE,
                        dots = TRUE,
                        summary = "median",
                        axisSize = 10,
                        axisLabelSize = 10,
                        dotSize = 0.5,
                        transparency = 1,
                        defaultTheme = TRUE,
                        gridLine = FALSE,
                        title = geneset,
                        titleSize = NULL
){ 
    availResults <- getPathwayResultNames(inSCE, stopIfNone = TRUE)
    if (!resultName %in% availResults) {
        stop("Specified resultName not identified as a pathway analysis ", 
             "result from singleCellTK. Use `getPathwayResultNames(inSCE)` ", 
             "to see available options.")
    } else {
        if (!resultName %in% SingleCellExperiment::reducedDimNames(inSCE)) {
            stop("Pathway analysis result identified, but missing in ", 
                 "reducedDims(inSCE)")
        }
    }
    genesetCollectionName <- strsplit(resultName,"_")[[1]][2]
    geneSetCollection <- .getGeneSetCollection(inSCE, genesetCollectionName)
    if(!geneset %in% names(geneSetCollection)) {
        stop('"', geneset, '" is not a geneset from the geneSetCollection: "',
             genesetCollectionName, '". Try `getGenesetNamesFromCollection(inSCE, "',
             genesetCollectionName, '")` for available options.')
    }
    xlab <- ifelse(is.null(groupBy), "Sample", groupBy)
    plotSCEViolin(inSCE, 
                  slotName = "reducedDims", 
                  itemName = resultName, 
                  dimension = geneset , 
                  xlab = xlab, 
                  ylab = resultName, 
                  sample = NULL, 
                  groupBy = groupBy, 
                  boxplot = boxplot,
                  dots = dots,
                  violin = violin, 
                  summary = summary,
                  axisSize = axisSize,
                  axisLabelSize = axisLabelSize,
                  dotSize = dotSize,
                  transparency = transparency,
                  defaultTheme = defaultTheme,
                  gridLine = gridLine,
                  title = title,
                  titleSize = titleSize)
}
