#' Calculate Differential Abundance with FET
#' @details This function will calculate the cell counting and fraction by
#' dividing all cells to groups specified by the arguments, together with
#' statistical summary by performing Fisher Exact Tests (FET).
#' @param inSCE A \code{\link[SingleCellExperiment]{SingleCellExperiment}}
#' object.
#' @param cluster A single \code{character}, specifying the name to store the
#' cluster label in \code{\link{colData}}.
#' @param variable A single \code{character}, specifying the name to store the
#' phenotype labels in \code{\link{colData}}.
#' @param control \code{character}. Specifying one or more categories that can
#' be found in the vector specified by \code{variable}.
#' @param case \code{character}. Specifying one or more categories that can
#' be found in the vector specified by \code{variable}.
#' @param analysisName A single \code{character}. Will be used for naming the
#' result table, which will be saved in metadata slot.
#' @return The original \code{\link[SingleCellExperiment]{SingleCellExperiment}}
#' object with \code{metadata(inSCE)} updated with a list
#' \code{diffAbundanceFET}, containing a new \code{data.frame} for the analysis
#' result, named by \code{analysisName}. The \code{data.frame} contains columns
#' for number and fraction of cells that belong to different cases, as well as
#' "Odds_Ratio", "PValue" and "FDR".
#' @export
#' @examples
#' data("mouseBrainSubsetSCE", package = "singleCellTK")
#' mouseBrainSubsetSCE <- diffAbundanceFET(inSCE = mouseBrainSubsetSCE,
#'                                                 cluster = "tissue",
#'                                                 variable = "level1class",
#'                                                 case = "oligodendrocytes",
#'                                                 control = "microglia",
#'                                                 analysisName = "diffAbundFET")
diffAbundanceFET <- function(inSCE, cluster, variable, control, case,
                             analysisName) {
  if (!inherits(inSCE, "SingleCellExperiment")) {
    stop("'inSCE' should be a SingleCellExperiment object.")
  }
  if (is.null(cluster) ||
      !cluster %in% names(SummarizedExperiment::colData(inSCE))) {
    stop("Argument 'cluster' should be a column name in colData(inSCE).")
  }
  if (is.null(variable) ||
      !variable %in% names(SummarizedExperiment::colData(inSCE))) {
    stop("Argument 'variable' should be a column name in colData(inSCE).")
  }

  cluster <- SummarizedExperiment::colData(inSCE)[, cluster]
  variable <- SummarizedExperiment::colData(inSCE)[, variable]
  if (!all(control %in% variable)) {
    nf1 <- control[which(!control %in% variable)]
    stop("Given 'control' not all found in 'variable': ",
         paste(nf1, collapse = ", "))
  }
  if (!all(case %in% variable)) {
    nf2 <- case[which(!case %in% variable)]
    stop("Given 'case' not all found in 'variable': ",
         paste(nf2, collapse = ", "))
  }
  if (any(control %in% case)) {
    warning("Overlapping variables found in 'control' and 'case'")
  }

  control.lab <- paste(control, collapse=",")
  case.lab <- paste(case, collapse=",")
  label <- rep(NA, length(variable))
  label[variable %in% case] <- case.lab
  label[variable %in% control] <- control.lab
  label <- factor(label, levels=c(control.lab, case.lab))

  res <- matrix(NA, nrow=length(unique(cluster)), ncol=10)
  cluster.label <- sort(unique(cluster))
  for(i in seq_along(cluster.label)) {
    cluster.factor <- factor(ifelse(cluster %in% cluster.label[i],
                                    cluster.label[i], "Other"),
                             levels=c(i, "Other"))

    ta <- table(cluster.factor, label)
    ta.p <- prop.table(ta, 2)
    fet <- stats::fisher.test(ta)
    res[i,] <- c(ta, ta.p, fet$estimate, fet$p.value)
  }
  colnames(res) <- c(paste0("Number of cells in cluster and in ",
                            control.lab),
                     paste0("Number of cells NOT in cluster and in ",
                            control.lab),
                     paste0("Number of cells in cluster and in ",
                            case.lab),
                     paste0("Number of cells NOT in cluster and in ",
                            case.lab),
                     paste0("Fraction of cells in cluster and in ",
                            control.lab),
                     paste0("Fraction of cells NOT in cluster and in ",
                            control.lab),
                     paste0("Fraction of cells in cluster and in ",
                            case.lab),
                     paste0("Fraction of cells NOT in cluster and in ",
                            case.lab),
                     "Odds_Ratio", "Pvalue")
  res <- data.frame(Cluster=cluster.label, res,
                    FDR=stats::p.adjust(res[,"Pvalue"], 'fdr'),
                    check.names=FALSE)
  getDiffAbundanceResults(inSCE, analysisName = analysisName) <- res
  return(inSCE)
}

#' Get/Set diffAbundanceFET result table
#' @rdname getDiffAbundanceResults
#' @param x A \code{\link[SingleCellExperiment]{SingleCellExperiment}}
#' object.
#' @param analysisName A single character string specifying an analysis 
#' performed with \code{\link{diffAbundanceFET}}
#' @param value The output table of \code{\link{diffAbundanceFET}}
#' @return The differential abundance table for getter method, or update the SCE
#' object with new result for setter method. 
#' @export
#' @examples 
#' data("mouseBrainSubsetSCE", package = "singleCellTK")
#' mouseBrainSubsetSCE <- diffAbundanceFET(inSCE = mouseBrainSubsetSCE,
#'                                                 cluster = "tissue",
#'                                                 variable = "level1class",
#'                                                 case = "oligodendrocytes",
#'                                                 control = "microglia",
#'                                                 analysisName = "diffAbund")
#' result <- getDiffAbundanceResults(mouseBrainSubsetSCE, "diffAbund")
setGeneric("getDiffAbundanceResults", signature = "x",
           function(x, analysisName) {
             standardGeneric("getDiffAbundanceResults")
           }
)

#' @rdname getDiffAbundanceResults
#' @export
setMethod("getDiffAbundanceResults", signature(x = "SingleCellExperiment"), 
          function(x, analysisName) {
  result.names <- names(S4Vectors::metadata(x)$sctk$diffAbundanceFET)
  if(!analysisName %in% result.names) {
    stop("The analysis '", analysisName, "' was not found in the results ", 
         "for tool 'diffAbundanceFET'")
  }
  
  results <- S4Vectors::metadata(x)$sctk$diffAbundanceFET[[analysisName]]
  return(results)
})

#' @rdname getDiffAbundanceResults
#' @export
setGeneric("getDiffAbundanceResults<-", function(x, analysisName, value) {
  standardGeneric("getDiffAbundanceResults<-")
  })

#' @rdname getDiffAbundanceResults
#' @export
setReplaceMethod("getDiffAbundanceResults", 
                 signature(x = "SingleCellExperiment"), 
                 function(x, analysisName, value) {
  S4Vectors::metadata(x)$sctk$diffAbundanceFET[[analysisName]] <- value
  return(x)
})

#' Plot the differential Abundance
#' @details This function will visualize the differential abundance in two given
#' variables, by making bar plots that presents the cell counting and fraction
#' in different cases.
#' @param inSCE A \code{\link[SingleCellExperiment]{SingleCellExperiment}}
#' object.
#' @param cluster A single \code{character}, specifying the name to store the
#' cluster label in \code{\link{colData}}.
#' @param variable A single \code{character}, specifying the name to store the
#' phenotype labels in \code{\link{colData}}.
#' @param combinePlot Must be either "all" or "none". "all" will combine all 
#' plots into a single \code{\link[ggplot2]{ggplot}} object. Default 
#' \code{"all"}.
#' @return When \code{combinePlot = "none"}, a \code{list} with 4 
#' \code{\link[ggplot2]{ggplot}} objects; when \code{combinePlot = "all"}, a 
#' single \code{\link[ggplot2]{ggplot}} object with for subplots. 
#' @export
#' @examples
#' data("mouseBrainSubsetSCE", package = "singleCellTK")
#' plotClusterAbundance(inSCE = mouseBrainSubsetSCE,
#'                      cluster = "tissue",
#'                      variable = "level1class")
plotClusterAbundance <- function(inSCE, cluster, variable, 
                                 combinePlot = c("all", "none")) {
  if (!inherits(inSCE, "SingleCellExperiment")) {
    stop("'inSCE' should be a SingleCellExperiment object.")
  }
  if (is.null(cluster) ||
      !cluster %in% names(SummarizedExperiment::colData(inSCE))) {
    stop("Argument 'cluster' should be a column name in colData(inSCE).")
  }
  if (is.null(variable) ||
      !variable %in% names(SummarizedExperiment::colData(inSCE))) {
    stop("Argument 'variable' should be a column name in colData(inSCE).")
  }
  combinePlot <- match.arg(combinePlot)
  cluster <- SummarizedExperiment::colData(inSCE)[, cluster]
  color_palette <- distinctColors(length(unique(cluster)))

  label <- SummarizedExperiment::colData(inSCE)[, variable]

  cluster.color <- color_palette
    df <- data.frame(Cluster=as.factor(cluster), Sample=as.factor(label))

  g1 <- ggplot2::ggplot(df, ggplot2::aes_string("Cluster")) +
    ggplot2::geom_bar(ggplot2::aes_string(fill="Sample")) +
    ggplot2::theme_classic() +
    ggplot2::scale_y_continuous(expand = c(0, 0))

  g2 <- ggplot2::ggplot(df, ggplot2::aes_string("Cluster")) +
    ggplot2::geom_bar(ggplot2::aes_string(fill="Sample"), position="fill") +
    ggplot2::theme_classic() +
    ggplot2::scale_y_continuous(expand = c(0, 0)) +
    ggplot2:: ylab("Fraction")

  g3 <- ggplot2::ggplot(df, ggplot2::aes_string("Sample")) +
    ggplot2::geom_bar(ggplot2::aes_string(fill="Cluster")) +
    ggplot2::theme_classic() +
    ggplot2::scale_y_continuous(expand = c(0, 0)) +
    ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45, hjust = 1)) +
    ggplot2::scale_fill_manual(values=cluster.color)

  g4 <- ggplot2::ggplot(df, ggplot2::aes_string("Sample")) +
    ggplot2::geom_bar(ggplot2::aes_string(fill="Cluster"), position="fill") +
    ggplot2::theme_classic() +
    ggplot2::scale_y_continuous(expand = c(0, 0)) +
    ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45, hjust = 1)) +
    ggplot2::scale_fill_manual(values=cluster.color) +
    ggplot2::ylab("Fraction")
  g <- list(g1, g2, g3, g4)
  if (combinePlot == "all") {
    g <- cowplot::plot_grid(g1, g3, g2, g4, ncol = 2, 
                            align = "hv", axis = "blr")
  } else if (combinePlot == "none") {
    g <- list(g1, g2, g3, g4)
  }
  return(g)
}
