.myenv <- new.env(parent = emptyenv())

#' Run GSVA analysis on a \linkS4class{SingleCellExperiment} object
#'
#' @param inSCE Input \linkS4class{SingleCellExperiment} object.
#' @param useAssay Indicate which assay to use. The default is "logcounts"
#' @param pathwaySource The pathway source if "Manual Input", the pathwayNames
#' should be rowData annotations that are (0,1) vectors. If, "MSigDB c2 (Human,
#' Entrez ID only)", the pathwayNames should be pathways from MSigDB c2 or "ALL"
#' to run on all available pathways.
#' @param pathwayNames List of pathway names to run, depending on pathwaySource
#' parameter.
#' @param ... Parameters to pass to gsva()
#'
#' @return gsvaSCE(): A data.frame of pathway activity scores from GSVA.
#' @export
gsvaSCE <- function(inSCE, useAssay = "logcounts", pathwaySource,
                    pathwayNames, ...){
  if (pathwaySource == "Manual Input"){
    #expecting logical vector
    biomarker <- lapply(pathwayNames, function(x) rownames(inSCE)[
      SingleCellExperiment::rowData(inSCE)[, x] == 1])
    gsvaRes <- GSVA::gsva(SummarizedExperiment::assay(inSCE, useAssay),
                          biomarker, ...)
    rownames(gsvaRes) <- pathwayNames
  } else if (pathwaySource == "MSigDB c2 (Human, Entrez ID only)") {
    utils::data("c2BroadSets", package = "GSVAdata", envir = .myenv)
    c2BroadSets <- .myenv$c2BroadSets
    #expecting some genes in list are in the rownames
    if ("ALL" %in% pathwayNames) {
      gsvaRes <- GSVA::gsva(SummarizedExperiment::assay(inSCE, useAssay),
                            c2BroadSets, ...)
    } else {
      c2sub <- c2BroadSets[base::setdiff(pathwayNames, "ALL")]
      gsvaRes <- GSVA::gsva(SummarizedExperiment::assay(inSCE, useAssay),
                            c2sub, ...)
    }
  } else{
    stop("ERROR: Unsupported gene list source ", pathwaySource)
  }
  return(gsvaRes)
}

#' @describeIn gsvaSCE Plot GSVA results.
#'
#' Plot GSVA Results
#'
#' @param gsvaData GSVA data to plot. Required.
#' @param plotType The type of plot to use, "Violin" or "Heatmap". Required.
#' @param condition The condition(s) to use for the Violin plot, or the
#' condition(s) to add as color bars above the Heatmap. Required for Violin,
#' optional for Heatmap.
#' @param show_row_names Display the row labels on the heatmap. The default
#' is TRUE.
#' @param show_column_names Display the column labels on the heatmap. The
#' default is TRUE
#' @param text_size Text size for plots. The default is 12
#'
#' @return gsvaPlot(): The requested plot of the GSVA results.
#'
#' @export
gsvaPlot <- function(inSCE, gsvaData, plotType, condition=NULL,
                     show_column_names = TRUE, show_row_names = TRUE,
                     text_size = 12){
  if (plotType == "Violin"){
    if (nrow(gsvaData) > 49){
      stop("Too Many results for Violin Plot. Try Heatmap.")
    }
    if (!(length(condition) > 0)){
      stop("You must specify a condition for Violin plot")
    }
    gsvaResT <- data.frame(t(gsvaData))
    cond <- apply(SingleCellExperiment::colData(inSCE)[, condition,
                                                         drop = FALSE],
                  1, paste, collapse = "_")
    gsvaResT[, paste(condition, collapse = "_")] <- cond
    gsvaResFlat <- reshape2::melt(gsvaResT, id.vars = paste(condition,
                                                            collapse = "_"),
                                    variable.name = "pathway")
    ggbase <- ggplot2::ggplot(gsvaResFlat, ggplot2::aes_string(
      x = paste(condition, collapse = "_"), y = "value",
      color = paste(condition, collapse = "_"))) +
      ggplot2::geom_violin() +
      ggplot2::geom_jitter() +
      ggplot2::facet_wrap(~pathway, scale = "free_y",
                          ncol = ceiling(sqrt(nrow(gsvaData)))) +
      ggplot2::theme(strip.text.x = ggplot2::element_text(size = text_size))
    return(ggbase)
  } else if (plotType == "Heatmap"){
    topha <- NULL
    if (length(condition) > 0){
      cond <- apply(SingleCellExperiment::colData(inSCE)[, condition,
                                                           drop = FALSE],
                    1, paste, collapse = "_")
      condLevels <- unique(cond)
      if (length(condLevels) > 9){
        colors <- distinctColors(length(condLevels))
      } else {
        colors <- RColorBrewer::brewer.pal(9, "Set1")
      }
      col <- list()
      col[[paste(condition, collapse = "_")]] <-
        stats::setNames(colors[seq_along(condLevels)], condLevels)
      conddf <- data.frame(cond, row.names = colnames(gsvaData))
      colnames(conddf) <- paste(condition, collapse = "_")
      topha <- ComplexHeatmap::HeatmapAnnotation(df = conddf,
                                                 col = col)
    }
    #TODO make more customizable
    return(ComplexHeatmap::draw(
      ComplexHeatmap::Heatmap(gsvaData,
                              row_names_gp = grid::gpar(fontsize = text_size),
                              top_annotation = topha,
                              show_column_names = show_column_names,
                              show_row_names = show_row_names,
                              name = "GSVA\nscore"),
      heatmap_legend_side = "left", annotation_legend_side = "bottom"
    ))
  } else {
    stop("ERROR: Unsupported plot type")
  }
}
