.myenv <- new.env(parent = emptyenv())

#' GSVA_sce
#'
#' Run GSVA analysis on a SCESet object.
#'
#' @param SCEdata SCESet object
#' @param use_assay The assay to use for the MAST calculations. The default is
#' "logcounts"
#' @param pathway_source The pathway source if "Manual Input", the pathway_names
#' should be rowData annotations that are (0,1) vectors. If, "MSigDB c2 (Human,
#' Entrez ID only)",
#' the pathway_names should be pathways from MSigDB c2 or "ALL" to run on all
#' available pathways.
#' @param pathway_names List of pathway names to run, depending on
#' pathway_source
#' parameter.
#' @param ... Parameters to pass to gsva
#' 
#' @return A data.frame of pathway activity scores from GSVA.
#' 
#' @export
#'
GSVA_sce <- function(SCEdata, use_assay = "logcounts", pathway_source, pathway_names, ...){
  if (pathway_source == "Manual Input"){
    #expecting logical vector
    biomarker <- lapply(pathway_names, function(x) rownames(SCEdata)[SingleCellExperiment::rowData(SCEdata)[, x] == 1])
    gsva_res <- GSVA::gsva(SummarizedExperiment::assay(SCEdata, use_assay), biomarker, ...)
    rownames(gsva_res) <- pathway_names
  } else if (pathway_source == "MSigDB c2 (Human, Entrez ID only)") {
    utils::data("c2BroadSets", package = "GSVAdata", envir = .myenv)
    c2BroadSets <- .myenv$c2BroadSets
    #expecting some genes in list are in the rownames
    if ("ALL" %in% pathway_names) {
      gsva_res <- GSVA::gsva(SummarizedExperiment::assay(SCEdata, use_assay), c2BroadSets, ...)
    } else {
      c2sub <- c2BroadSets[base::setdiff(pathway_names, "ALL")]
      gsva_res <- GSVA::gsva(SummarizedExperiment::assay(SCEdata, use_assay), c2sub, ...)
    }
  } else{
    stop("ERROR: Unsupported gene list source ", pathway_source)
  }
  return(gsva_res)
}

#' GSVA_plot
#'
#' Plot GSVA results.
#'
#' @param SCEdata SCESet object
#' @param gsva_data GSVA data to plot. Required.
#' @param plot_type The type of plot to use, "Violin" or "Heatmap". Required.
#' @param condition The condition(s) to use for the Violin plot, or the condition(s)
#' to add as color bars above the Heatmap. Required for Violin, optional for Heatmap.
#' 
#' @return The requested plot of the GSVA results.
#' 
#' @export
#'
GSVA_plot <- function(SCEdata, gsva_data, plot_type, condition=NULL){
  if (plot_type == "Violin"){
    if (nrow(gsva_data) > 49){
      stop("TOO Many results for Violin Plot. Try Heatmap.")
    }
    if (!(length(condition) > 0)){
      stop("You must specify a condition for Violin plot")
    }
    gsva_res_t <- data.frame(t(gsva_data))
    cond <- apply(SingleCellExperiment::colData(SCEdata)[, condition, drop = FALSE], 1, paste, collapse = "_")
    gsva_res_t[, paste(condition, collapse = "_")] <- cond
    gsva_res_flat <- reshape2::melt(gsva_res_t, id.vars = paste(condition, collapse = "_"),
                                    variable.name = "pathway")
    ggbase <- ggplot2::ggplot(gsva_res_flat, ggplot2::aes_string(x = paste(condition, collapse = "_"),
                                                                 y = "value",
                                                                 color = paste(condition, collapse = "_"))) +
      ggplot2::geom_violin() +
      ggplot2::geom_jitter() +
      ggplot2::facet_wrap(~pathway, scale = "free_y", ncol = ceiling(sqrt(nrow(gsva_data))))
    return(ggbase)
  } else if (plot_type == "Heatmap"){
    topha <- NULL
    if (length(condition) > 0){
      colors <- RColorBrewer::brewer.pal(8, "Set1")
      cond <- apply(SingleCellExperiment::colData(SCEdata)[, condition, drop = FALSE], 1, paste, collapse = "_")
      cond_levels <- unique(cond)
      if (length(cond_levels) < 8){
        col <- list()
        col[[paste(condition, collapse = "_")]] <- stats::setNames(colors[1:length(cond_levels)], cond_levels)
        conddf <- data.frame(cond, row.names = colnames(gsva_data))
        colnames(conddf) <- paste(condition, collapse = "_")
        topha <- ComplexHeatmap::HeatmapAnnotation(df = conddf,
                                                   col = col)
      } else {
        stop("Too many levels in selected condition(s)")
      }
    }
    ComplexHeatmap::Heatmap(gsva_data, top_annotation = topha)
  } else {
    stop("ERROR: Unsupported plot type")
  }
}
