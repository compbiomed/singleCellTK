#' visPlot
#'
#' Given a plotting method with condition and gene list, return
#' the respective visualization plot(s).
#'
#' @param inSCE Input \linkS4class{SingleCellExperiment} object.
#' @param useAssay The assay to use in the visualization plot. Required
#' @param method Visualization method. Available options are boxplot,
#' scatterplot, or heatmap. Required
#' @param condition colData annotation of the experiment. Required
#' @param glist selected genes for visualization. Maximum 25 genes. Required
#' @param facetWrap facet wrap according to genes for boxplot, scatterplot and
#' barplot. Default is FALSE. Optional
#' @param scaleHMap scale heatmap expression values. Default is TRUE. Optional
#' @param convertFactor If the condition is not a factor, convert it to a factor
#' before plotting. The default is FALSE
#'
#' @return A visualization plot
#'
#' @export
#' @examples
#' visPlot(mouseBrainSubsetSCE, "logcounts", "boxplot", "level1class", "C1qa")
#' visPlot(mouseBrainSubsetSCE, "counts", "scatterplot", "age", "Cmtm5")
#' visPlot(mouseBrainSubsetSCE, "counts", "heatmap", "level1class",
#'         c("Cmtm5", "C1qa"))
visPlot <- function(inSCE, useAssay, method, condition = NULL, glist,
                    facetWrap = TRUE, scaleHMap = TRUE, convertFactor = FALSE) {
  if (!inherits(inSCE, "SingleCellExperiment")){
    stop("Please use a SingleCellExperiment object")
  }
  #test for assay existing
  if (!all(useAssay %in% names(assays(inSCE)))){
    stop("assay '", useAssay, "' does not exist.")
  }
  #test for valid plot name
  if (!(method == "boxplot" | method == "scatterplot" | method == "heatmap" | method == "barplot")){
    stop("method '", method, "' is not a valid method.")
  }
  #test for gene list existing
  if (!all(glist %in% rownames(inSCE))){
    stop("Gene in gene list not found in input object.")
  }
  #test that only one condition is given
  if (length(condition) > 1){
    stop("Only 1 condition allowed")
  }

  #Main condition: check if the gene list is provided
  if (is.null(glist)){
    stop("gene list is required")
  } else {
    countsData <- data.frame(SummarizedExperiment::assay(inSCE, useAssay)[glist, , drop = FALSE])
    if (!is.null(condition)){
      if (inherits(inSCE, "SummarizedExperiment")) {
        annotData <- data.frame(SummarizedExperiment::colData(inSCE)[, condition, drop = FALSE])
      } else {
        annotData <- data.frame(SingleCellExperiment::colData(inSCE)[, condition, drop = FALSE])
      }
      if (convertFactor) {
        annotData[, condition] <- as.factor(annotData[, condition])
      }
      if (any(is.na(annotData[, condition]))){
        if (method == "heatmap" && is.factor(annotData[, condition])){
          #change NA to empty string for heatmap
          annotData[, condition] <- as.character(annotData[, condition])
          annotData[, condition][is.na(annotData[, condition])] <- ''
          annotData[, condition] <- factor(annotData[, condition])
        } else {
          stop("Annotation data has NA values. Filter them to continue.")
        }
      }
    } else{
      #condition required for boxplot or scatterplot
      if (!(method == "barplot" | method == "heatmap")){
        stop("Please supply a condition")
      }
    }
     #transforming the data into a data.frame
    if (!(method == "heatmap" | (method == "barplot" & is.null(condition)))){
      expDF <- cbind.data.frame(t(countsData), annotData)
      expDF$sample <- rownames(expDF)
      meltDF <- reshape2::melt(expDF, id.vars = c(condition, "sample"),
                               variable.name = "Genes", value.name = "assay")
    }
    #common for all except heatmap -- coloring each gene uniquely
    if (!method == 'heatmap'){
      if (length(glist) > 9){
        scale_values <- distinctColors(length(glist))
      } else {
        scale_values <- RColorBrewer::brewer.pal(9, "Set1")
      }
    }

    if (method == "boxplot"){
      if (length(glist) <= 25){
        if (is.factor(annotData[, condition])){
          ggplotObj <- ggplot2::ggplot(meltDF, ggplot2::aes_string(x = condition, y = "assay")) +
            ggplot2::geom_violin(ggplot2::aes_string(fill = condition)) +
            ggplot2::geom_boxplot(width = .1) +
            ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 90, hjust = 1)) +
            ggplot2::xlab(condition) +
            ggplot2::ylab(useAssay) +
            ggplot2::scale_fill_manual(values = scale_values,
                                       guide = FALSE)
          if (facetWrap) {
            ggplotObj + ggplot2::facet_wrap("Genes", scales = 'free')
          } else {
            ggplotObj
          }
        } else{
          stop("Boxplot requires condition to be a factor, use scatterplot instead")
        }
      } else{
        stop("Maximum limit of genes reached. Please enter 25 or less genes.")
      }
    } else if (method == "scatterplot"){
      if (length(glist) <= 25){
        if (!is.factor(annotData[, condition])){
          ggplotObj <- ggplot2::ggplot(meltDF, ggplot2::aes_string(x = condition, y = "assay")) +
            ggplot2::geom_point(ggplot2::aes_string(col = "Genes")) +
            ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 90, hjust = 1)) +
            ggplot2::xlab(condition) +
            ggplot2::ylab(useAssay) +
            ggplot2::scale_fill_manual(values = scale_values,
                                       guide = FALSE)
          if (facetWrap) {
            ggplotObj + ggplot2::facet_wrap("Genes", scales = 'free')
          } else {
            ggplotObj
          }
        } else {
          stop("Scatterplot requires a condition to be continuous, use boxplot as an alternative")
        }
      } else {
        stop("Maximum limit of genes reached. Please enter 25 or less genes.")
      }
    } else if (method == "barplot"){
      if (length(glist) <= 25){
        if (!is.null(condition)){
          ggplotObj <- ggplot2::ggplot(meltDF, ggplot2::aes_string(x = condition, y = "assay"))
          } else {
          scDF <- data.frame(t(countsData))
          scDF$sample <- rownames(scDF)
          meltDF <- reshape2::melt(scDF, id.vars = "sample", variable.name = "Genes",
                                   value.name = "assay")
          ggplotObj <- ggplot2::ggplot(meltDF, ggplot2::aes_string(x = 'sample', y = "assay"))
        }
        ggplotMain <- ggplotObj +
          ggplot2::geom_bar(stat = "identity", position = "dodge", ggplot2::aes_string(fill = "Genes")) +
          ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 90, hjust = 1)) +
          ggplot2::xlab("Sample") +
          ggplot2::ylab(useAssay) +
          ggplot2::scale_fill_manual(values = scale_values, guide = FALSE)
        if (facetWrap){
          ggplotMain + ggplot2::facet_grid("Genes", scales = 'free')
        } else {
          ggplotMain + ggplot2::scale_fill_manual(values = scale_values, guide = ggplot2::guide_legend())
        }
      } else{
        stop("Maximum limit of genes reached. Please enter 25 or less genes.")
      }
    } else if (method == "heatmap"){
      zeroSum <- which(matrixStats::rowSds(
        SummarizedExperiment::assay(inSCE, useAssay)[glist, , drop = FALSE]) == 0)
      if (length(zeroSum) != 0){
        stop("Gene ", paste(glist[zeroSum], collapse = ","), " has zero variance, please filter and continue.")
      }
      if (length(glist) > 0 & is.null(condition)){
        ComplexHeatmap::Heatmap(t(scale(t(countsData[glist, ]))),
                                name = "Expression",
                                column_title = "Differential Expression")
      } else{
        annData <- annotData[, condition]
        condLevels <- unique(annotData[, condition])
        if (length(condLevels) > 9){
          colors <- distinctColors(length(condLevels))
        } else {
          colors <- RColorBrewer::brewer.pal(9, "Set1")
        }
        if (is.factor(annData)){
          col <- list()
          col[[condition]] <- stats::setNames(colors[seq_along(condLevels)],
                                              condLevels)
          topha <- ComplexHeatmap::HeatmapAnnotation(
            df = annotData[, condition, drop = FALSE], col = col)
        } else if (is.numeric(annData)) {
          col <- list()
          col[[condition]] <- circlize::colorRamp2(c(min(annotData[, condition]),
                                                     max(annotData[, condition])),
                                                   c("white", "darkgreen"))
          topha <- ComplexHeatmap::HeatmapAnnotation(
            df = annotData[, condition, drop = FALSE],
            col = col)
        } else {
          stop("Data doesn't appear to be a factor or numeric type. Verify the",
               "annotation data.")
        }
        heatmapkey <- useAssay
        heatdata <- SummarizedExperiment::assay(inSCE, useAssay)[glist, ]
        if (scaleHMap){
          heatdata <- t(scale(t(heatdata)))
          heatmapkey <- paste("Scaled", heatmapkey, sep = "\n")
        }
        ComplexHeatmap::draw(
          ComplexHeatmap::Heatmap(
            heatdata,
            name = heatmapkey, column_title = "Differential Expression",
            cluster_rows = TRUE, cluster_columns = TRUE, top_annotation = topha,
            show_row_names = TRUE, show_column_names = TRUE,
            show_row_dend = TRUE,
            show_column_dend = TRUE
          ), heatmap_legend_side = "right", annotation_legend_side = "bottom"
        )
      }
    }
  }
}
