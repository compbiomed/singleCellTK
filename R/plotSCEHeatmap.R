#' Plot heatmap of using data stored in SingleCellExperiment Object
#' @rdname plotSCEHeatmap
#' @param inSCE \linkS4class{SingleCellExperiment} inherited object.
#' @param useAssay character. A string indicating the assay name that
#' provides the expression level to plot. Only for \code{plotSCEHeatmap}.
#' @param useReducedDim character. A string indicating the reducedDim name that
#' provides the expression level to plot. Only for \code{plotSCEDimReduceHeatmap}.
#' @param doLog Logical scalar. Whether to do \code{log(assay + 1)}
#' transformation on the assay indicated by \code{useAssay}. Default
#' \code{FALSE}.
#' @param featureIndex A vector that can subset the input SCE object by rows
#' (features). Alternatively, it can be a vector identifying features in
#' another feature list indicated by \code{featureIndexBy}. Default \code{NULL}.
#' @param cellIndex A vector that can subset the input SCE object by columns
#' (cells). Alternatively, it can be a vector identifying cells in another
#' cell list indicated by \code{featureIndexBy}. Default \code{NULL}.
#' @param scale Whether to perform z-score or min-max scaling on each row.Choose from \code{"zscore"}, \code{"min-max"} or default
#' \code{TRUE} or \code{FALSE}
#' @param trim A 2-element numeric vector. Values outside of this range will be
#' trimmed to their nearst bound. Default \code{c(-2, 2)}
#' @param featureIndexBy A single character specifying a column name of
#' \code{rowData(inSCE)}, or a vector of the same length as \code{nrow(inSCE)},
#' where we search for the non-rowname feature indices. Not applicable for
#' \code{plotSCEDimReduceHeatmap}. Default \code{"rownames"}.
#' @param cellIndexBy A single character specifying a column name of
#' \code{colData(inSCE)}, or a vector of the same length as \code{ncol(inSCE)},
#' where we search for the non-rowname cell indices. Default \code{"rownames"}.
#' @param cluster_columns A logical scalar that turns on/off 
#' clustering of columns. Default \code{FALSE}. Clustering columns should be turned off when using reduced dim 
#' for plotting as it will be sorted by PCs
#' @param cluster_rows A logical scalar that turns on/off clustering of rows. 
#' Default \code{FALSE}.
#' @param rowDataName character. The column name(s) in \code{rowData} that need
#' to be added to the annotation. Not applicable for
#' \code{plotSCEDimReduceHeatmap}. Default \code{NULL}.
#' @param colDataName character. The column name(s) in \code{colData} that need
#' to be added to the annotation. Default \code{NULL}.
#' @param aggregateRow Feature variable for aggregating the heatmap by row. Can
#' be a vector or a \code{rowData} column name for feature variable. Multiple
#' variables are allowed. Default \code{NULL}.
#' @param aggregateCol Cell variable for aggregating the heatmap by column. Can
#' be a vector or a \code{colData} column name for cell variable. Multiple
#' variables are allowed. Default \code{NULL}.
#' @param featureAnnotations \code{data.frame}, with \code{rownames} containing
#' all the features going to be plotted. Character columns should be factors.
#' Default \code{NULL}.
#' @param cellAnnotations \code{data.frame}, with \code{rownames} containing
#' all the cells going to be plotted. Character columns should be factors.
#' Default \code{NULL}.
#' @param featureAnnotationColor A named list. Customized color settings for
#' feature labeling. Should match the entries in the \code{featureAnnotations}
#' or \code{rowDataName}. For each entry, there should be a list/vector of
#' colors named with categories. Default \code{NULL}.
#' @param cellAnnotationColor A named list. Customized color settings for
#' cell labeling. Should match the entries in the \code{cellAnnotations} or
#' \code{colDataName}. For each entry, there should be a list/vector of colors
#' named with categories. Default \code{NULL}.
#' @param palette Choose from \code{"ggplot"}, \code{"celda"} or \code{"random"}
#' to generate unique category colors.
#' @param heatmapPalette Choose from \code{"sequential"}, \code{"diverging"} or supply custom palette with colorScheme
#' to generate unique category colors. Default is \code{"sequential"}
#' @param addCellSummary Add summary barplots to column annotation. Supply the name of the column in colData as a character. This option will add summary for categorical variables 
#' as stacked barplots.
#' @param rowSplitBy character. Do semi-heatmap based on the grouping of
#' this(these) annotation(s). Should exist in either \code{rowDataName} or
#' \code{names(featureAnnotations)}. Default \code{NULL}.
#' @param colSplitBy character. Do semi-heatmap based on the grouping of
#' this(these) annotation(s). Should exist in either \code{colDataName} or
#' \code{names(cellAnnotations)}. Default \code{NULL}.
#' @param rowLabel Use a logical for whether to display all the feature names,
#' a single character to display a column of \code{rowData(inSCE)} annotation,
#' a vector of the same length as full/subset \code{nrow(inSCE)} to display
#' customized info. Default \code{FALSE}.
#' @param colLabel Use a logical for whether to display all the cell names, a
#' single character to display a column of \code{colData(inSCE)} annotation,
#' a vector of the same length as full/subset \code{ncol(inSCE)} to display
#' customized info. Default \code{FALSE}.
#' @param rowLabelSize A number for the font size of feature names. Default
#' \code{8}
#' @param colLabelSize A number for the font size of cell names. Default
#' \code{8}
#' @param rowDend Whether to display row dendrogram. Default \code{TRUE}.
#' @param colDend Whether to display column dendrogram. Default \code{TRUE}.
#' @param title The main title of the whole plot. Default \code{NULL}.
#' @param rowTitle The subtitle for the rows. Default \code{"Genes"}.
#' @param colTitle The subtitle for the columns. Default \code{"Cells"}.
#' @param rowGap A numeric value or a \code{\link[grid]{unit}} object. For the
#' gap size between rows of the splitted heatmap. Default
#' \code{grid::unit(0, 'mm')}.
#' @param colGap A numeric value or a \code{\link[grid]{unit}} object. For the
#' gap size between columns of the splitted heatmap. Default
#' \code{grid::unit(0, 'mm')}.
#' @param border A logical scalar. Whether to show the border of the heatmap or
#' splitted heatmaps. Default \code{TRUE}.
#' @param colorScheme function. A function that generates color code by giving
#' a value. Can be generated by \code{\link[circlize]{colorRamp2}}.
#' Default \code{NULL}.
#' @param ... Other arguments passed to \code{\link[ComplexHeatmap]{Heatmap}}.
#' @examples
#' data(scExample, package = "singleCellTK")
#' plotSCEHeatmap(sce[1:3,1:3], useAssay = "counts")
#' @return A \code{\link[ggplot2]{ggplot}} object.
#' @export
#' @author Yichen Wang
#' @importFrom scuttle aggregateAcrossCells aggregateAcrossFeatures
#' @importFrom SingleCellExperiment SingleCellExperiment
#' @importFrom SummarizedExperiment colData assayNames<-
#' @importFrom stringr str_replace_all str_c
#' @importFrom stats prcomp quantile
#' @importFrom dplyr select arrange group_by count ungroup mutate one_of desc
#' @importFrom tidyr spread unite 
#' @importFrom tibble column_to_rownames remove_rownames
#' @importFrom grid gpar
#' @importFrom ComplexHeatmap anno_barplot
#' @importFrom rlang .data
#' @importFrom S4Vectors DataFrame
#' 
#' 
plotSCEHeatmap <- function(inSCE, useAssay = 'logcounts', useReducedDim = NULL,
                           doLog = FALSE, featureIndex = NULL, cellIndex = NULL,
                           scale = TRUE, trim = c(-2,2),
                           featureIndexBy = 'rownames',
                           cellIndexBy = 'rownames',
                           cluster_columns = FALSE,
                           cluster_rows = FALSE,
                           rowDataName = NULL, colDataName = NULL,
                           aggregateRow = NULL, aggregateCol = NULL,
                           featureAnnotations = NULL, cellAnnotations = NULL,
                           featureAnnotationColor = NULL,
                           cellAnnotationColor = NULL,
                           palette = c("ggplot", "celda", "random"),
                           heatmapPalette = c("sequential","diverging"),
                           addCellSummary = NULL,
                           rowSplitBy = NULL, colSplitBy = NULL,
                           rowLabel = FALSE, colLabel = FALSE,
                           rowLabelSize = 6, colLabelSize = 6,
                           rowDend = TRUE, colDend = TRUE,
                           title = NULL, rowTitle = 'Features',
                           colTitle = 'Cells',
                           rowGap = grid::unit(0, 'mm'),
                           colGap = grid::unit(0, 'mm'),
                           border = FALSE, colorScheme = NULL, ...){
  palette<-match.arg(palette)
  heatmapPalette<-match.arg(heatmapPalette)
  # STAGE 1: Create clean SCE object with only needed information ####
  ## .selectSCEMatrix, .manageCellVar and .manageFeatureVar perform checks
  useMat <- .selectSCEMatrix(inSCE, useAssay = useAssay,
                             useReducedDim = useReducedDim,
                             returnMatrix = TRUE, cellAsCol = TRUE)
  useAssay <- useMat$names$useAssay
  useReducedDim <- useMat$names$useReducedDim
  useData <- ifelse(!is.null(useAssay), useAssay, useReducedDim)
  ### cell annotation
  oldColData <- colData(inSCE)
  colDataName <- unique(c(colDataName, aggregateCol))
  colDataAnns <- lapply(colDataName, function(x) .manageCellVar(inSCE, x))
  if (length(colDataName) > 0)
    colDataAnns <- data.frame(colDataAnns, row.names = colnames(inSCE))
  else
    colDataAnns <- data.frame(row.names = colnames(inSCE))
  colnames(colDataAnns) <- colDataName
  cellAnnotations <- .mergeAnnotationDF(colDataAnns, cellAnnotations)
  if (!is.null(colSplitBy) &&
      any(!colSplitBy %in% colnames(cellAnnotations)))
    stop('Specified `colSplitBy` variables not found.')
  if (isTRUE(colLabel)) {
    colLabelName <- colnames(inSCE)
  } else if (isFALSE(colLabel)) {
    colLabelName <- NULL
  } else {
    colLabelName <- .manageCellVar(inSCE, colLabel)
    colLabel <- TRUE
  }
  ### feature annotation
  rowDataAnns <- data.frame(row.names = rownames(useMat$mat))
  if (!is.null(useAssay)) {
    # When using reducedDim, no rowData can be applied
    rowDataName <- unique(c(rowDataName, aggregateRow))
    rowDataAnns <- lapply(rowDataName, function(x) .manageFeatureVar(inSCE, x))
    if (length(rowDataName) > 0)
      rowDataAnns <- data.frame(rowDataAnns, row.names = rownames(inSCE))
    else
      rowDataAnns <- data.frame(row.names = rownames(inSCE))
    colnames(rowDataAnns) <- rowDataName
  }
  # But customized featureAnnotations should work
  featureAnnotations <- .mergeAnnotationDF(rowDataAnns, featureAnnotations)
  if (!is.null(rowSplitBy) &&
      any(!rowSplitBy %in% colnames(featureAnnotations)))
    stop('Specified `rowSplitBy` variables not found.')
  if (isTRUE(rowLabel)) {
    rowLabelName <- rownames(useMat$mat)
  } else if (isFALSE(rowLabel)) {
    rowLabelName <- NULL
  } else {
    if (!is.null(useAssay)) {
      rowLabelName <- .manageFeatureVar(inSCE, rowLabel)
      rowLabel <- TRUE
    } else {
      # Using customized rowLabel for reducedDim
      if (length(rowLabel) != nrow(useMat$mat))
        stop("Length of `rowLabel` does not match nrow of specified ",
             "`useReducedDim`")
      rowLabelName <- rowLabel
      rowLabel <- TRUE
    }
  }
  ### create SCE object
  SCE <- SingleCellExperiment(assay = list(useMat$mat),
                              colData = cellAnnotations,
                              rowData = featureAnnotations)
  assayNames(SCE) <- useData
  
  .minmax<-function(mat){
    min_max<- function(x) {
      new_x =  (x - min(x))/ (max(x) - min(x))
      return(new_x)}
    new_mat<-as.matrix(apply(mat,FUN = min_max,MARGIN = 2))
    return(new_mat)
    }
  
  # STAGE 2: Subset SCE object as needed ####
  # Manage cell subsetting
  if(is.null(cellIndex)){
    cellIndex <- seq(ncol(SCE))
  } else if (is.character(cellIndex)) {
    # cellIndexBy not necessarily included in new "SCE"
    cellIndex <- retrieveSCEIndex(inSCE, cellIndex, axis = "col",
                                  by = cellIndexBy)
  } else if (is.logical(cellIndex)) {
    if (length(cellIndex) != ncol(inSCE)) {
      stop("Logical index length does not match ncol(inSCE)")
    }
    cellIndex <- which(cellIndex)
  }
  # Manage feature subsetting
  if(is.null(featureIndex)){
    featureIndex <- seq(nrow(SCE))
  } else if (is.character(featureIndex)) {
    if (!is.null(useAssay))
      featureIndex <- retrieveSCEIndex(inSCE, featureIndex, axis = "row",
                                       by = featureIndexBy)
    else
      # When using reducedDim, can only go with "PC" names
      # or customized "by"
      featureIndex <- retrieveSCEIndex(SCE, featureIndex, axis = "row",
                                       by = featureIndexBy)
  } else if (is.logical(featureIndex)) {
    if (length(featureIndex) != nrow(SCE)) {
      stop("Logical index length does not match nrow(inSCE)")
    }
    featureIndex <- which(featureIndex)
  }
  if(is.null(colLabelName)){
    colnames(SCE) <- NULL
  }
  else{
    colnames(SCE) <- colLabelName
  }
 
  if(is.null(rowLabelName)){
    rownames(SCE) <- NULL
  }
  else{
    rownames(SCE) <- rowLabelName
  }
  
  SCE <- SCE[featureIndex, cellIndex]
  ### Scaling should be done before aggregating
  if (isTRUE(doLog)) assay(SCE) <- log1p(assay(SCE))
  if(isTRUE(scale)) scale <- "zscore"
  if ((scale == "zscore")) {
    assay(SCE) <- as.matrix(base::scale(assay(SCE)))
  } else if (scale ==  "min_max") {
    assay(SCE) <- as.matrix(.minmax(assay(SCE)))
  }    
  
  
  # STAGE 3: Aggregate As needed ####
  if (!is.null(aggregateCol)) {
    # TODO: whether to also aggregate numeric variable that users want
    # Might need to use "coldata.merge" in aggregate function
    colIDS <- colData(SCE)[, aggregateCol]
    origRowData <- rowData(SCE)
    SCE <- aggregateAcrossCells(SCE, ids = colIDS,
                                use.assay.type = useData,
                                store.number = NULL, statistics = "mean")
    # TODO: `aggregateAcrossCells` produce duplicated variables in colData
    # and unwanted "ncell" variable even if I set `store.number = NULL`.
   #colData(SCE) <- colData(SCE)[,c(aggregateCol),drop=FALSE] ##change
    
    temp_df<-as.data.frame(colData(SCE)[,c(aggregateCol),drop=FALSE]) %>% 
      unite("new_colnames",dplyr::everything(),sep = "_",remove = FALSE) %>% 
      remove_rownames() %>% 
    #  mutate(aggregated_column = new_colnames) %>%
    #  dplyr::select(new_colnames, aggregated_column) %>%
      column_to_rownames("new_colnames")

    colData(SCE)<-DataFrame(temp_df)
    rowData(SCE) <- origRowData
  }
  if (!is.null(aggregateRow)) {
    # `aggregateAcrossFeatures` doesn't work by with multi-var
    # Remake one single variable vector
    rowIDS <- rowData(SCE)[, aggregateRow, drop = FALSE]
    rowIDS <- do.call(paste, c(rowIDS, list(sep = "_")))
    origColData <- colData(SCE)
    SCE <- aggregateAcrossFeatures(SCE, ids = rowIDS, average = TRUE,
                                   use.assay.type = useData)
    colData(SCE) <- origColData
  }
  # STAGE 4: Other minor preparation for plotting ####
 
  # Create a function that sorts the matrix by PC1
  .orderMatrix<-function(mat){
    # Adding extra character to rownames because presence of some char gets a "." if I don't
    mat2<-data.frame(t(mat))
    rownames(mat2)<-stringr::str_c("K_",rownames(mat2))
    pca_mat<-stats::prcomp(mat2,center = TRUE, scale. = FALSE)
    kl<-dplyr::arrange(data.frame(pca_mat$x)["PC1"],desc(data.frame(pca_mat$x)["PC1"]))
    mat<-data.frame(t(mat2)) %>% dplyr::select(rownames(kl))
    colnames(mat)<-stringr::str_replace_all(colnames(mat),"K_","")
    return(as.matrix(mat))
  }
  
  # Prepare
  
  if(!is.null(useReducedDim)){
    mat <- assay(SCE)
    mat <- .orderMatrix(mat)
    
  } else{
    
    if(class(assay(SCE))[1] == "dgCMatrix"){
      mat<- as.matrix(assay(SCE))
    }
    else{
     mat <- assay(SCE)
    }
  }
   
  
  if (!is.null(trim) & scale == "zscore") {
    assay(SCE) <- trimCounts(assay(SCE), trim)  
  }    
  
  
  if (is.null(colorScheme)) {
    if (isFALSE(scale)){
      if (heatmapPalette == "sequential"){
        colorScheme <- circlize::colorRamp2(quantile(mat,na.rm=TRUE),
                                            c('white', "#fecc5c",'#fdae61',"#f03b20","#bd0026"))
      }
      else if (heatmapPalette == "diverging"){
      colorScheme <- circlize::colorRamp2(c(min(mat),
                                            (max(mat) + min(mat))/2,
                                            max(mat)),
                                          c('blue', 'white', 'red'))
      }
    }
    else if (scale == "zscore"){
      colorScheme <- circlize::colorRamp2(quantile(assay(SCE), na.rm = TRUE),
                                          c('#2c7bb6','#abd9e9','#ffffbf','#fdae61','#d7191c'))
    }
    else if (scale == "min_max"){
      if(heatmapPalette == "sequential"){
        colorScheme <- circlize::colorRamp2(c(0,0.3,0.6,0.8,1),
                                            c('white', "#fecc5c",'#fdae61',"#f03b20","#bd0026")) 
        
      }
      else if (heatmapPalette == "diverging") {
      colorScheme <- circlize::colorRamp2(c(0,0.3,0.6,0.8,1),
                                          c('#2c7bb6','#abd9e9','#ffffbf','#fdae61','#d7191c'))     
      }
    }
  } else {
    if (!is.function(colorScheme))
      stop('`colorScheme` must be a function generated by ',
           'circlize::colorRamp2')
    breaks <- attr(colorScheme, 'breaks')
    if (breaks[1] != min(trim) || breaks[length(breaks)] != max(trim))
      stop('Breaks of `colorScheme` do not match with `trim`.')
  }
  
  # Avoid documentation error by setting these values to NULL
  # removed NOTE about namespaces.
  n<-value<-NULL
  
  ### Generate HeatmapAnnotation object
  ca <- NULL
  cellAnnotationColor <- .heatmapAnnColor(SCE, slot = "colData",
                                          custom = cellAnnotationColor,
                                          palette = palette)
  if(dim(cellAnnotations)[2] > 0)
    if(is.null(addCellSummary)){
      ca <- ComplexHeatmap::HeatmapAnnotation(df = as.data.frame(colData(SCE)),
                                              col = cellAnnotationColor)
    }
  else if (!addCellSummary %in% colnames(oldColData)){
    stop(addCellSummary,
         "' not found in colData")
  }
  else if (addCellSummary %in% colnames(oldColData)){
    oldColData %>%
      as.data.frame() %>%
      group_by(!!!rlang::syms(aggregateCol),!!!rlang::syms(addCellSummary)) %>%
      count() %>%
      ungroup() %>%
      group_by(!!! rlang::syms(aggregateCol)) %>%
      mutate(sum = sum(n)) %>%
      mutate(value = n/sum) %>%
      dplyr::select(-n,sum) %>%
      spread(one_of(addCellSummary),value) %>%
      ungroup() %>%
      dplyr::select(-one_of(aggregateCol),-sum) -> boxdata
    
    
    boxdata[is.na(boxdata)]  <- 0
    boxdata<-as.matrix(boxdata)
    ca <- ComplexHeatmap::HeatmapAnnotation(addCellSummary = anno_barplot(boxdata,
                                                                         gp = gpar(fill = 2:5)),
                                            annotation_label = addCellSummary,
                                            col = cellAnnotationColor)
  }
  ra <- NULL
  featureAnnotationColor <- .heatmapAnnColor(SCE, slot = "rowData",
                                             custom = featureAnnotationColor,
                                             palette = palette)
  if(ncol(rowData(SCE)) > 0)
    ra <- ComplexHeatmap::rowAnnotation(df = rowData(SCE),
                                        col = featureAnnotationColor)
  ### Set split variable
  cs <- NULL
  if (!is.null(colSplitBy)) cs <- colData(SCE)[colSplitBy]
  rs <- NULL
  if (!is.null(rowSplitBy)) rs <- rowData(SCE)[rowSplitBy]
  ###
  if (!is.null(colGap)) {
    if (!inherits(colGap, "unit"))
      stop("`colGap` has to be 'unit' object. Try `grid::unit(", colGap,
           ", 'mm')`.")
  }
  else colGap <- grid::unit(0, 'mm')
  if (!is.null(rowGap)) {
    if (!inherits(rowGap, "unit"))
      stop("`rowGap` has to be 'unit' object. Try `grid::unit(", rowGap,
           ", 'mm')`.")
  }
  else rowGap <- grid::unit(0, 'mm')
  
  if (!is.null(useAssay)) name <- useAssay
  else name <- useReducedDim
  hm <- ComplexHeatmap::Heatmap(mat, name = name, left_annotation = ra,
                                top_annotation = ca, col = colorScheme,
                                row_split = rs, column_split = cs,
                                row_title = rowTitle, column_title = colTitle,
                                show_row_names = rowLabel,
                                row_names_gp = grid::gpar(fontsize = rowLabelSize),
                                show_row_dend = rowDend,
                                show_column_dend = colDend,
                                row_dend_reorder = TRUE,
                                cluster_columns = cluster_columns,
                                cluster_rows = cluster_rows,
                                show_column_names = colLabel,
                                column_names_gp = grid::gpar(fontsize = colLabelSize),
                                row_gap = rowGap, column_gap = colGap,
                                border = border,
                                ...)
  return(hm)
}

.mergeAnnotationDF <- function(origin, external) {
  if (!is.null(external)) {
    external <- external[match(rownames(origin), rownames(external)), ,drop = FALSE]
    origin <- cbind(origin, external)
  }
  return(origin)
}

.heatmapAnnColor <- function(inSCE, slot = c("colData", "rowData"),
                             custom = NULL, palette = palette) {
  slot <- match.arg(slot)
  if (!is.null(custom) && !is.list(custom))
    stop("'cellAnnotationColor' or 'featureAnnotationColor' must be a list.")
  if (is.null(custom)) custom <- list()
  if (slot == "colData") data <- SummarizedExperiment::colData(inSCE)
  if (slot == "rowData") data <- SummarizedExperiment::rowData(inSCE)
  todoNames <- colnames(data)
  todoNames <- todoNames[!todoNames %in% names(custom)]
  newColor <- lapply(todoNames, function(n) {
    var <- data[[n]]
    if (is.factor(var)) categories <- levels(var)
    else categories <- unique(var)
    colors <- discreteColorPalette(length(categories), palette = palette)
    names(colors) <- categories
    return(colors)
  })
  names(newColor) <- todoNames
  custom <- c(custom, newColor)
  return(custom)
}
# Test
#logcounts(sceBatches) <- log1p(counts(sceBatches))
#plotSCEHeatmap2(sceBatches, "logcounts",
#                featureIndex = c("GCG1", "COX11", "INS1", "ND41"),
#                featureIndexBy = rowData(sceBatches)$feature_name,
#                cellIndex = c("reads.16087_", "Sample_1073_",
#                              "reads.29330_", "Sample_801_"),
#                cellIndexBy = paste0(colnames(sceBatches), "_"),
#                rowLabel = "feature_name", rowDend = FALSE,
#                cluster_rows = FALSE, colLabel = TRUE, cluster_columns = FALSE,
#                colDataName = c("batch", "cell_type"), aggregateCol = c("cell_type", "batch"))
#sce <-plotSCEHeatmap2(sceBatches, aggregateCol = "batch")
#plotSCEHeatmap2(sceBatches, aggregateCol = c("cell_type", "batch"))
#plotFindMarkerHeatmap(sce, log2fcThreshold = 0, minClustExprPerc = 0.4,
#                      maxCtrlExprPerc = 0.5)
#plotFindMarkerHeatmap(sce, log2fcThreshold = 0, minClustExprPerc = 0.4,
#                      maxCtrlExprPerc = 0.5,
#                      aggregateRow = "marker")
#plotSCEDimReduceColData(sce, "cluster", "UMAP")
CellVarColor <- function(inSCE, var,
                         palette = c("ggplot", "random", "celda"),
                         seed = 12345, ...) {
  var <- .manageCellVar(inSCE, var = var)
  palette <- match.arg(palette)
  if (is.factor(var)) uniqVar <- levels(var)
  else uniqVar <- unique(var)
  colors <- discreteColorPalette(length(uniqVar), palette = palette, seed = seed, ...)
  names(colors) <- uniqVar
  return(colors)
}

