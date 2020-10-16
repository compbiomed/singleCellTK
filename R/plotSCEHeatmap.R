#' Extract columns from row/colData and transfer to factors
#' @param inSCE \linkS4class{SingleCellExperiment} inherited object.
#' @param axis Choose from \code{"col"} or \code{"row"}.
#' @param columns character vector. The columns needed to be extracted. If
#' \code{NULL}, will return an empty \code{data.frame} with matched row
#' names. Default \code{NULL}.
#' @param index Valid index to subset the col/row.
#' @return A \code{data.frame} object.
.extractSCEAnnotation <- function(inSCE, axis = NULL, columns = NULL,
                                  index = NULL){
    if(is.null(axis) || !axis %in% c('col', 'row')){
        stop("axis should be 'col' or 'row'.")
    } else if(axis == 'col'){
        data <- SummarizedExperiment::colData(inSCE)
    } else if(axis == 'row'){
        data <- SummarizedExperiment::rowData(inSCE)
    }
    if(!is.null(index)){
        data <- data[index, , drop = FALSE]
    }
    if(is.null(columns)){
        return(data.frame(row.names = rownames(data)))
    } else {
        df <- data[, columns, drop = FALSE]
        for(i in colnames(df)){
            if(is.character(df[[i]]) || is.logical(df[[i]])){
                # Only converting character and logical columns, but not integer
                # cluster labels..
                df[[i]] <- as.factor(df[[i]])
            }
        }
        return(df)
    }
}

#' Generate distinct colors for all categorical col/rowData entries.
#' Character columns will be considered as well as all-integer columns. Any
#' column with all-distinct values will be excluded.
#' @param inSCE \linkS4class{SingleCellExperiment} inherited object.
#' @param axis Choose from \code{"col"} or \code{"row"}.
#' @param colorGen A function that generates color code vector by giving an
#' integer for the number of colors. Alternatively,
#' \code{\link[grDevices]{rainbow}}. Default \code{\link{distinctColors}}.
#' @return A \code{list} object containing distinct colors mapped to all
#' possible categorical entries in \code{rowData(inSCE)} or
#' \code{colData(inSCE)}.
#' @author Yichen Wang
dataAnnotationColor <- function(inSCE, axis = NULL,
                                colorGen = distinctColors){
    if(!is.null(axis) && axis == 'col'){
        data <- SummarizedExperiment::colData(inSCE)
    } else if(!is.null(axis) && axis == 'row'){
        data <- SummarizedExperiment::rowData(inSCE)
    } else {
        stop('please specify "col" or "row"')
    }
    nColor <- 0
    for(i in names(data)){
        if(length(grep('counts', i)) > 0){
            next
        }
        column <- stats::na.omit(data[[i]])
        if(is.numeric(column)){
            if(!all(as.integer(column) == column)){
                # Temporarily the way to tell whether numeric categorical
                next
            }
        }
        if(is.factor(column)){
            uniqLevel <- levels(column)
        } else {
            uniqLevel <- unique(column)
        }
        if(!length(uniqLevel) == nrow(data)){
            # Don't generate color for all-uniq annotation (such as IDs/symbols)
            nColor <- nColor + length(uniqLevel)
        }
    }
    allColors <- colorGen(nColor)
    nUsed <- 0
    allColorMap <- list()
    for(i in names(data)){
        if(length(grep('counts', i)) > 0){
            next
        }
        column <- stats::na.omit(data[[i]])
        if(is.numeric(column)){
            if(!all(as.integer(column) == column)){
                # Temporarily the way to tell whether numeric categorical
                next
            }
        }
        if(is.factor(column)){
            uniqLevel <- levels(column)
        } else {
            uniqLevel <- unique(column)
        }
        if(!length(uniqLevel) == nrow(data)){
            subColors <- allColors[(nUsed+1):(nUsed+length(uniqLevel))]
            names(subColors) <- uniqLevel
            allColorMap[[i]] <- subColors
            nUsed <- nUsed + length(uniqLevel)
        }
    }
    return(allColorMap)
}

#' Plot heatmap of using data stored in SingleCellExperiment Object
#' @param inSCE \linkS4class{SingleCellExperiment} inherited object.
#' @param useAssay character. A string indicating the assay name that
#' provides the expression level to plot.
#' @param featureIndex A vector that can subset the input SCE object by rows
#' (features). Alternatively, it can be a vector identifying features in
#' another feature list indicated by \code{featureIndexBy}. Default \code{NULL}.
#' @param cellIndex A vector that can subset the input SCE object by columns
#' (cells). Alternatively, it can be a vector identifying cells in another
#' cell list indicated by \code{featureIndexBy}. Default \code{NULL}.
#' @param featureIndexBy A single character specifying a column name of
#' \code{rowData(inSCE)}, or a vector of the same length as \code{nrow(inSCE)},
#' where we search for the non-rowname feature indices. Default
#' \code{"rownames"}.
#' @param cellIndexBy A single character specifying a column name of
#' \code{colData(inSCE)}, or a vector of the same length as \code{ncol(inSCE)},
#' where we search for the non-rowname cell indices. Default \code{"rownames"}.
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
#' @param rowDataName character. The column name(s) in \code{rowData} that need
#' to be added to the annotation. Default \code{NULL}.
#' @param colDataName character. The column name(s) in \code{colData} that need
#' to be added to the annotation. Default \code{NULL}.
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
#' @param rowDend Whether to display row dendrogram. Default \code{TRUE}.
#' @param colDend Whether to display column dendrogram. Default \code{TRUE}.
#' @param scale Whether to perform z-score scaling on each row. Default
#' \code{TRUE}.
#' @param trim A 2-element numeric vector. Values outside of this range will be
#' trimmed to their nearst bound. Default \code{c(-2, 2)}
#' @param title The main title of the whole plot. Default \code{"SCE Heatmap"}
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
#' @return A \code{\link[ComplexHeatmap]{Heatmap}} object
#' @export
#' @author Yichen Wang
plotSCEHeatmap <- function(inSCE, useAssay = 'logcounts', featureIndex = NULL,
    cellIndex = NULL, featureIndexBy = 'rownames', cellIndexBy = 'rownames',
    featureAnnotations = NULL, cellAnnotations = NULL,
    featureAnnotationColor = NULL, cellAnnotationColor = NULL,
    rowDataName = NULL, colDataName = NULL, rowSplitBy = NULL,
    colSplitBy = NULL, rowLabel = FALSE, colLabel = FALSE, rowDend = TRUE,
    colDend = TRUE, scale = TRUE, trim = c(-2, 2),
    title = 'SCE Heatmap', rowTitle = 'Genes', colTitle = 'Cells',
    rowGap = grid::unit(0, 'mm'), colGap = grid::unit(0, 'mm'), border = FALSE,
    colorScheme = NULL, ...){
    # Check input
    if(!inherits(inSCE, "SingleCellExperiment")){
        stop('Input object is not a valid SingleCellExperiment object.')
    }
    if(!useAssay %in% SummarizedExperiment::assayNames(inSCE)){
        stop('Specified assay does not exist in input SCE object')
    }
    if(!all(rowDataName %in% names(SummarizedExperiment::rowData(inSCE)))){
        notIn <- !rowDataName %in% names(SummarizedExperiment::rowData(inSCE))
        notIn <- rowDataName[notIn]
        stop('rowDataName - Specified columns: ', paste(notIn, collapse = ', '),
             ', not found. ')
    }
    if(!all(colDataName %in% names(SummarizedExperiment::colData(inSCE)))){
        notIn <- !colDataName %in% names(SummarizedExperiment::colData(inSCE))
        notIn <- colDataName[notIn]
        stop('colDataName - Specified columns: ', paste(notIn, collapse = ', '),
             ', not found. ')
    }
    if(!is.null(rowSplitBy) &&
       any(!rowSplitBy %in% c(rowDataName, names(featureAnnotations)))){
        notIn <- !rowSplitBy %in% c(names(SummarizedExperiment::rowData(inSCE)),
                                    featureAnnotations)
        notIn <- rowSplitBy[notIn]
        stop('rowSplitBy - Specified columns: ', paste(notIn, collapse = ', '),
             ', not found. ')
    }
    if(!is.null(colSplitBy) &&
       any(!colSplitBy %in% c(colDataName, names(cellAnnotations)))){
        notIn <- !colSplitBy %in% c(names(SummarizedExperiment::colData(inSCE)),
                                    cellAnnotations)
        notIn <- colSplitBy[notIn]
        stop('colSplitBy - Specified columns: ', paste(notIn, collapse = ', '),
             ', not found. ')
    }
    if(is.null(featureIndex)){
        featureIndex <- seq_len(nrow(inSCE))
    } else {
        if(is.character(featureIndexBy) && length(featureIndexBy) == 1){
            if(!featureIndexBy == 'rownames'){
                # Search by a column in rowData
                featureIndex <- celda::retrieveFeatureIndex(featureIndex,
                                                            inSCE,
                                                            featureIndexBy)
            }
        } else if(length(featureIndexBy) == nrow(inSCE)){
            # featureIndexBy is vector or single-col/row matrix
            featureIndex <- celda::retrieveFeatureIndex(featureIndex,
                                                        featureIndexBy,
                                                        '')
        } else {
            stop('Given "featureIndexBy" not valid. Please give a single ',
                 'character to specify a column in rowData(inSCE) or a vector ',
                 'as long as nrow(inSCE) where you search for "featureIndex".')
        }
    }
    ### Force index as numeric
    if(is.character(featureIndex)){
        featureIndex <- which(rownames(inSCE) %in% featureIndex)
    } else if(is.logical(featureIndex)){
        featureIndex <- which(featureIndex)
    }
    if(is.null(cellIndex)){
        cellIndex <- seq_len(ncol(inSCE))
    } else {
        if(is.character(cellIndexBy) && length(cellIndexBy) == 1){
            if(!cellIndexBy == 'rownames'){
                # Search by a column in colData
                if(!cellIndexBy %in%
                   names(SummarizedExperiment::colData(inSCE))){
                    stop('"cellIndexBy": ', cellIndexBy, ' is not a column of ',
                         'colData(inSCE)')
                }
                searchIn <- SummarizedExperiment::colData(inSCE)[[cellIndexBy]]
                cellIndex <- celda::retrieveFeatureIndex(cellIndex,
                                                            searchIn,
                                                            '')
            }
        } else if(length(cellIndexBy) == ncol(inSCE)){
            # featureIndexBy is vector or single-col/row matrix
            cellIndex <- celda::retrieveFeatureIndex(cellIndex,
                                                        cellIndexBy,
                                                        '')
        } else {
            stop('Given "cellIndexBy" not valid. Please give a single ',
                 'character to specify a column in colData(inSCE) or a vector ',
                 'as long as ncol(inSCE) where you search for "cellIndex".')
        }
    }
    ### Force index as numeric
    if(is.character(cellIndex)){
        cellIndex <- which(colnames(inSCE) %in% cellIndex)
    } else if (is.logical(cellIndex)){
        cellIndex <- which(cellIndex)
    }
    ## Customized row text labeling
    rowLabelText <- rownames(inSCE)[featureIndex]
    if(!is.logical(rowLabel)){
        if(is.character(rowLabel) && length(rowLabel) == 1){
            if(!rowLabel %in% names(SummarizedExperiment::rowData(inSCE))){
                stop('"rowLabel": ', rowLabel, ' is not a column of ',
                     'rowData(inSCE).')
            }
            rowLabelText <- SummarizedExperiment::rowData(inSCE)[featureIndex,
                                                                 rowLabel]
            rowLabel <- TRUE
        } else if(length(rowLabel) == nrow(inSCE)){
            rowLabelText <- rowLabel[featureIndex]
            rowLabel <- TRUE
        } else if(length(rowLabel) == length(featureIndex)){
            rowLabelText <- rowLabel
            rowLabel <- TRUE
        } else {
            stop('Invalid "rowLabel". Use TRUE/FALSE, a column name of ',
                 'rowData(inSCE), or a vector as the same length of ',
                 'nrow(inSCE) or the subsetted number of features.')
        }
    }
    ## Customized col text labeling
    colLabelText <- colnames(inSCE)[cellIndex]
    if(!is.logical(colLabel)){
        if(is.character(colLabel) && length(colLabel) == 1){
            if(!colLabel %in% names(SummarizedExperiment::colData(inSCE))){
                stop('"colLabel": ', colLabel, ' is not a column of ',
                     'colData(inSCE).')
            }
            colLabelText <- SummarizedExperiment::colData(inSCE)[cellIndex,
                                                                 colLabel]
            colLabel <- TRUE
        } else if(length(colLabel) == ncol(inSCE)){
            colLabelText <- colLabel[cellIndex]
            colLabel <- TRUE
        } else if(length(colLabel) == length(cellIndex)){
            colLabelText <- colLabel
            colLabel <- TRUE
        } else {
            stop('Invalid "colLabel". Use TRUE/FALSE, a column name of ',
                 'colData(inSCE), or a vector as the same length of ',
                 'ncol(inSCE) or the subsetted number of cells.')
        }
    }

    inSCE <- inSCE[featureIndex, cellIndex]

    if (!is.null(trim) && length(trim) != 2) {
        stop("'trim' should be a 2 element vector specifying the lower",
             " and upper boundaries")
    }
    if(0 %in% dim(inSCE)){
        stop('Given indices specified 0-dim')
    }
    if(!is.null(featureAnnotations)){
        if(!all(rownames(inSCE) %in% rownames(featureAnnotations))){
            stop('Incomplete feature names in `featureAnnotations')
        } else {
            featureAnnotations <-
                featureAnnotations[rownames(inSCE), , drop = FALSE]
        }
    }
    if(!is.null(cellAnnotations)){
        if(!all(colnames(inSCE) %in% rownames(cellAnnotations))){
            stop('Incomplete cell names in cellAnnotations')
        } else {
            cellAnnotations <- cellAnnotations[colnames(inSCE), , drop = FALSE]
        }
    }

    # Extract
    mat <- as.matrix(SummarizedExperiment::assay(inSCE, useAssay))
    ## rowData info
    rowDataExtract <- .extractSCEAnnotation(inSCE, 'row', rowDataName)
    rowDataColor <- dataAnnotationColor(inSCE, 'row')
    if(is.null(rowDataName)){
        rowDataColor <- NULL
    } else {
        # Have to do an extraction because continuous values won't be in
        # rowDataColor
        rowDataColor <- rowDataColor[rowDataName[rowDataName %in%
                                                     names(rowDataColor)]]
    }
    if(!is.null(featureAnnotationColor)){
        add <- setdiff(names(rowDataColor), names(featureAnnotationColor))
        featureAnnotationColor <- c(rowDataColor[add], featureAnnotationColor)
    } else {
        featureAnnotationColor <- rowDataColor
    }
    ## colData info
    colDataExtract <- .extractSCEAnnotation(inSCE, 'col', colDataName)
    colDataColor <- dataAnnotationColor(inSCE, 'col')
    if(is.null(colDataName)){
        colDataColor <- NULL
    } else {
        colDataColor <- colDataColor[colDataName[colDataName %in%
                                                     names(colDataColor)]]
    }
    if(!is.null(cellAnnotationColor)){
        add <- setdiff(names(colDataColor), names(cellAnnotationColor))
        cellAnnotationColor <- c(colDataColor[add], cellAnnotationColor)
    } else {
        cellAnnotationColor <- colDataColor
    }
    ## Merge with extra annotations
    if(is.null(featureAnnotations)){
        featureAnnotations <- rowDataExtract
    } else {
        featureAnnotations <- data.frame(rowDataExtract, featureAnnotations)
    }
    if(is.null(cellAnnotations)){
        cellAnnotations <- colDataExtract
    } else {
        cellAnnotations <- data.frame(colDataExtract, cellAnnotations)
    }
    # Data process
    if(isTRUE(scale)){
        mat <- as.matrix(computeZScore(mat))
    }
    if (!is.null(trim)) {
        mat <- trimCounts(mat, trim)
    }
    # Plot
    if(is.null(colorScheme)){
        if(!is.null(trim)){
            colorScheme <- circlize::colorRamp2(c(trim[1], 0, trim[2]),
                                                c('blue', 'white', 'red'))
        } else {
            colorScheme <- circlize::colorRamp2(c(min(mat),
                                                  (max(mat) + min(mat))/2,
                                                  max(mat)),
                                                c('blue', 'white', 'red'))
        }

    } else {
        if(!is.function(colorScheme)){
            stop('`colorScheme` must be a function generated by ',
                 'circlize::colorRamp2')
        }
        breaks <- attr(colorScheme, 'breaks')
        if(breaks[1] != min(trim) || breaks[length(breaks)] != max(trim)){
            stop('Breaks of `colorScheme` do not match with `trim`.')
        }
    }
    if(dim(featureAnnotations)[2] > 0){
        ra <- ComplexHeatmap::rowAnnotation(df = featureAnnotations,
                                            col = featureAnnotationColor)
    } else {
        ra <- NULL
    }
    if(dim(cellAnnotations)[2] > 0){
        ca <- ComplexHeatmap::HeatmapAnnotation(df = cellAnnotations,
                                                col = cellAnnotationColor)
    } else {
        ca <- NULL
    }
    if(!is.null(rowSplitBy)){
        rs <- featureAnnotations[rowSplitBy]
    } else {
        rs <- NULL
    }
    if(!is.null(colSplitBy)){
        cs <- cellAnnotations[colSplitBy]
    } else {
        cs <- NULL
    }
    if(!is.null(rowGap)) {
      if(inherits(rowGap, "unit")){
        rowGap <- rowGap
      } else if (is.numeric(rowGap)) {
        warning("rowGap is given a numeric value. Using 'mm' as the unit")
        rowGap <- grid::unit(rowGap, 'mm')
      } else {
        stop("Given value for 'rowGap' not understandable.")
      }
    } else {
      rowGap <- grid::unit(0, 'mm')
    }
    if(!is.null(colGap)) {
      if (inherits(colGap, "unit")) {
        colGap <- colGap
      } else if(is.numeric(colGap)){
        warning("colGap is given a numeric value. Using 'mm' as the unit")
        colGap <- grid::unit(colGap, 'mm')
      } else {
        stop("Given value for 'colGap' not understandable.")
      }
    } else {
      colGap <- grid::unit(0, 'mm')
    }
    rownames(mat) <- rowLabelText
    colnames(mat) <- colLabelText
    hm <- ComplexHeatmap::Heatmap(mat, name = useAssay, left_annotation = ra,
                                  top_annotation = ca, col = colorScheme,
                                  row_split = rs, column_split = cs,
                                  row_title = rowTitle, column_title = colTitle,
                                  show_row_names = rowLabel,
                                  show_row_dend = rowDend,
                                  show_column_names = colLabel,
                                  show_column_dend = colDend,
                                  row_gap = rowGap, column_gap = colGap,
                                  border = border,
                                  ...)
    #HM <- ComplexHeatmap::draw(hm, column_title = title)
    return(hm)
}
