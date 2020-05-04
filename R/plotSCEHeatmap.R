#' Extract columns from row/colData and transfer to factors
#' @param inSCE SingleCellExperiment object.
#' @param axis Choose from "col" or "row"
#' @param columns character vector, default NULL (return NULL). The columns
#' needed to be extracted.
#' @param index Valid index to subset the col/row.
#' @return A data.frame object
extractData <- function(inSCE, axis = NULL, columns = NULL, index = NULL){
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
#' @param inSCE SingleCellExperiment object.
#' @param axis Choose from "col" or "row"
#' @param colorGen function, default `rainbow`. A function that generates color
#' code vector by giving an integer for the number of colors. Alternatively,
#' `celda::distinctColors`.
#' @return An list object containing distinct colors mapped to all possible
#' categorical entries in row/colData.
dataAnnotationColor <- function(inSCE, axis = NULL, 
                                colorGen = grDevices::rainbow){
    if(!is.null(axis) && axis == 'col'){
        data <- SummarizedExperiment::colData(inSCE)
    } else if(!is.null(axis) && axis == 'row'){
        data <- SummarizedExperiment::rowData(inSCE)
    } else {
        stop('please specify "col" or "row"')
    }
    nColor <- 0
    for(i in names(data)){
        column <- data[[i]] # Double bracket to extract the vector
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
        nColor <- nColor + length(uniqLevel)
    }
    allColors <- colorGen(nColor)
    nUsed <- 0
    allColorMap <- list()
    for(i in names(data)){
        column <- data[[i]] # Double bracket to extract the vector
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
        subColors <- allColors[(nUsed+1):(nUsed+length(uniqLevel))]
        names(subColors) <- uniqLevel
        allColorMap[[i]] <- subColors
        nUsed <- nUsed + length(uniqLevel)
    }
    return(allColorMap)
}

#' Plot heatmap of using data stored in SingleCellExperiment Object
#' @param inSCE SingleCellExperiment object. 
#' @param useAssay character. A string indicating the assay name that 
#' provides the expression level to plot.
#' @param featureIndex vector, default `NULL`, that can subset the input SCE
#' object by rows (features)
#' @param cellIndex vector, default `NULL`, that can subset the input SCE
#' object by columns (cells)
#' @param featureAnnotations data.frame, with `rownames` containing all the
#' features going to be plotted, default NULL. Character columns should be
#' factors.
#' @param cellAnnotations data.frame, with `rownames` containing all the
#' cells going to be plotted, default NULL. Character columns should be
#' factors.
#' @param featureAnnotationColor A named list, default `NULL`. Customized color
#' setting for features. Should match the entries in the `featureAnnotations`
#' or `rowDataName`. For each entry, there should be a list/vector of colors
#' named with classes.
#' @param cellAnnotationColor A named list, default `NULL`. Customized color
#' setting for cells. Should match the entries in the `cellAnnotations` or
#' `colDataName`. For each entry, there should be a list/vector of colors named
#' with classes.
#' @param rowDataName character, default `NULL`. The column name(s) in
#' `rowData` that need to be added to the annotation.
#' @param colDataName character, default `NULL`. The column name(s) in
#' `colData` that need to be added to the annotation
#' @param rowSplitBy A single character, default `NULL`. Do semi-heatmap based
#' on the grouping of this annotation. Should exist in either `rowDataName` or
#' `names(featureAnnotations)`.
#' @param colSplitBy A single character, default `NULL`. Do semi-heatmap based
#' on the grouping of this annotation. Should exist in either `colDataName` or
#' `names(cellAnnotations)`.
#' @param rowLabel A single logical, default `FALSE`. Whether to display all
#' the feature names.
#' @param colLabel A single logical, default `FALSE`. Whether to display all
#' the cell names.
#' @param rowDend A single logical, default `True`. Whether to display row
#' dendrogram.
#' @param colDend A single logical, default `True`. Whether to display column
#' dendrogram.
#' @param scaleRow function, default `scale`. A function to be applied to each
#' row to scale the data.
#' @param trim A 2-element numeric vector, default `c(-2, 2)`. Values outside
#' of this range will be their nearst bound.
#' @param title character, default `"SCE Heatmap"`. The main title of the whole 
#' plot.
#' @param colorScheme function, default NULL. A function that generates color 
#' code by giving a value. Can be generated by `colorRamp2()`.
#' @param ... Arguments passed to ComplexHeatmap::Heatmap()
#' @return An ComplexHeatmap::Heatmap object
#' @export
plotSCEHeatmap <- function(inSCE, useAssay = 'logcounts', featureIndex = NULL,
    cellIndex = NULL, featureAnnotations = NULL, cellAnnotations = NULL,
    featureAnnotationColor = NULL, cellAnnotationColor = NULL,
    rowDataName = NULL, colDataName = NULL, rowSplitBy = NULL,
    colSplitBy = NULL, rowLabel = FALSE, colLabel = FALSE, rowDend = TRUE,
    colDend = TRUE, scaleRow = scale, trim = c(-2, 2), title = 'SCE Heatmap',
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
        featureIndex <- 1:nrow(inSCE)
    }
    if(is.null(cellIndex)){
        cellIndex <- 1:nrow(inSCE)
    }
    if (!is.null(scaleRow) && !is.function(scaleRow)) {
        stop("'scaleRow' needs to be of class 'function'")
    }
    if (!is.null(trim) && length(trim) != 2) {
        stop("'trim' should be a 2 element vector specifying the lower",
             " and upper boundaries")
    }
    inSCE <- inSCE[featureIndex, cellIndex]
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
    if(is.null(colorScheme)){
        colorScheme <- circlize::colorRamp2(c(trim[1], 0, trim[2]),
                                       c('blue', 'white', 'red'))
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
    # Extract
    mat <- SummarizedExperiment::assay(inSCE, useAssay)
    ## rowData info
    rowDataExtract <- extractData(inSCE, 'row', rowDataName)
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
    colDataExtract <- extractData(inSCE, 'col', colDataName)
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
    if (!is.null(scaleRow)) {
        cn <- colnames(mat)
        mat <- t(apply(mat, 1, scaleRow))
        colnames(mat) <- cn
    }
    if (!is.null(trim)) {
        trim <- sort(trim)
        mat[mat < trim[1]] <- trim[1]
        mat[mat > trim[2]] <- trim[2]
    }
    # Plot
    if(!is.null(featureAnnotations)){
        ra <- ComplexHeatmap::rowAnnotation(df = featureAnnotations,
                                            col = featureAnnotationColor)
    } else {
        ra <- NULL
    }
    if(!is.null(cellAnnotations)){
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
    hm <- ComplexHeatmap::Heatmap(mat, name = useAssay, left_annotation = ra,
                                  top_annotation = ca, col = colorScheme,
                                  row_split = rs, column_split = cs,
                                  row_gap = grid::unit(1, 'mm'),
                                  column_gap = grid::unit(1, 'mm'),
                                  row_title = 'Genes', column_title = 'Cells',
                                  show_row_names = rowLabel,
                                  show_row_dend = rowDend,
                                  show_column_names = colLabel,
                                  show_column_dend = colDend, ...)
    HM <- ComplexHeatmap::draw(hm, column_title = title)
    return(HM)
}
