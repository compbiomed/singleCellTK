#' Summarize an assay in a \linkS4class{SingleCellExperiment}
#'
#' Creates a table of summary metrics from an input
#' \linkS4class{SingleCellExperiment}
#'
#' @param inSCE Input SingleCellExperiment object.
#' @param useAssay Indicate which assay to summarize. If \code{NULL}, then the
#' first assay in \code{inSCE} will be used. Default \code{NULL}.
#' @param sampleVariableName Variable name in \code{colData} denoting which
#' sample each cell belongs to. If \code{NULL}, all cells will be assumed
#' to come from the same sample. Default \code{"sample"}.
#'
#' @return A data.frame object of summary metrics.
#' @export
#' @examples
#' data("mouseBrainSubsetSCE")
#' summarizeSCE(mouseBrainSubsetSCE, sample = NULL)
#' @importFrom SummarizedExperiment assays colData
summarizeSCE <- function(inSCE, useAssay = NULL, sampleVariableName = NULL){

  if(is.null(useAssay)) {
    useAssay <- names(assays(inSCE))[1]
  }

  if(is.null(sampleVariableName)) {
    sampleVariable <- rep("Sample", ncol(inSCE))
  } else {
    if(!(sampleVariableName %in% colnames(colData(inSCE)))) {
      stop("'", sampleVariableName, "' was not found in the 'colData' of ",
           "'inSCE'.")
    }
    sampleVariable <- colData(inSCE)[,sampleVariableName]
  }

  numCells <- table(sampleVariable)
  var <- colSums(SummarizedExperiment::assay(inSCE, useAssay))
  meanCounts <- stats::aggregate(var, by = list(sampleVariable), FUN = mean)
  medianCounts <- stats::aggregate(var, by = list(sampleVariable),
                                   FUN = stats::median)
  var2 <- colSums(SummarizedExperiment::assay(inSCE, useAssay) > 0)
  meanDetected <- stats::aggregate(var2, by = list(sampleVariable), FUN = mean)
  medianDetected <- stats::aggregate(var2, by = list(sampleVariable),
                                     FUN = stats::median)

  df <- data.frame("Sample" = names(numCells),
                   "Number of Cells" = as.integer(round(as.numeric(numCells))),
                   "Mean counts per cell" = as.integer(round(meanCounts[,2])),
                   "Median counts per cell" = as.integer(round(medianCounts[,2])),
                   "Mean features detected per cell" = as.integer(round(meanDetected[,2])),
                   "Median features detected per cell" = as.integer(round(medianDetected[,2])),
                   stringsAsFactors = FALSE, check.names = FALSE)
  return(df)
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

#' Generate given number of color codes
#'
#' @description Three different generation methods are wrapped, including
#' \code{\link[celda]{distinctColors}},
#' [randomcoloR](SCTK_PerformingQC_Cell_V3.Rmd) and the \code{ggplot}
#' default color generation.
#' @param n An integer, the number of color codes to generate.
#' @param palette A single character string. Select the method, available
#' options are \code{"ggplot"}, \code{"celda"} and \code{"random"}. Default
#' \code{"random"}.
#' @param seed An integer. Set the seed for random process that happens only in
#' "random" generation. Default \code{12345}.
#' @param ... Other arguments that are passed to the internal function,
#' according to the method selected.
#' @return A character vector of \code{n} hex color codes.
#' @examples
#' discreteColorPalette(n = 3)
#' @export
discreteColorPalette <- function(n, palette = c("random", "ggplot", "celda"),
                                 seed = 12345, ...) {
  palette <- match.arg(palette)
  # Generate the colors
  if (palette == "random") {
    withr::with_seed(seed, {
      # Using randomcoloR original codes here
      km <- stats::kmeans(colorSpace, n, iter.max = 20)
      colors <- unname(colorspace::hex(colorspace::LAB(km$centers)))
      colors <- colors[order(colors)]
    })
  } else if (palette == "ggplot") {
    hues <- seq(15, 375, length = n + 1)
    colors <- grDevices::hcl(h = hues, l = 65, c = 100)[seq_len(n)]
  } else if (palette == "celda") {
    colors <- distinctColors(n, ...)
  } else {
    stop("Given palette '", palette, "' is not supported. Please choose from ",
         "'ggplot', 'celda' or 'random'.")
  }
  return(colors)
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

## Convert a matrix to a sparse matrix and preserve column/row names
.convertToMatrix <- function(x) {
  cn <- colnames(x)
  rn <- rownames(x)
  limit <- (2^32/2-1)
  dimN <- dim(x)
  chuS <- floor(floor(limit/dimN[1])) # size of chunk
  chuN <- ceiling(dimN[2]/chuS) # number of chunks
  Mat <- list()

  for (i in seq_len(chuN)) {
    start <- (i-1)*chuS + 1
    end <- min(i*chuS, dimN[2])
    if (methods::is(x, 'DelayedMatrix')) {
      # Efficient way to convert DelayedArray to dgCMatrix
      Mat[[i]] <- methods::as(x[, start:end], "Matrix")
    } else {
      # Convert dgTMatrix to dgCMatrix
      Mat[[i]] <- methods::as(x[, start:end], "dgCMatrix")
    }
  }
  x <- do.call(base::cbind, Mat)
  colnames(x) <- cn
  rownames(x) <- rn

  return(x)
}

#' Resolve duplicated feature names in a matrix
#'
#' Adds '-1', '-2', ... '-i' to multiple duplicated feature names
#' @param countmat matrix, with row names as feature names that need to be
#' resolved.
#' @return The same matrix as input with rowname duplication resolved.
featureNameDedup <- function(countmat){
    if(!inherits(rownames(countmat), "character")){
        stop("No character feature name found.")
    }
    gene.table <- table(rownames(countmat))
    gene.duplicates <- gene.table[gene.table > 1]
    gene.duplicates.names <- names(gene.duplicates)
    empty <- character(0)
    if (!identical(empty, gene.duplicates.names)){
        for (genename in gene.duplicates.names){
            genename <- gsub(" (1 of many)", "", genename, fixed=TRUE)
            indices <- which(grepl(genename, rownames(countmat)))
            for (i in seq_along(indices)){
                rownames(countmat)[indices[i]] <- paste0(genename, "-", i)
            }
        }
    }
    return(countmat)
}


#' Retrieve cell/feature index by giving identifiers saved in col/rowData
#'
#' @description Originally written in \code{\link[celda]{retrieveFeatureIndex}}.
#' Modified for also retrieving cell indices and only working for
#' \linkS4class{SingleCellExperiment} object. This will return indices of
#' features among the \code{rowData}/\code{colData}. Partial matching (i.e.
#' grepping) can be used.
#' @param inSCE Input \linkS4class{SingleCellExperiment} object. Required
#' @param IDs Character vector of identifiers for features or cells to find in
#' \code{rowData} or \code{colData} of \code{inSCE}
#' @param axis A character scalar to specify whether to search for features or
#' cells. Use \code{"row"}, \code{"feature"} or \code{"gene"} for features;
#' \code{"col"} or \code{"cell"} for cells.
#' @param by Character. In which column to search for features/cells in
#' \code{rowData}/\code{colData}. Default \code{NULL} for search the
#' \code{rownames}/\code{colnames}
#' @param exactMatch A logical scalar. Whether to only identify exact matches
#' or to identify partial matches using \code{\link{grep}}. Default \code{TRUE}
#' @param firstMatch A logical scalar. Whether to only identify the first
#' matches or to return all plausible matches. Default \code{TRUE}
#' @examples
#' data(scExample, package = "singleCellTK")
#' retrieveSCEIndex(inSCE = sce, IDs = "ENSG00000205542",
#'  axis = "row")
#' @return A unique, non-NA numeric vector of indices for the matching
#' features/cells in \code{inSCE}.
#' @author Yusuke Koga, Joshua Campbell
#' @export
retrieveSCEIndex <- function(inSCE, IDs, axis, by = NULL,
                             exactMatch = TRUE, firstMatch = TRUE){
  if(!inherits(inSCE, "SingleCellExperiment")){
    stop("`inSCE` should inherits from a SingleCellExperiment object.")
  }
  if(!axis %in% c('row', 'col', 'cell', 'feature', 'gene')){
    stop("Invalid axis specification")
  }
  if(axis %in% c('row', 'feature', 'gene')){
    data <- SummarizedExperiment::rowData(inSCE)
  } else {
    data <- SummarizedExperiment::colData(inSCE)
  }

  if(!is.null(by)){
    if(!by %in% colnames(data)){
      stop('"', by, '" annotation not found for "', axis, '".')
    }
    search <- data[[by]]
  } else {
    if (is.null(rownames(data))) {
      stop("'rownames' of 'inSCE' are 'NULL'. Please set 'rownames' or change",
           " 'by' to search a different annotation of 'inSCE'.")
    }
    search <- rownames(data)
  }
  Indices <- numeric()
  notFound <- numeric()
  if (!isTRUE(exactMatch)) {
    for(i in seq_along(IDs)){
      g <- grep(IDs[i], search)
      if(length(g) == 0){
        notFound <- c(notFound, i)
      } else if (length(g) == 1){
        Indices <- c(Indices, g)
      } else if(length(g) > 1){
        if(isTRUE(firstMatch)){
          Indices <- c(Indices, g[1])
        } else {
          Indices <- c(Indices, g)
        }
      }
    }
    dupMatched <- search[unique(Indices[duplicated(Indices)])]
  } else {
    if(isTRUE(firstMatch)){
      Indices <- match(IDs, search)
      notFound <- which(is.na(Indices))
      Indices <- Indices[!is.na(Indices)]
      dupMatched <- search[unique(Indices[duplicated(Indices)])]
    } else {
      Indices <- which(search %in% IDs)
      notFound <- which(!IDs %in% search)
      dupMatched <- unique(IDs[duplicated(IDs[IDs %in% search])])
    }
  }
  if(length(notFound) == length(IDs)){
    if (isTRUE(exactMatch)) {
      warning("None of the provided features had matching items in '", by,
              "'. Check the spelling or try setting `exactMatch` to FALSE.")
    } else {
      warning("None of the provided features had matching items in '", by,
              "'. Check the spelling and make sure `by` is set to the ",
              "appropriate annotation.")
    }
  } else if(length(notFound) > 0){
    warning("The following IDs were not present in specified annotation: \n'",
            paste(IDs[notFound], collapse = "', '"), "'")
  }
  if(length(dupMatched) > 0){
    warning("Each of the following entries from '", by, "' was matched by ",
            "multiple queries in 'IDs': \n'",
            paste(dupMatched, collapse = "', '"), "'")
  }
  Indices <- unique(Indices)
  return(Indices)
}


# backup or restore 'factor' columns in a dataframe
# (for use in col/row annotation editor)
.manageFactor <- function(df, operation = "backup"){
  if(operation == "backup"){
    data <- list()
    data$data_type <- list()
    data$df <- df
    for (i in seq_along(colnames(data$df))) {
      data$data_type[[colnames(data$df)[i]]] <- c(typeof(data$df[,i]),
                                                  is.factor(data$df[,i]))
    }
    data$df <- .convertFactorToCharacter(df)
  }
  else if(operation == "restore"){
    data <- df
    for (i in seq(length(colnames(data$df)))) {
      if(data$data_type[[i]][2] == TRUE){
        data$df[,i] <- as.factor(data$df[,i])
      }
    }
  }
  return(data)
}

#converts the columns of a dataframe from factor to character
#for use with .manageFactor
.convertFactorToCharacter <- function(df){
  for (i in seq(length(colnames(df)))) {
    if (is.factor(df[, i])) {
      df[,i] = as.character(df[,i])
    }
  }
  return(df)
}
