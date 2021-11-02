#' computeHeatmap
#' The computeHeatmap method computes the heatmap visualization for a set
#'  of features against a set of dimensionality reduction components. This
#'  method uses the heatmap computation algorithm code from \code{Seurat} but
#'  plots the heatmap using \code{ComplexHeatmap} and \code{cowplot} libraries.
#' @param inSCE Input \code{SingleCellExperiment} object.
#' @param useAssay The assay to use for heatmap computation.
#' @param dims Specify the number of dimensions to use for heatmap. Default
#' \code{10}.
#' @param nfeatures Specify the number of features to use for heatmap. Default
#' is \code{30}.
#' @param cells Specify the samples/cells to use for heatmap computation.
#' Default is \code{NULL} which will utilize all samples in the assay.
#' @param reduction Specify the reduction slot in the input object. Default
#' is \code{"pca"}.
#' @param disp.min Specify the minimum dispersion value to use for floor
#' clipping of assay values. Default is \code{-2.5}.
#' @param disp.max Specify the maximum dispersion value to use for ceiling
#' clipping of assay values. Default is \code{2.5}.
#' @param balanced Specify if the number of of up-regulated and down-regulated
#' features should be balanced. Default is \code{TRUE}.
#' @param nCol Specify the number of columns in the output plot. Default
#' is \code{NULL} which will auto-compute the number of columns.
#' @param externalReduction Specify an external reduction if not present in
#' the input object. This external reduction should be created
#' using \code{CreateDimReducObject} function.
#'
#' @return Heatmap plot object.
#' @export
computeHeatmap <- function(inSCE,
                        useAssay,
                        dims = 10,
                        nfeatures = 30,
                        cells = NULL,
                        reduction = 'pca',
                        disp.min = -2.5,
                        disp.max = 2.5,
                        balanced = TRUE,
                        nCol = NULL,
                        externalReduction = NULL){
  object <- convertSCEToSeurat(inSCE, scaledAssay = useAssay)
  slot <- "scale.data"
  assays <- NULL
  ncol <- NULL
  projected <- FALSE
  dims <- seq(dims)

  if(!is.null(externalReduction)){
    if(reduction == "pca"){
      object@reductions <- list(pca = externalReduction)
    }
    else{
      object@reductions <- list(ica = externalReduction)
    }
  }

  if(length(x = dims) > 2){
    ncol <- 3
  }
  else{
    ncol <- length(x = dims)
  }
  #ncol <- ncol %||% ifelse(test = length(x = dims) > 2, yes = 3,
  #no = length(x = dims))

  #empty list = number of dims
  plots <- vector(mode = 'list', length = length(x = dims))

  #set rna
  assays <- Seurat::DefaultAssay(object = object)
  # assays <- assays %||% Seurat::DefaultAssay(object = object)

  #set disp.max
  # disp.max <- disp.max %||% ifelse(
  #   test = slot == 'scale.data',
  #   yes = 2.5,
  #   no = 6
  # )

  #no. of cells
  cells <- ncol(x = object)
  # cells <- cells %||% ncol(x = object)

  #for each dim get top cells
  cells <- lapply(
    X = dims,
    FUN = function(x) {
      cells <- Seurat::TopCells(
        object = object[[reduction]],
        dim = x,
        ncells = cells,
        balanced = balanced
      )
      if (balanced) {
        cells$negative <- rev(x = cells$negative)
      }
      cells <- unlist(x = unname(obj = cells))
      return(cells)
    }
  )

  #get top features against each dim
  features <- lapply(
    X = dims,
    FUN = Seurat::TopFeatures,
    object = object[[reduction]],
    nfeatures = nfeatures,
    balanced = balanced,
    projected = projected
  )

  #unlist all features in one vector
  features.all <- unique(x = unlist(x = features))


  #assays == 1
  features.keyed <- features.all

  #set default assay
  Seurat::DefaultAssay(object = object) <- assays

  #convert (_) to (-) as required by FetchData function below
  cells <- .convertToHyphen(cells)

  features.keyed <- .convertToHyphen(features.keyed)

  for (i in seq_len(length(dims))){
    features[[i]] <- .convertToHyphen(features[[i]])
  }

  # get assay data with only selected features (all dims) and
  # selected cells (all)
  data.all <- Seurat::FetchData(
    object = object,
    vars = unique(x = unlist(x = features.keyed)),
    cells = unique(x = unlist(x = cells)),
    slot = slot
  )

  #clip off values for heatmap
  data.all <- Seurat::MinMax(data = data.all, min = disp.min, max = disp.max)
  data.limits <- c(min(data.all), max(data.all))

  #draw heatmap for each dim
  for (i in seq_along(dims)) {
    dim.features <- c(features[[i]][[2]], rev(x = features[[i]][[1]]))
    dim.features <- rev(x = unlist(x = lapply(
      X = dim.features,
      FUN = function(feat) {
        return(grep(pattern = paste0(feat, '$'),
                    x = features.keyed, value = TRUE))
      }
    )))
    dim.cells <- cells[[i]]
    data.plot <- data.all[dim.cells, dim.features]
    hm <- suppressMessages(ComplexHeatmap::Heatmap(t(data.plot),
            show_row_dend = FALSE,
            show_column_dend = FALSE,
            cluster_rows = FALSE,
            cluster_columns = FALSE,
            show_column_names = FALSE,
            row_names_side = "left",
            show_heatmap_legend = FALSE))

      plots[[i]] <- grid::grid.grabExpr(suppressMessages(ComplexHeatmap::draw(hm)))
  }

  return(plots)
}
