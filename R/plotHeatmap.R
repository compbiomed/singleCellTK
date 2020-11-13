#' plotHeatmap
#'
#' @param inSCE 
#' @param useAssay 
#' @param dims 
#' @param nfeatures 
#' @param cells 
#' @param reduction 
#' @param disp.min 
#' @param disp.max 
#' @param balanced 
#' @param projected 
#' @param ncol 
#' @param fast 
#' @param raster 
#' @param slot 
#' @param assays 
#' @param combine 
#'
#' @return
#' @export
plotHeatmap <- function(inSCE,
                        useAssay,
                        dims = 1:2,
                        nfeatures = 30,
                        cells = NULL,
                        reduction = 'pca',
                        disp.min = -2.5,
                        disp.max = NULL,
                        balanced = TRUE,
                        projected = FALSE,
                        ncol = NULL,
                        fast = FALSE,
                        raster = TRUE,
                        slot = 'scale.data',
                        assays = NULL,
                        combine = TRUE){
  object <- convertSCEToSeurat(inSCE, scaledAssay = useAssay)
  
  #ncol is 3 now
  ncol <- ncol %||% ifelse(test = length(x = dims) > 2, yes = 3, no = length(x = dims))
  
  #empty 4 list = number of dims
  plots <- vector(mode = 'list', length = length(x = dims))
  
  #set rna
  assays <- assays %||% DefaultAssay(object = object)
  
  #set disp.max
  disp.max <- disp.max %||% ifelse(
    test = slot == 'scale.data',
    yes = 2.5,
    no = 6
  )
  
  #no. of cells = 2700
  cells <- cells %||% ncol(x = object)
  
  #for each dim get top cells
  cells <- lapply(
    X = dims,
    FUN = function(x) {
      cells <- TopCells(
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
    FUN = TopFeatures,
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
  DefaultAssay(object = object) <- assays
  
  #get assay data with only selected features (all dims) and selected cells(all)
  data.all <- FetchData(
    object = object,
    vars = features.keyed,
    cells = unique(x = unlist(x = cells)),
    slot = slot
  )
  
  #clip off values for heatmap
  data.all <- MinMax(data = data.all, min = disp.min, max = disp.max)
  data.limits <- c(min(data.all), max(data.all))
  
  #draw heatmap for each dim
  for (i in 1:length(x = dims)) {
    dim.features <- c(features[[i]][[2]], rev(x = features[[i]][[1]]))
    dim.features <- rev(x = unlist(x = lapply(
      X = dim.features,
      FUN = function(feat) {
        return(grep(pattern = paste0(feat, '$'), x = features.keyed, value = TRUE))
      }
    )))
    dim.cells <- cells[[i]]
    data.plot <- data.all[dim.cells, dim.features]
    if (fast) {
      SingleImageMap(
        data = data.plot,
        title = paste0(Key(object = object[[reduction]]), dims[i]),
        order = dim.cells
      )
    } else {
      plots[[i]] <- heatmap3(t(as.matrix(data.plot)), 
                             Rowv = NA, 
                             Colv = NA, 
                             scale = "none", 
                             margins = c(3,3), 
                             balanceColor = TRUE)
    }
  }
  return(plots[[1]])
}