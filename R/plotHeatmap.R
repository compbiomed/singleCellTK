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
#' @param externalReduction
#'
#' @return
#' @export
plotHeatmap <- function(inSCE,
                        useAssay,
                        dims = 1:dims,
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
                        combine = TRUE,
                        externalReduction = NULL){
  object <- convertSCEToSeurat(inSCE, scaledAssay = useAssay)
  
  if(!is.null(externalReduction)){
    if(reduction == "pca"){
      object@reductions <- list(pca = externalReduction)
    }
    else{
      object@reductions <- list(ica = externalReduction)
    }
  }
  
  ncol <- ncol %||% ifelse(test = length(x = dims) > 2, yes = 3, no = length(x = dims))
  
  #empty list = number of dims
  plots <- vector(mode = 'list', length = length(x = dims))
  
  #set rna
  assays <- assays %||% DefaultAssay(object = object)
  
  #set disp.max
  disp.max <- disp.max %||% ifelse(
    test = slot == 'scale.data',
    yes = 2.5,
    no = 6
  )
  
  #no. of cells
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
  
  #convert (_) to (-) as required by FetchData function below
  cells <- lapply(
    X = cells, 
    FUN = function(t) gsub(
      pattern = "_", 
      replacement = "-", 
      x = t, 
      fixed = TRUE)
    )
  
  #get assay data with only selected features (all dims) and selected cells (all)
  data.all <- FetchData(
    object = object,
    vars = features.keyed,
    cells = unique(x = unlist(x = cells)),
    slot = slot
  )
  
  #
  
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
    # if (fast) {
    #   #need new ggplot fast method here = todo
    #   
    #   # SingleImageMap(
    #   #   data = data.plot,
    #   #   title = paste0(Key(object = object[[reduction]]), dims[i]),
    #   #   order = dim.cells
    #   # )
    #   
    # } else {
    #   data <- data.plot
    #   data$samples <- rownames(data)
    # 
    #   #data for only for ggplot convert to tibble
    #   exp.long <- pivot_longer(data = data,
    #                            cols = -c(samples),
    #                            names_to = "gene",
    #                            values_to = "expression")
    #   
    #   #to reorder samples and features
    #   exp.long$samples <- factor(x = exp.long$samples,
    #                                  levels = dim.cells)
    #   exp.long$gene <- factor(x = exp.long$gene,
    #                               levels = dim.features)
    #   #plot
    #   exp.heatmap <- ggplot(data = exp.long, mapping = aes(x = samples,
    #                                                        y = gene,
    #                                                        fill = expression)) +
    #     geom_tile() +
    #     xlab(label = "samples") + 
    #     theme(axis.title.y = element_blank(), 
    #           axis.text.x = element_text(angle = 90, size = 1))
    
    hm <- ComplexHeatmap::Heatmap(t(data.plot), 
            show_row_dend = FALSE, 
            show_column_dend = FALSE, 
            cluster_rows = FALSE, 
            cluster_columns = FALSE, 
            show_column_names = FALSE, 
            row_names_side = "left")
      
      plots[[i]] <- hm
    #}
  }
  
  ht_list = NULL
  for(i in seq(plots)) {
    ht_list = ht_list + plots[[i]]
  }
  
  # nCol <- floor(sqrt(length(plots)))
  # plot <- do.call("grid.arrange", c(plots, ncol=nCol))
  # plot <- as_ggplot(plot)
  ht_list <- ComplexHeatmap::draw(ht_list, auto_adjust = FALSE)
  return(ht_list)
}