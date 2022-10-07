#' @title Deconvolution of RNASeq data using single cell data
#' @description A wrapper that performs deconvolution and clustering using MuSiC tool and 
#' SingleCellExperiment object
#' @param inSCE A \link[SingleCellExperiment]{SingleCellExperiment} object.
#' @param analysisType Character. Specify which function to run 
#'  Available options are  "EstCellProp","PreGroupedClustProp","SingleCellClust"
#' @param analysisName Character. User-defined analysis name. 
#' This will be used as the slot name and results can be stored and retrived from SCE object using this name
#' @param markers List. list of gene names. Same as group.markers option from MuSiC package. The list include differential expressed genes within groups. 
#' List name must be the same as `clusters`. Default is NULL
#' @param clusters character, the colData of single cell dataset used as clusters; Default is "cellType"
#' @param samples . Default is sampleID.
#' groups = NULL, 
#' @param selectCt vector of cell types, default as NULL. If NULL, then use all cell types provided by single cell dataset; NULL, #same as select.ct
#' @param cellSize 	data.frame of cell sizes.same as cell_size; data.frame of cell sizes. 1st column contains the names of cell types, 2nd column has the cell sizes per cell type. Default as NULL. If NULL, then estimate cell size from data;
#' @param ctCov logical. If TRUE, use the covariance across cell types; #same as ctCov in MuSiC
#' @param preClusterlist 	list of cell types. The list identify groups of similar cell types.
#' @param verbose logical, default as TRUE.
#' @param iter.max 	numeric, maximum iteration number. Default 1000
#' @param nu  regulation parameter, take care of weight when taking reciprocal 1e-04,
#' @param eps Threshold of convergence. Default 0.01,
#' @param centered logic, subtract avg of Y and D. Default FALSE
#' @param normalize logic, divide Y and D by their standard deviation. Default FALSE
#' @return SingleCellExperiment object containing the outputs of the
#'  specified algorithms in the \link{colData}
#' of \code{inSCE}.
#' @examples
#' data(scExample, package = "singleCellTK")
#' Add bulk data here
#' \dontrun{
#' sce <- runMusic(sce,bulkdata, analysisType = "EstCellProp",analysisName = "test")
#' }
#' @export



#' @rdname s4_methods
setGeneric("getMusicResults", signature = c("x","y"),
           function(x,y) {standardGeneric("getMusicResults")}
)

setClass("y", representation(name = "character"))


#' @rdname s4_methods
setMethod("getMusicResults", signature = c(x = "SingleCellExperiment"), function(x,y){
  results <- S4Vectors::metadata(x)$sctk$music[[y]]
  if(is.null(results)) {
    stop("No results from 'Music' are found. Please run `runMusic` first.") 
  }    
  return(results)
})

#' @rdname s4_methods
setGeneric("getMusicResults<-", function(x,y, value) standardGeneric("getMusicResults<-"))

#' @rdname s4_methods
setReplaceMethod("getMusicResults", signature(x = "SingleCellExperiment",y = "character"), function(x,y, value) {
  S4Vectors::metadata(x)$sctk$music[[y]] <- value
  return(x)
})


runMusic<-function(inSCE, 
                   bulkData, #camelcase
                   analysisName = "NULL",
                   analysisType = c("EstCellProp","PreGroupedClustProp","SingleCellClust"), 
                   markers = NULL, 
                   clusters = "cellType", # not a default -- user input 
                   samples = "sampleID", #sample is the default and not sampleID but keeping this as sampleID for testing purpose
                   preClusterlist = NULL,
                   DEmarkers = NULL,
                   groups = NULL, 
                   selectCt = NULL, #same as select.ct
                   cellSize = NULL, #same as cell_size 
                   ctCov = FALSE, #same as ctCov
                   verbose = TRUE,
                   iterMax = 1000, #same as iter.max 
                   nu = 1e-04,
                   eps = 0.01,
                   centered = FALSE,
                   normalize = FALSE,
                   nonZero = TRUE # sane as non.zero
){
  
  
  #####################################################################
  ## Extract the data
  ####################################################################
  # convert inSCE into eset
  
  assay = assay(inSCE)
  pheno = colData(inSCE)
  meta= new('AnnotatedDataFrame', data = data.frame(pheno))
  eset = Biobase::ExpressionSet(assayData = assay, phenoData = meta)
  
  
  # Estimate cell type proportions 
  
  
  .musicProp <-function(bulkData, 
                        eset, 
                        analysisType, 
                        markers, 
                        clusters, 
                        samples, 
                        selectCt, 
                        cellSize,
                        ctCov,
                        iterMax,
                        nu,
                        eps,
                        centered,
                        normalize){
    # Can also supply list of marker genes here as an input
    est_prop =  music_prop(bulk = bulkData,
                           sc.eset = eset,
                           markers = markers, 
                           samples = samples, 
                           clusters = clusters, 
                           select.ct= selectCt, 
                           cell_size = cellSize,
                           ct.cov = ctCov,
                           iter.max = iterMax,
                           nu= nu,
                           eps = eps,
                           centered = centered,
                           normalize = normalize)
    est_prop$analysisType = analysisType
    return(est_prop)
  }
  
  
  
  
  
  ###############################################################################
  ##### Estimation of cell types with pre-grouping of cell types 
  ###############################################################################
  
  
  
  .musicBase<- function(eset,
                        clusters,
                        samples, 
                        markers,
                        selectCt,
                        nonZero,
                        cellSize,
                        ctCov
  ){
    
    basis_object = music_basis(x = eset, 
                               clusters = clusters, 
                               samples = samples,
                               markers = markers,
                               select.ct = selectCt,
                               non.zero = nonZero, 
                               cell_size = cellSize,
                               ct.cov = ctCov,
    )
    basis_object$analysisType = analysisType
    # putting things back to sce
    # should this be a new S4? How to create a new slot?
    return(basis_object)
    
  }
  
  .musicPropCluster<- function(bulkData,
                               eset, 
                               clusters, 
                               groups,
                               preClusterlist, # list of list cluster groups
                               DEmarkers, # group names should be same as cluster names
                               samples,
                               iterMax,
                               nu,
                               eps,
                               centered,
                               normalize){
    
    
    # Preprocess cluster labels
    preClusterlist = clusters.type
    data<-pData(eset)
    clusterExclude = levels(factor(unique(data[[clusters]][data[[clusters]] %in% unlist(preClusterlist) == FALSE])))
    mergeall<-append(preClusterlist,clusterExclude)
    names(mergeall)<-c(names(preClusterlist),clusterExclude)
    cluster_new<-data.frame(do.call(cbind,mergeall)) %>% gather() %>% unique() %>% dplyr::rename(!!clusters:="value", !!groups:="key")
    
    
    # adding cluster labels to phenodata
    data %>% 
      left_join(cluster_new) %>%
      transform(clusterType = as.factor(clusterType)) -> data
    rownames(data)<-rownames(pData(eset))
    pData(eset)<-data
    
    
    prop_clust  = music_prop.cluster(bulk.eset = bulkData, 
                                     sc.eset = eset, 
                                     group.markers = DEmarkers, 
                                     clusters = clusters, 
                                     groups = groups, 
                                     samples = samples, 
                                     clusters.type = preClusterlist,
                                     iter.max = iterMax,
                                     nu = nu,
                                     eps = eps,
                                     centered = centered,
                                     normalize = normalize,)
    
    return(prop_clust)
  }
  
  
  #####################################################################
  ## Run the tool
  #####################################################################
  
  # Estimate cell type proportions
  
  if(analysisType == "EstCellProp"){
    
    temp_result<- .musicProp(bulkData, 
                             eset, 
                             analysisType, 
                             markers, 
                             clusters, 
                             samples, 
                             selectCt, 
                             cellSize,
                             ctCov,
                             iterMax,
                             nu,
                             eps,
                             centered,
                             normalize)
    temp_result$analysisType = analysisType
  }
  
  # Clustering of single cell data 
  
  else if (analysisType == "SingleCellClust"){
    
    temp_result<- .musicBase(eset,
                             clusters,
                             samples, 
                             markers,
                             selectCt,
                             nonZero,
                             cellSize,
                             ctCov)
    temp_result$analysisType = analysisType
    
  }
  
  # Bulk tissue type estimation
  
  else if(analysisType == "PreGroupedClustProp") {
    
    if(class(preClusterlist) == "list"){
      
      temp_result = .musicPropCluster(bulkData = bulkData, 
                                      eset = eset, 
                                      DEmarkers = IEmarkers, 
                                      clusters = clusters, 
                                      groups = groups, 
                                      samples = samples, 
                                      preClusterlist = preClusterlist,
                                      iterMax = iterMax,
                                      nu = nu,
                                      eps = eps,
                                      centered = centered,
                                      normalize = normalize)
      
      temp_result$analysisType = analysisType
      
    }
    
    
    
  }
  
  
  temp_result[["params"]]<-c(as.list(environment()))
  
  
  if(length(inSCE@metadata$sctk$music)>0){
    getMusicResults(x = inSCE, y = analysisName) <- temp_result
   # metadata(inSCE)$sctk$music[[analysisName]]<-temp_result
  }
  else{
    new_list<-c()
  #  metadata(inSCE)$sctk$music<-new_list
   # metadata(inSCE)$sctk$music[[analysisName]]<-temp_result
    getMusicResults(x = inSCE, y = analysisName)<-temp_result
  }
  
  
  
  return(inSCE)
  
}



