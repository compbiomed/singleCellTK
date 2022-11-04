#' @title Plotting function of runMusic.R
#' @description A wrapper that plots heatmap and cluster plots for results from runMusic.R
#' @param inSCE A \link[SingleCellExperiment]{SingleCellExperiment} object.
#' @param analysisType Character. Specify which function to run 
#'  Available options are  "EstCellProp","PreGroupedClustProp","SingleCellClust"
#' @param heatmapTitle Character. Title for heatmap; Default is NULL
#' @param analysisName Character. User-defined analysis name.
#' This will be used as the slot name and results can be stored and retrived from SCE object using this name
#' @return SingleCellExperiment object containing the outputs of the
#'  specified algorithms in the \link{colData}
#' of \code{inSCE}.

library(singleCellTK)


plotMusicResults<- function(inSCE, 
                            analysisType = c("EstCellProp","PreGroupedClustProp","SingleCellClust"),
                            heatmapTitle = NULL,
                            analysisName = NULL,
                            useAssay = NULL){
  
  
  # Plot clusters
  
  .musicPlotclusters<- function(inSCE,analysisName){
    basisObject = metadata(inSCE)[["sctk"]][["music"]][[analysisName]][["Disgn.mtx"]]
    calc_dmat_dist = dist(t(log(basisObject + 1e-6)), method = "euclidean")
    hc1<- hclust(calc_dmat_dist, method = "complete")
    basisObject2 =  metadata(inSCE)[["sctk"]][["music"]][[analysisName]][["M.theta"]]
    calc_RA_dist = dist(t(log(basisObject2 + 1e-8)), method = "euclidean")
    hc2<- hclust(calc_RA_dist,method = "complete")
    p1<-plot(hc1,cex =0.6,hang = -1, main = 'Cluster log(Design Matrix)')
    p2<-plot(hc2,cex =0.6,hang = -1, main = 'Cluster log(Mean of RA)')
    return(list(hc1,hc2))
    
  }
  
  
  # plot heatmap  # use the ones from lab inSCEheatmap? (uses complexheatmap)
  
  .plotHeatmap<- function(inSCE, # 
                          analysisName,
                          heatmapTitle = NULL, 
                          useAssay = useAssay){
    
    
    
    
    testBulk<-metadata(inSCE)$sctk$music[[analysisName]][[useAssay]]
    bulkinSCE<-SingleCellExperiment(assays = list(EstProps = testBulk))
    names(assays(bulkinSCE))<-useAssay
    heatmap<-plotSCEHeatmap(bulkinSCE,  # Check version
                            useAssay = useAssay, 
                            rowLabel = T, 
                            colLabel = T, 
                            rowTitle = "Subjects",
                            colTitle = "CellType",
                            title = heatmapTitle)
    return(heatmap)
  }
  
  # sce -> extract bulk -> new_sce -> plotsce
  
  
  
  if(analysisType == "SingleCellClust"){
    plots<-.musicPlotclusters(inSCE,analysisName = analysisName)
  #  metadata(inSCE)$sctk$music[[analysisName]][["Clusters"]]<- temp_results 
  }
  else if(analysisType == "EstCellProp"){

    plots<-.plotHeatmap(inSCE,heatmapTitle = heatmapTitle,useAssay = "Est.prop.weighted" ,analysisName = analysisName )
    # Do if else for pulling the assay data and have one single call for the sce heatmap
    
    #metadata(inSCE)$sctk$music[[analysisName]][["Heatmap"]]<- temp_results 
    
  }
    
  else if(analysisType == "PreGroupedClustProp"){
        
    plots<- .plotHeatmap(inSCE,heatmapTitle = heatmapTitle,useAssay = "Est.prop.weighted.cluster",analysisName = analysisName )
    #metadata(inSCE)$sctk$music[[analysisName]][["Heatmap"]]<- temp_results 
    
  }


  return(plots)
} 

