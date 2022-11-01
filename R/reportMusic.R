
reportMusic <- function(inSCE, output_file = NULL,
                            output_dir = NULL) {
  
  if (is.null(output_dir)){
    output_dir<- getwd()
  }
  
  rmarkdown::render(system.file("rmarkdown/reportMusicRun.Rmd", package = "singleCellTK"),
                    params = list(object = inSCE, 
                                  bulkData = bulkData,
                                  analysisName = analysisName, 
                                  analysisType = c("EstCellProp","PreGroupedClustProp","SingleCellClust")),
                    output_file = output_file,
                    output_dir = output_dir,
                    intermediates_dir = output_dir,
                    knit_root_dir = output_dir)
  
}

reportMusicresults