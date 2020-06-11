#' @title Get runDropletQC .html report 
#' @description A  function to generate .html Rmarkdown report containing the visualizations of the runDropletQC function output
#' @param inSCE A \link[SingleCellExperiment]{SingleCellExperiment} object containing
#' the full droplet count matrix with the output from runDropletQC function
#' @param output_file name of the generated file. If NULL/default then the output file name will be based on the name of the Rmarkdown template 
#' @param output_dir name of the output directory to save the rendered file. If NULL/default the file is stored to the current working directory
#' @return .html file
#' @examples
#' \donttest{
#' data(emptyDropsSceExample, package = "singleCellTK")
#' sce <- runDropletQC(emptyDropsSceExample)
#' reportDropletQC(inSCE = sce)
#' }
#' @export
reportDropletQC <- function(inSCE, output_file = NULL,
                                   output_dir = NULL) {
  
  if (is.null(output_dir)){
    output_dir<- getwd()
    }
 
  rmarkdown::render(system.file("rmarkdown/qc/DropletQC.Rmd", package="singleCellTK"), params = list(
    object=inSCE), output_file = output_file, output_dir = output_dir )
 }


#' @title Get runCellQC .html report 
#' @description A  function to generate .html Rmarkdown report containing the visualizations of the runCellQC function output
#' @param inSCE A \link[SingleCellExperiment]{SingleCellExperiment} object containing
#' the filtered count matrix with the output from runCellQC function
#' @return .html file
#' @param output_file name of the generated file. If NULL/default then the output file name will be based on the name of the Rmarkdown template. 
#' @param output_dir name of the output directory to save the rendered file. If NULL/default the file is stored to the current working directory
#' @return .html file
#' @examples
#' \donttest{
#' data(sceQCExample, package = "singleCellTK")
#' sceQCExample <- runCellQC(sceQCExample)
#' reportCellQC(inSCE = sceQCExample)
#' }
#' @export
reportCellQC <- function(inSCE, output_file = NULL,
                                output_dir = NULL) {
  if (is.null(output_dir)){
    output_dir<- getwd()
  }
  rmarkdown::render(system.file("rmarkdown/qc/CellQC.Rmd", package="singleCellTK"), params = list(object=inSCE), output_file = output_file, output_dir = output_dir)
}


#' @title Get .html report of the output of the selected QC algorithm
#' @description A  function to generate .html Rmarkdown report for the specified QC algorithm output
#' @param inSCE A \link[SingleCellExperiment]{SingleCellExperiment} object containing
#' the count matrix (full droplets or filtered matrix, depends on the selected QC algorithm) with the output from at least one of these functions:
#' runQCMetrics, runScrublet, runDoubletCells, runCxds, runBcds, runCxdsBcdsHybrid, runDecontX, runBarcodeRankDrops, runEmptyDrops
#' @param algorithm Character. Specifies which QC algorithm report to generate.
#'  Available options are "BarcodeRankDrops", "EmptyDrops", "QCMetrics", "Scrublet", "DoubletCells", "Cxds", "Bcds", "CxdsBcdsHybrid", "DoubletFinder"  and "DecontX".
#' @param output_file name of the generated file. If NULL/default then the output file name will be based on the name of the selected QC algorithm name . 
#' @param output_dir name of the output directory to save the rendered file. If NULL/default the file is stored to the current working directory
#' @return .html file
#' @examples
#' \donttest{
#' data(sceQCExample, package = "singleCellTK")
#' sceQCExample <- runDecontX(sceQCExample)
#' reportQCTool(inSCE = sceQCExample, algorithm = "DecontX")
#' }
#' @export
reportQCTool <- function(inSCE, algorithm=c("BarcodeRankDrops",
                                            "EmptyDrops",
                                            "QCMetrics",
                                            "Scrublet",
                                            "DoubletCells",
                                            "Cxds",
                                            "Bcds",
                                            "CxdsBcdsHybrid",
                                            "DoubletFinder",
                                            "DecontX"),
                         output_file = NULL,
                            output_dir = NULL) {
  
  algorithm <- match.arg(algorithm)
  
  if (is.null(output_dir)){
    output_dir<- getwd()
  }
  
  if (algorithm =="BarcodeRankDrops"){
    rmarkdown::render(system.file("rmarkdown/qc/BarcodeRankDrops.Rmd", package="singleCellTK"), params = list(
      object=inSCE), output_file = output_file, output_dir = output_dir)
  }
  if (algorithm =="EmptyDrops"){
    rmarkdown::render(system.file("rmarkdown/qc/EmptyDrops.Rmd", package="singleCellTK"), params = list(
      object=inSCE), output_file = output_file, output_dir = output_dir)
  }
  if (algorithm =="Cxds"){
    rmarkdown::render(system.file("rmarkdown/qc/Cxds.Rmd", package="singleCellTK"), params = list(
      object=inSCE), output_file = output_file, output_dir = output_dir)
  }
  if (algorithm =="Bcds"){
    rmarkdown::render(system.file("rmarkdown/qc/Bcds.Rmd", package="singleCellTK"), params = list(
      object=inSCE), output_file = output_file, output_dir = output_dir)
  }
  if (algorithm =="CxdsBcdsHybrid"){
    rmarkdown::render(system.file("rmarkdown/qc/CxdsBcdsHybrid.Rmd", package="singleCellTK"), params = list(
      object=inSCE), output_file = output_file, output_dir = output_dir)
  }
  if (algorithm =="DecontX"){
    rmarkdown::render(system.file("rmarkdown/qc/DecontX.Rmd", package="singleCellTK"), params = list(
      object=inSCE), output_file = output_file, output_dir = output_dir)
  }
  if (algorithm =="DoubletCells"){
    rmarkdown::render(system.file("rmarkdown/qc/DoubletCells.Rmd", package="singleCellTK"), params = list(
      object=inSCE), output_file = output_file, output_dir = output_dir)
  }
  if (algorithm =="QCMetrics"){
    rmarkdown::render(system.file("rmarkdown/qc/QCMetrics.Rmd", package="singleCellTK"), params = list(
      object=inSCE), output_file = output_file, output_dir = output_dir)
  }
  if (algorithm =="Scrublet"){
    rmarkdown::render(system.file("rmarkdown/qc/Scrublet.Rmd", package="singleCellTK"), params = list(
      object=inSCE), output_file = output_file, output_dir = output_dir)
  }
  if (algorithm =="DoubletFinder"){
    rmarkdown::render(system.file("rmarkdown/qc/DoubletFinder.Rmd", package="singleCellTK"), params = list(
      object=inSCE), output_file = output_file, output_dir = output_dir)
  }
 }


    
