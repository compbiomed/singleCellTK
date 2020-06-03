#' @title Get runDropletQC .html report 
#' @description A  function to generate .html Rmarkdown report containing the visualizations of the runDropletQC function output
#' @param inSCE A \link[SingleCellExperiment]{SingleCellExperiment} object containing
#' the full droplet count matrix with the output from runDropletQC function
#' @param output_file name of the generated file. If NULL/default then the output file name will be based on the name of the Rmarkdown template 
#' @param output_dir name of the output directory to save the rendered file. If NULL/default the file is stored to the current working directory
#' @return .html file
#' @examples
#' inSCE <- runDropletQC(inSCE)
#' getHtmlReportDropletQC(inSCE = inSCE)
#' @export
getHtmlReportDropletQC <- function(inSCE, output_file = NULL,
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
#' inSCE <- runCellQC(inSCE)
#' getHtmlReportCellQC(inSCE = inSCE)
#' @export
getHtmlReportCellQC <- function(inSCE, output_file = NULL,
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
#' inSCE <- runDecontX(inSCE)
#' getHtmlReportPerQC(inSCE = outSCE, QCtype = "DecontX")
#' @export
getHtmlReportPerQC <- function(inSCE, QCtype="DecontX", output_file = NULL,
                            output_dir = NULL) {
  
  if (is.null(output_dir)){
    output_dir<- getwd()
  }
  
  if (QCtype =="BarcodeRankDrops"){
    rmarkdown::render(system.file("rmarkdown/qc/BarcodeRankDrops.Rmd", package="singleCellTK"), params = list(
      object=inSCE), output_file = output_file, output_dir = output_dir)
  }
  if (QCtype =="EmptyDrops"){
    rmarkdown::render(system.file("rmarkdown/qc/EmptyDrops.Rmd", package="singleCellTK"), params = list(
      object=inSCE), output_file = output_file, output_dir = output_dir)
  }
  if (QCtype =="Cxds"){
    rmarkdown::render(system.file("rmarkdown/qc/Cxds.Rmd", package="singleCellTK"), params = list(
      object=inSCE), output_file = output_file, output_dir = output_dir)
  }
  if (QCtype =="Bcds"){
    rmarkdown::render(system.file("rmarkdown/qc/Bcds.Rmd", package="singleCellTK"), params = list(
      object=inSCE), output_file = output_file, output_dir = output_dir)
  }
  if (QCtype =="CxdsBcdsHybrid"){
    rmarkdown::render(system.file("rmarkdown/qc/CxdsBcdsHybrid.Rmd", package="singleCellTK"), params = list(
      object=inSCE), output_file = output_file, output_dir = output_dir)
  }
  if (QCtype =="DecontX"){
    rmarkdown::render(system.file("rmarkdown/qc/DecontX.Rmd", package="singleCellTK"), params = list(
      object=inSCE), output_file = output_file, output_dir = output_dir)
  }
  if (QCtype =="DoubletCells"){
    rmarkdown::render(system.file("rmarkdown/qc/DoubletCells.Rmd", package="singleCellTK"), params = list(
      object=inSCE), output_file = output_file, output_dir = output_dir)
  }
  if (QCtype =="QCMetrics"){
    rmarkdown::render(system.file("rmarkdown/qc/QCMetrics.Rmd", package="singleCellTK"), params = list(
      object=inSCE), output_file = output_file, output_dir = output_dir)
  }
  if (QCtype =="Scrublet"){
    rmarkdown::render(system.file("rmarkdown/qc/Scrublet.Rmd", package="singleCellTK"), params = list(
      object=inSCE), output_file = output_file, output_dir = output_dir)
  }
  if (QCtype =="DoubletFinder"){
    rmarkdown::render(system.file("rmarkdown/qc/DoubletFinder.Rmd", package="singleCellTK"), params = list(
      object=inSCE), output_file = output_file, output_dir = output_dir)
  }
 }


    
