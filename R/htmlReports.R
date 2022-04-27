#' @title Get runDropletQC .html report
#' @description A  function to generate .html Rmarkdown report containing the visualizations of the runDropletQC function output
#' @param inSCE A \link[SingleCellExperiment]{SingleCellExperiment} object containing
#' the full droplet count matrix with the output from runDropletQC function
#' @param subTitle subtitle of the QC HTML report. Default is NULL.
#' @param studyDesign description of the data set and experiment design. It would be shown at the top of QC HTML report. Default is NULL.
#' @param output_file name of the generated file. If NULL/default then the output file name will be based on the name of the Rmarkdown template
#' @param output_dir name of the output directory to save the rendered file. If NULL/default the file is stored to the current working directory
#' @return .html file
#' @examples
#' data(scExample, package = "singleCellTK")
#' \dontrun{
#' sce <- runDropletQC(sce)
#' reportDropletQC(inSCE = sce)
#' }
#' @export
reportDropletQC <- function(inSCE, output_file = NULL,
                                   output_dir = NULL,
                                   subTitle = NULL,
                                   studyDesign = NULL) {

  if (is.null(output_dir)){
    output_dir<- getwd()
    }

  #report_path <- tempfile(fileext = ".Rmd")
  #file.copy(system.file("rmarkdown/qc/DropletQC.Rmd", package = "singleCellTK"), report_path, overwrite = TRUE)

  ## create temp Rmd file to bypass permission issue on server
  rmarkdown::render(system.file("rmarkdown/qc/DropletQC.Rmd", package = "singleCellTK"),
    params = list(object = inSCE, subTitle = subTitle, studyDesign = studyDesign),
    output_file = output_file,
    output_dir = output_dir,
    intermediates_dir = output_dir,
    knit_root_dir = output_dir)
 }


#' @title Get runCellQC .html report
#' @description A  function to generate .html Rmarkdown report containing the visualizations of the runCellQC function output
#' @param inSCE A \link[SingleCellExperiment]{SingleCellExperiment} object containing
#' the filtered count matrix with the output from runCellQC function
#' @param useReducedDim Character. The name of the saved dimension reduction slot including cells  
#' from all samples in then\linkS4class{SingleCellExperiment} object, Default is NULL
#' @param subTitle subtitle of the QC HTML report. Default is NULL.
#' @param studyDesign Character. The description of the data set and experiment design. It would be shown at the top of QC HTML report. Default is NULL.
#' @param output_file Character. The name of the generated file. If NULL/default then the output file name will be based on the name of the Rmarkdown template.
#' @param output_dir Character. The name of the output directory to save the rendered file. If NULL/default the file is stored to the current working directory
#' @return .html file
#' @examples
#' data(scExample, package = "singleCellTK")
#' sce <- subsetSCECols(sce, colData = "type != 'EmptyDroplet'")
#' \dontrun{
#' sce <- runCellQC(sce)
#' reportCellQC(inSCE = sce)
#' }
#' @export
reportCellQC <- function(inSCE, output_file = NULL,
                                output_dir = NULL,
                                subTitle = NULL,
                                studyDesign = NULL,
                                useReducedDim = NULL) {
  if (is.null(output_dir)){
    output_dir<- getwd()
  }
  #report_path <- tempfile(fileext = ".Rmd")
  #file.copy(system.file("rmarkdown/qc/CellQC.Rmd", package = "singleCellTK"), report_path, overwrite = TRUE)

  rmarkdown::render(system.file("rmarkdown/qc/CellQC.Rmd", package = "singleCellTK"),
    params = list(object = inSCE, subTitle = subTitle, studyDesign = studyDesign,
    reducedDimName = useReducedDim),
    output_file = output_file,
    output_dir = output_dir,
    intermediates_dir = output_dir,
    knit_root_dir = output_dir)
}


#' @title Get .html report of the output of the selected QC algorithm
#' @description A  function to generate .html Rmarkdown report for the specified QC algorithm output
#' @param inSCE A \link[SingleCellExperiment]{SingleCellExperiment} object containing
#' the count matrix (full droplets or filtered matrix, depends on the selected QC algorithm) with the output from at least one of these functions:
#' runQCMetrics, runScrublet, runScDblFinder, runCxds, runBcds, runCxdsBcdsHybrid, runDecontX, runBarcodeRankDrops, runEmptyDrops
#' @param algorithm Character. Specifies which QC algorithm report to generate.
#'  Available options are "BarcodeRankDrops", "EmptyDrops", "QCMetrics", "Scrublet", "ScDblFinder", "Cxds", "Bcds", "CxdsBcdsHybrid", "DoubletFinder", "DecontX" and "SoupX".
#' @param output_file name of the generated file. If NULL/default then the output file name will be based on the name of the selected QC algorithm name .
#' @param output_dir name of the output directory to save the rendered file. If NULL/default the file is stored to the current working directory
#' @return .html file
#' @examples
#' data(scExample, package = "singleCellTK")
#' sce <- subsetSCECols(sce, colData = "type != 'EmptyDroplet'")
#' \dontrun{
#' sce <- runDecontX(sce)
#' sce <- getUMAP(sce)
#' reportQCTool(inSCE = sce, algorithm = "DecontX")
#' }
#' @export
reportQCTool <- function(inSCE, algorithm=c("BarcodeRankDrops",
                                            "EmptyDrops",
                                            "QCMetrics",
                                            "Scrublet",
                                            "ScDblFinder",
                                            "Cxds",
                                            "Bcds",
                                            "CxdsBcdsHybrid",
                                            "DoubletFinder",
                                            "DecontX",
                                            "SoupX"),
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
  if (algorithm =="SoupX"){
    rmarkdown::render(system.file("rmarkdown/qc/SoupX.Rmd", package="singleCellTK"), params = list(
      object=inSCE), output_file = output_file, output_dir = output_dir)
  }
  if (algorithm =="ScDblFinder"){
    rmarkdown::render(system.file("rmarkdown/qc/ScDblFinder.Rmd", package="singleCellTK"), params = list(
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

#' @title Get runDEAnalysis .html report
#' @description A  function to generate .html Rmarkdown report containing the
#' visualizations of the \code{\link{runDEAnalysis}} function output
#' @param inSCE A \code{\link[SingleCellExperiment]{SingleCellExperiment}}
#' object containing the output from \code{\link{runDEAnalysis}} function
#' @param study The specific analysis to visualize, used as \code{analysisName}
#' argument when running differential expression.
#' @param output_file name of the generated file. If \code{NULL} then the output
#' file name will be based on the name of the Rmarkdown template. Default
#' \code{NULL}.
#' @param output_dir name of the output directory to save the rendered file. If
#' \code{NULL} the file is stored to the current working directory.
#' Default \code{NULL}.
#' @return Saves the HTML report in the specified output directory.
#' @export
reportDiffExp <- function(inSCE, study,
                          output_file = NULL,
                          output_dir = NULL) {

  if (is.null(output_dir)){
    output_dir <- getwd()
  }
  if (!study %in% names(S4Vectors::metadata(inSCE)$diffExp)) {
    stop("Specified study not found in given SCE object")
  }
  rmarkdown::render(system.file("rmarkdown/de/DifferentialExpression.Rmd",
                                package="singleCellTK"),
                    params = list(object=inSCE, study=study),
                    output_file = output_file,
                    output_dir = output_dir )
}

#' @title Get findMarkerDiffExp .html report
#' @description A  function to generate .html Rmarkdown report containing the
#' visualizations of the \code{\link{findMarkerDiffExp}} function output
#' @param inSCE A \code{\link[SingleCellExperiment]{SingleCellExperiment}}
#' object containing the output from \code{\link{findMarkerDiffExp}} function
#' @param output_file name of the generated file. If \code{NULL} then the output
#' file name will be based on the name of the Rmarkdown template. Default
#' \code{NULL}.
#' @param output_dir name of the output directory to save the rendered file. If
#' \code{NULL} the file is stored to the current working directory.
#' Default \code{NULL}.
#' @return An HTML file of the report will be generated at the path specified
#' in the arguments.
#' @export
reportFindMarker <- function(inSCE, output_file = NULL, output_dir = NULL) {

  if (is.null(output_dir)){
    output_dir <- getwd()
  }
  if (!"findMarker" %in% names(S4Vectors::metadata(inSCE))) {
    stop("Find marker result not presented in input SCE object. Run ",
         "findMarkerDiffExp() first. ")
  }
  att <- names(attributes(S4Vectors::metadata(inSCE)$findMarker))
  if (!"useAssay" %in% att) {
    stop("Can't identify the structure of find marker result. Run ",
         "findMarkerDiffExp() first. ")
  }
  rmarkdown::render(system.file("rmarkdown/de/FindMarker.Rmd",
                                package="singleCellTK"),
                    params = list(object = inSCE),
                    output_file = output_file,
                    output_dir = output_dir )
}


#' Generates an HTML report for Seurat Run (including Normalization, 
#'  Feature Selection, Dimensionality Reduction & Clustering) and returns the 
#'  SCE object with the results computed and stored inside the object.
#' @param inSCE Input \code{\link[SingleCellExperiment]{SingleCellExperiment}}
#'  object.
#' @param biological.group A character value that specifies the name of the 
#'  \code{colData()} column to use as the main biological group in the Seurat 
#'  report for tSNE & UMAP visualization.
#' @param phenotype.groups A character value that specifies the name of the 
#'  \code{colData()} column to use as additional phenotype variables in the 
#'  Seurat report for tSNE & UMAP visualization.
#' @param variable.features A numeric value indicating the number of top 
#'  variable genes to identify in the report. Default is \code{2000}.
#' @param pc.count 	A numeric value indicating the number of principal 
#'  components to use in the analysis workflow. Default is \code{50}.
#' @param runHVG A logical value indicating if feature selection should be run
#'  in the report. Default \code{TRUE}.
#' @param plotHVG A logical value indicating if the top variable genes should
#'  be visualized through a mean-to-variance plot. Default is \code{TRUE}.
#' @param runDimRed A logical value indicating if PCA should be computed in the
#'  report. Default is \code{TRUE}.
#' @param plotJackStraw A logical value indicating if the JackStraw plot should
#'  be visualized for the principal components. Default is \code{FALSE}.
#' @param plotElbowPlot A logical value indicating if the ElbowPlot should be
#'  visualized for the principal components. Default is \code{FALSE}.
#' @param plotHeatmaps A logical value indicating if the Heatmaps should be 
#'  visualized for the principal components. Default is \code{FALSE}.
#' @param runClustering A logical value indicating if Clustering should be
#'  run over multiple resolutions as defined by the \code{minResolution} and
#'  \code{maxResolution} parameters. Default is \code{TRUE}.
#' @param plotTSNE A logical value indicating if TSNE plot should be visualized
#'  for clusters. Default is \code{TRUE}.
#' @param plotUMAP A logical value indicating if UMAP plot should be visualized
#'  for clusters. Default is \code{TRUE}.
#' @param minResolution A numeric value indicating the minimum resolution to use
#'  for clustering. Default \code{0.3}.
#' @param maxResolution A numeric value indicating the maximum resolution to use
#'  for clustering. Default \code{1.5}.
#' @param outputFile Specify the name of the generated output HTML file. 
#'  If \code{NULL} then the output file name will be based on the name of the 
#'  Rmarkdown template. Default \code{NULL}.
#' @param outputDir Specify the name of the output directory to save the 
#'  rendered HTML file. If \code{NULL} the file is stored to the current 
#'  working directory. Default \code{NULL}.
#' @param subtitle A character value specifying the subtitle to use in the 
#'  report. Default \code{NULL}.
#' @param authors A character value specifying the names of the authors to use 
#'  in the report. Default \code{NULL}.
#' @param showSession 	A logical value indicating if session information 
#'  should be displayed or not. Default is \code{FALSE}.
#' @param pdf A logical value indicating if a pdf should also be generated for 
#'  each figure in the report. Default is \code{FALSE}.
#' @param forceRun A logical value indicating if all computations previously
#'  computed should be re-calculated regardless if these computations are
#'  available in the input object. Default is \code{FALSE}.
#'
#' @return A \code{\link[SingleCellExperiment]{SingleCellExperiment}} object
#'  with computations stored.
#' @export
reportSeuratRun <- function(inSCE,
                            biological.group = NULL,
                            phenotype.groups = NULL,
                            variable.features = 2000,
                            pc.count = 50,
                            runHVG = TRUE,
                            plotHVG = TRUE,
                            runDimRed = TRUE,
                            plotJackStraw = FALSE,
                            plotElbowPlot = TRUE,
                            plotHeatmaps = TRUE,
                            runClustering = TRUE,
                            plotTSNE = TRUE,
                            plotUMAP = TRUE,
                            minResolution = 0.3,
                            maxResolution = 1.5,
                            outputFile = NULL,
                            outputDir = NULL,
                            subtitle = NULL,
                            authors =  NULL,
                            showSession = FALSE,
                            pdf = FALSE,
                            forceRun = FALSE){

  if(is.null(biological.group)){
    stop("Must specify atleast one biological.group that is present in the colData of input object.")
  }

  if(!biological.group %in% names(colData(inSCE))){
    stop(biological.group, " not found in the colData of input object.")
  }

  if(!is.null(phenotype.groups)){
    if(!all(phenotype.groups %in% names(colData(inSCE)))){
      stop(phenotype.groups, " not found in the colData of input object.")
    }
  }

  if(is.null(outputDir)){
    outputDir <- getwd()
    message("No output directory defined, using current working directory ", outputDir, " instead.")
  }

  data <- inSCE
  
  rmarkdown::render(system.file("rmarkdown/seurat/reportSeuratRun.Rmd",
                                package="singleCellTK"),
                    params = list(
                      subtitle = subtitle,
                      authors = authors,
                      sce = data,
                      biological.group = biological.group,
                      phenotype.groups = phenotype.groups,
                      variable.features = variable.features,
                      pc.count = pc.count,
                      outputPath = outputDir,
                      showSession = showSession,
                      pdf = pdf,
                      runHVG = runHVG,
                      plotHVG = plotHVG,
                      runDimRed = runDimRed,
                      plotJackStraw = plotJackStraw,
                      plotElbowPlot = plotElbowPlot,
                      plotHeatmaps = plotHeatmaps,
                      runClustering = runClustering,
                      plotTSNE = plotTSNE,
                      plotUMAP = plotUMAP,
                      minResolution = minResolution,
                      maxResolution = maxResolution,
                      forceRun = forceRun
                    ),
                    output_file = outputFile,
                    output_dir = outputDir,
                    intermediates_dir = getwd(),
                    knit_root_dir = getwd())
  
  path <- paste0(outputDir, "SCE_SeuratRun", "-", gsub(" ", "_", Sys.Date()), ".rds")
  saveRDS(data, path)
  message("Output SCE object stored as ", paste0("SCE_SeuratRun", "-", gsub(" ", "_", Sys.Date()), ".rds"), " in ", outputDir, ".")
  message("Output HTML file stored as ", outputFile, " in ", outputDir, ".")

  return(data)
}


#' Generates an HTML report for Seurat Results (including Clustering & Marker
#'  Selection) and returns the SCE object with the results computed and stored
#'  inside the object.
#' @param inSCE Input \code{\link[SingleCellExperiment]{SingleCellExperiment}}
#'  object previously passed through \code{reportSeuratRun()}. 
#' @param biological.group A character value that specifies the name of the 
#'  \code{colData()} column to use as the main biological group in the Seurat 
#'  report for marker selection and grouping.
#' @param phenotype.groups A character vector that specifies the names of the 
#'  \code{colData()} columns to use for differential expression in addition to 
#'  the \code{biological.group} parameter.
#' @param selected.markers A character vector containing the user-specified 
#'  gene symbols or feature names of marker genes that be used to generate 
#'  gene plots in addition to the gene markers computed from 
#'  differential expression.
#' @param clustering.resolution A numeric value indicating the user-specified 
#'  final resolution to use with clustering. Default is \code{0.8}.
#' @param pc.count A numeric value indicating the number of principal components
#'  to use in the analysis workflow. Default is \code{50}.
#' @param plotTSNE A logical value indicating if TSNE plots should be visualized
#'  in the clustering section of the report. Default is \code{TRUE}.
#' @param plotUMAP A logical value indicating if UMAP plots should be visualized
#'  in the clustering section of the report. Default is \code{TRUE}.
#' @param runClustering A logical value indicating if Clustering should be run
#'  or not in the report. Default is \code{TRUE}. If \code{FALSE}, parameters
#'   \code{plotTSNE} and \code{plotUMAP} are also set to \code{FALSE}.
#' @param runMSClusters A logical value indicating if the marker selection
#'  section for identifying marker genes between clusters should be run and
#'  visualized in the report. Default \code{TRUE}.
#' @param runMSBioGroup A logical value indicating if the marker selection
#'  section for identifying marker genes between the \code{biological.group} 
#'  parameter should be run and visualized in the report. Default \code{TRUE}.
#' @param outputFile Specify the name of the generated output HTML file. 
#'  If \code{NULL} then the output file name will be based on the name of the 
#'  Rmarkdown template. Default \code{NULL}.
#' @param outputDir Specify the name of the output directory to save the 
#'  rendered HTML file. If \code{NULL} the file is stored to the current 
#'  working directory. Default \code{NULL}.
#' @param subtitle A character value specifying the subtitle to use in the 
#'  report. Default \code{NULL}.
#' @param authors A character value specifying the names of the authors to use 
#'  in the report. Default \code{NULL}.
#' @param showSession 	A logical value indicating if session information 
#'  should be displayed or not. Default is \code{FALSE}.
#' @param pdf A logical value indicating if a pdf should also be generated for 
#'  each figure in the report. Default is \code{FALSE}.
#' @param forceRun A logical value indicating if all computations previously
#'  computed should be re-calculated regardless if these computations are
#'  available in the input object. Default is \code{FALSE}.
#'
#' @return A \code{\link[SingleCellExperiment]{SingleCellExperiment}} object
#'  with computations stored.
#' @export
reportSeuratResults <- function(inSCE,
                                biological.group = NULL,
                                phenotype.groups = NULL,
                                selected.markers = NULL,
                                clustering.resolution = 0.8,
                                pc.count = 50,
                                plotTSNE = TRUE,
                                plotUMAP = TRUE,
                                runClustering = TRUE,
                                runMSClusters = TRUE,
                                runMSBioGroup = TRUE,
                                outputFile = NULL,
                                outputDir = NULL,
                                subtitle = NULL,
                                authors =  NULL,
                                showSession = FALSE,
                                pdf = FALSE,
                                forceRun = FALSE){
  
  if(is.null(biological.group)){
    stop("Must specify atleast one biological.group that is present in the colData of input object.")
  }
  
  if(!biological.group %in% names(colData(inSCE))){
    stop(biological.group, " not found in the colData of input object.")
  }
  
  if(!is.null(phenotype.groups)){
    if(!all(phenotype.groups %in% names(colData(inSCE)))){
      stop(phenotype.groups, " not found in the colData of input object.")
    }
  }
  
  if(is.null(outputDir)){
    outputDir <- getwd()
    message("No output directory defined, using current working directory ", outputDir, " instead.")
  }
  
  data <- inSCE
  
  rmarkdown::render(system.file("rmarkdown/seurat/reportSeuratResults.Rmd",
                                package="singleCellTK"),
                    params = list(
                      subtitle = subtitle,
                      authors = authors,
                      sce = data,
                      biological.group = biological.group,
                      phenotype.groups = phenotype.groups,
                      selected.markers = selected.markers,
                      clustering.resolution = clustering.resolution,
                      pc.count = pc.count,
                      outputPath = outputDir,
                      showSession = showSession,
                      pdf = pdf,
                      plotTSNE = plotTSNE,
                      plotUMAP = plotUMAP,
                      runClustering = runClustering,
                      runMSClusters = runMSClusters,
                      runMSBioGroup = runMSBioGroup,
                      forceRun = forceRun
                    ),
                    output_file = outputFile,
                    output_dir = outputDir,
                    intermediates_dir = getwd(),
                    knit_root_dir = getwd())
  
  path <- paste0(outputDir, "SCE_SeuratResults", "-", gsub(" ", "_", Sys.Date()), ".rds")
  saveRDS(data, path)
  message("Output SCE object stored as ", paste0("SCE_SeuratResults", "-", gsub(" ", "_", Sys.Date()), ".rds"), " in ", outputDir, ".")
  message("Output HTML file stored as ", outputFile, " in ", outputDir, ".")
  
  return(data)
}


#' Generates an HTML report for Seurat Dimensionality Reduction
#'  and returns the SCE object with the results computed and stored
#'  inside the object.
#' @param inSCE Input \code{\link[SingleCellExperiment]{SingleCellExperiment}}
#'  object.
#' @param pc.count A numeric value indicating the number of principal components
#'  to compute. Default is \code{50}.
#' @param runDimRed A logical value indicating if dimenionality reduction should
#'  be computed. Default \code{TRUE}.
#' @param plotJackStraw A logical value indicating if JackStraw plot should be
#'  visualized. Default \code{FALSE}.
#' @param plotElbowPlot A logical value indicating if ElbowPlot should be 
#'  visualized. Default \code{TRUE}.
#' @param plotHeatmaps A logical value indicating if heatmaps should be
#'  visualized. Default \code{TRUE}.
#' @param outputFile Specify the name of the generated output HTML file. 
#'  If \code{NULL} then the output file name will be based on the name of the 
#'  Rmarkdown template. Default \code{NULL}.
#' @param outputDir Specify the name of the output directory to save the 
#'  rendered HTML file. If \code{NULL} the file is stored to the current 
#'  working directory. Default \code{NULL}.
#' @param subtitle A character value specifying the subtitle to use in the 
#'  report. Default \code{NULL}.
#' @param authors A character value specifying the names of the authors to use 
#'  in the report. Default \code{NULL}.
#' @param showSession 	A logical value indicating if session information 
#'  should be displayed or not. Default is \code{FALSE}.
#' @param pdf A logical value indicating if a pdf should also be generated for 
#'  each figure in the report. Default is \code{FALSE}.
#' @param forceRun A logical value indicating if all computations previously
#'  computed should be re-calculated regardless if these computations are
#'  available in the input object. Default is \code{FALSE}.
#'
#' @return A \code{\link[SingleCellExperiment]{SingleCellExperiment}} object
#'  with computations stored.
#' @export
reportSeuratDimRed <- function(inSCE,
                               pc.count = 50,
                               runDimRed = TRUE,
                               plotJackStraw = FALSE,
                               plotElbowPlot = TRUE,
                               plotHeatmaps = TRUE,
                               outputFile = NULL,
                               outputDir = NULL,
                               subtitle = NULL,
                               authors =  NULL,
                               showSession = FALSE,
                               pdf = FALSE,
                               forceRun = FALSE){

  if(is.null(outputDir)){
    outputDir <- getwd()
    message("No output directory defined, using current working directory ", outputDir, " instead.")
  }
  
  data <- inSCE
  
  rmarkdown::render(system.file("rmarkdown/seurat/reportSeuratDimRed.Rmd",
                                package="singleCellTK"),
                    params = list(
                      subtitle = subtitle,
                      authors = authors,
                      sce = data,
                      pc.count = pc.count,
                      outputPath = outputDir,
                      showSession = showSession,
                      pdf = pdf,
                      runDimRed = runDimRed,
                      plotJackStraw = plotJackStraw,
                      plotElbowPlot = plotElbowPlot,
                      plotHeatmaps = plotHeatmaps,
                      forceRun = forceRun
                    ),
                    output_file = outputFile,
                    output_dir = outputDir,
                    intermediates_dir = outputDir,
                    knit_root_dir = outputDir)
  
  path <- paste0(outputDir, "SCE_SeuratDimRed", "-", gsub(" ", "_", Sys.Date()), ".rds")
  saveRDS(data, path)
  message("Output SCE object stored as ", paste0("SCE_SeuratDimRed", "-", gsub(" ", "_", Sys.Date()), ".rds"), " in ", outputDir, ".")
  message("Output HTML file stored as ", outputFile, " in ", outputDir, ".")
  
  return(data)
}


#' Generates an HTML report for Seurat Normalization 
#'  and returns the SCE object with the results computed and stored
#'  inside the object.
#' @param inSCE Input \code{\link[SingleCellExperiment]{SingleCellExperiment}}
#'  object previously passed through \code{reportSeuratRun()}. 
#' @param outputFile Specify the name of the generated output HTML file. 
#'  If \code{NULL} then the output file name will be based on the name of the 
#'  Rmarkdown template. Default \code{NULL}.
#' @param outputDir Specify the name of the output directory to save the 
#'  rendered HTML file. If \code{NULL} the file is stored to the current 
#'  working directory. Default \code{NULL}.
#' @param subtitle A character value specifying the subtitle to use in the 
#'  report. Default \code{NULL}.
#' @param authors A character value specifying the names of the authors to use 
#'  in the report. Default \code{NULL}.
#' @param showSession 	A logical value indicating if session information 
#'  should be displayed or not. Default is \code{FALSE}.
#' @param pdf A logical value indicating if a pdf should also be generated for 
#'  each figure in the report. Default is \code{FALSE}.
#' @param forceRun A logical value indicating if all computations previously
#'  computed should be re-calculated regardless if these computations are
#'  available in the input object. Default is \code{FALSE}.
#'
#' @return A \code{\link[SingleCellExperiment]{SingleCellExperiment}} object
#'  with computations stored.
#' @export
reportSeuratNormalization <- function(inSCE,
                               outputFile = NULL,
                               outputDir = NULL,
                               subtitle = NULL,
                               authors =  NULL,
                               showSession = FALSE,
                               pdf = FALSE,
                               forceRun = FALSE){
  
  if(is.null(outputDir)){
    outputDir <- getwd()
    message("No output directory defined, using current working directory ", outputDir, " instead.")
  }
  
  data <- inSCE
  
  rmarkdown::render(system.file("rmarkdown/seurat/reportSeuratNormalizeData.Rmd",
                                package="singleCellTK"),
                    params = list(
                      subtitle = subtitle,
                      authors = authors,
                      sce = data,
                      outputPath = outputDir,
                      showSession = showSession,
                      pdf = pdf,
                      forceRun = forceRun
                    ),
                    output_file = outputFile,
                    output_dir = outputDir,
                    intermediates_dir = getwd(),
                    knit_root_dir = getwd())
  
  path <- paste0(outputDir, "SCE_SeuratNormalization", "-", gsub(" ", "_", Sys.Date()), ".rds")
  saveRDS(data, path)
  message("Output SCE object stored as ", paste0("SCE_SeuratNormalization", "-", gsub(" ", "_", Sys.Date()), ".rds"), " in ", outputDir, ".")
  message("Output HTML file stored as ", outputFile, " in ", outputDir, ".")
  
  return(data)
}


#' Generates an HTML report for Seurat Feature Selection and returns the
#' SCE object with the results computed and stored inside the object.
#' @param inSCE Input \code{\link[SingleCellExperiment]{SingleCellExperiment}}
#'  object.
#' @param variable.features A numeric value indicating the number of top variable
#'  features to identify. Default \code{2000}.
#' @param runHVG A logical value indicating if the feature selection algorithm
#'  should be run or not. Default \code{TRUE}.
#' @param plotHVG A logical value indicating if the mean-to-variance plot
#'  of the top variable feature should be visualized or not. Default \code{TRUE}.
#' @param outputFile Specify the name of the generated output HTML file. 
#'  If \code{NULL} then the output file name will be based on the name of the 
#'  Rmarkdown template. Default \code{NULL}.
#' @param outputDir Specify the name of the output directory to save the 
#'  rendered HTML file. If \code{NULL} the file is stored to the current 
#'  working directory. Default \code{NULL}.
#' @param subtitle A character value specifying the subtitle to use in the 
#'  report. Default \code{NULL}.
#' @param authors A character value specifying the names of the authors to use 
#'  in the report. Default \code{NULL}.
#' @param showSession 	A logical value indicating if session information 
#'  should be displayed or not. Default is \code{FALSE}.
#' @param pdf A logical value indicating if a pdf should also be generated for 
#'  each figure in the report. Default is \code{FALSE}.
#' @param forceRun A logical value indicating if all computations previously
#'  computed should be re-calculated regardless if these computations are
#'  available in the input object. Default is \code{FALSE}.
#'
#' @return A \code{\link[SingleCellExperiment]{SingleCellExperiment}} object
#'  with computations stored.
#' @export
reportSeuratFeatureSelection <- function(inSCE,
                                         variable.features = 2000,
                                         runHVG = TRUE,
                                         plotHVG = TRUE,
                                         outputFile = NULL,
                                         outputDir = NULL,
                                         subtitle = NULL,
                                         authors =  NULL,
                                         showSession = FALSE,
                                         pdf = FALSE,
                                         forceRun = FALSE){
  if(is.null(outputDir)){
    outputDir <- getwd()
    message("No output directory defined, using current working directory ", outputDir, " instead.")
  }
  
  data <- inSCE
  
  rmarkdown::render(system.file("rmarkdown/seurat/reportSeuratFeatureSelection.Rmd",
                                package="singleCellTK"),
                    params = list(
                      subtitle = subtitle,
                      authors = authors,
                      sce = data,
                      variable.features = variable.features,
                      outputPath = outputDir,
                      showSession = showSession,
                      pdf = pdf,
                      forceRun = forceRun
                    ),
                    output_file = outputFile,
                    output_dir = outputDir,
                    intermediates_dir = getwd(),
                    knit_root_dir = getwd())
  
  path <- paste0(outputDir, "SCE_SeuratFeatureSelection", "-", gsub(" ", "_", Sys.Date()), ".rds")
  saveRDS(data, path)
  message("Output SCE object stored as ", paste0("SCE_SeuratFeatureSelection", "-", gsub(" ", "_", Sys.Date()), ".rds"), " in ", outputDir, ".")
  message("Output HTML file stored as ", outputFile, " in ", outputDir, ".")
  
  return(data)
}


#' Generates an HTML report for Seurat Scaling 
#'  and returns the SCE object with the results computed and stored
#'  inside the object.
#' @param inSCE Input \code{\link[SingleCellExperiment]{SingleCellExperiment}}
#'  object.
#' @param outputFile Specify the name of the generated output HTML file. 
#'  If \code{NULL} then the output file name will be based on the name of the 
#'  Rmarkdown template. Default \code{NULL}.
#' @param outputDir Specify the name of the output directory to save the 
#'  rendered HTML file. If \code{NULL} the file is stored to the current 
#'  working directory. Default \code{NULL}.
#' @param subtitle A character value specifying the subtitle to use in the 
#'  report. Default \code{NULL}.
#' @param authors A character value specifying the names of the authors to use 
#'  in the report. Default \code{NULL}.
#' @param showSession 	A logical value indicating if session information 
#'  should be displayed or not. Default is \code{FALSE}.
#' @param pdf A logical value indicating if a pdf should also be generated for 
#'  each figure in the report. Default is \code{FALSE}.
#' @param forceRun A logical value indicating if all computations previously
#'  computed should be re-calculated regardless if these computations are
#'  available in the input object. Default is \code{FALSE}.
#'
#' @return A \code{\link[SingleCellExperiment]{SingleCellExperiment}} object
#'  with computations stored.
#' @export
reportSeuratScaling <- function(inSCE,
                           outputFile = NULL,
                           outputDir = NULL,
                           subtitle = NULL,
                           authors =  NULL,
                           showSession = FALSE,
                           pdf = FALSE,
                           forceRun = FALSE){
  
  if(is.null(outputDir)){
    outputDir <- getwd()
    message("No output directory defined, using current working directory ", outputDir, " instead.")
  }
  
  data <- inSCE
  
  rmarkdown::render(system.file("rmarkdown/seurat/reportSeuratScaleData.Rmd",
                                package="singleCellTK"),
                    params = list(
                      subtitle = subtitle,
                      authors = authors,
                      sce = data,
                      outputPath = outputDir,
                      showSession = showSession,
                      pdf = pdf,
                      forceRun = forceRun
                    ),
                    output_file = outputFile,
                    output_dir = outputDir,
                    intermediates_dir = getwd(),
                    knit_root_dir = getwd())
  
  path <- paste0(outputDir, "SCE_SeuratScaleData", "-", gsub(" ", "_", Sys.Date()), ".rds")
  saveRDS(data, path)
  message("Output SCE object stored as ", paste0("SCE_SeuratScaleData", "-", gsub(" ", "_", Sys.Date()), ".rds"), " in ", outputDir, ".")
  message("Output HTML file stored as ", outputFile, " in ", outputDir, ".")
  
  return(data)
}


#' Generates an HTML report for Seurat Clustering and returns the SCE object 
#'  with the results computed and stored inside the object.
#' @param inSCE Input \code{\link[SingleCellExperiment]{SingleCellExperiment}}
#'  object.
#' @param biological.group A character value that specifies the name of the 
#'  \code{colData()} column to use as the main biological group in the Seurat 
#'  report for marker selection and grouping.
#' @param phenotype.groups A character vector that specifies the names of the 
#'  \code{colData()} columns to use for differential expression in addition to 
#'  the \code{biological.group} parameter.
#' @param runClustering A logical value indicating if Clustering should be run
#'  or not in the report. Default is \code{TRUE}. If \code{FALSE}, parameters
#'   \code{plotTSNE} and \code{plotUMAP} are also set to \code{FALSE}.
#' @param plotTSNE A logical value indicating if TSNE plots should be visualized
#'  in the clustering section of the report. Default is \code{TRUE}.
#' @param plotUMAP A logical value indicating if UMAP plots should be visualized
#'  in the clustering section of the report. Default is \code{TRUE}.
#' @param minResolution A numeric value indicating the minimum resolution to use
#'  for clustering. Default \code{0.3}.
#' @param maxResolution A numeric value indicating the maximum resolution to use
#'  for clustering. Default \code{1.5}.
#' @param numClusters temp (to remove)
#' @param significant_PC temp (change to pc.use)
#' @param outputFile Specify the name of the generated output HTML file. 
#'  If \code{NULL} then the output file name will be based on the name of the 
#'  Rmarkdown template. Default \code{NULL}.
#' @param outputDir Specify the name of the output directory to save the 
#'  rendered HTML file. If \code{NULL} the file is stored to the current 
#'  working directory. Default \code{NULL}.
#' @param subtitle A character value specifying the subtitle to use in the 
#'  report. Default \code{NULL}.
#' @param authors A character value specifying the names of the authors to use 
#'  in the report. Default \code{NULL}.
#' @param showSession 	A logical value indicating if session information 
#'  should be displayed or not. Default is \code{FALSE}.
#' @param pdf A logical value indicating if a pdf should also be generated for 
#'  each figure in the report. Default is \code{FALSE}.
#' @param forceRun A logical value indicating if all computations previously
#'  computed should be re-calculated regardless if these computations are
#'  available in the input object. Default is \code{FALSE}.
#'
#' @return A \code{\link[SingleCellExperiment]{SingleCellExperiment}} object
#'  with computations stored.
#' @export
reportSeuratClustering <- function(inSCE,
                                   biological.group = NULL,
                                   phenotype.groups = NULL,
                                   runClustering = TRUE,
                                   plotTSNE = TRUE,
                                   plotUMAP = TRUE,
                                   minResolution = 0.3,
                                   maxResolution = 1.5,
                                   numClusters = 10,
                                   significant_PC = 10,
                                   outputFile = NULL,
                                   outputDir = NULL,
                                   subtitle = NULL,
                                   authors =  NULL,
                                   showSession = FALSE,
                                   pdf = FALSE,
                                   forceRun = FALSE){
  
  if(is.null(biological.group)){
    stop("Must specify atleast one biological.group that is present in the colData of input object.")
  }
  
  if(!biological.group %in% names(colData(inSCE))){
    stop(biological.group, " not found in the colData of input object.")
  }
  
  if(!is.null(phenotype.groups)){
    if(!all(phenotype.groups %in% names(colData(inSCE)))){
      stop(phenotype.groups, " not found in the colData of input object.")
    }
  }
  
  data <- inSCE
  
  rmarkdown::render(system.file("rmarkdown/seurat/reportSeuratClustering.Rmd",
                                package="singleCellTK"),
                    params = list(
                      subtitle = subtitle,
                      authors = authors,
                      sce = data,
                      biological.group = biological.group,
                      phenotype.groups = phenotype.groups,
                      outputPath = outputDir,
                      showSession = showSession,
                      pdf = pdf,
                      runClustering = runClustering,
                      plotTSNE = plotTSNE,
                      plotUMAP = plotUMAP,
                      minResolution = minResolution,
                      maxResolution = maxResolution,
                      numClusters = numClusters,
                      significant_PC = significant_PC,
                      forceRun = forceRun
                    ),
                    output_file = outputFile,
                    output_dir = outputDir,
                    intermediates_dir = getwd(),
                    knit_root_dir = getwd())
  
  path <- paste0(outputDir, "SCE_SeuratClustering", "-", gsub(" ", "_", Sys.Date()), ".rds")
  saveRDS(data, path)
  message("Output SCE object stored as ", paste0("SCE_SeuratClustering", "-", gsub(" ", "_", Sys.Date()), ".rds"), " in ", outputDir, ".")
  message("Output HTML file stored as ", outputFile, " in ", outputDir, ".")
  
  return(data)
}


#' Generates an HTML report for Seurat Results (including Clustering & Marker
#'  Selection) and returns the SCE object with the results computed and stored
#'  inside the object.
#' @param inSCE Input \code{\link[SingleCellExperiment]{SingleCellExperiment}}
#'  object. 
#' @param biological.group A character value that specifies the name of the 
#'  \code{colData()} column to use as the main biological group in the Seurat 
#'  report for marker selection and grouping.
#' @param phenotype.groups A character vector that specifies the names of the 
#'  \code{colData()} columns to use for differential expression in addition to 
#'  the \code{biological.group} parameter.
#' @param selected.markers A character vector containing the user-specified 
#'  gene symbols or feature names of marker genes that be used to generate 
#'  gene plots in addition to the gene markers computed from 
#'  differential expression.
#' @param runMarkerSelection A logical value indicating if the marker selection
#'  computation should be run or not. Default \code{TRUE}.
#' @param plotMarkerSelection A logical value indicating if the gene marker
#'  plots should be visualized or not. Default \code{TRUE}.
#' @param countFeatures A numeric value indicating the number of top features
#'  to visualize in each group. Default \code{10}.
#' @param outputFile Specify the name of the generated output HTML file. 
#'  If \code{NULL} then the output file name will be based on the name of the 
#'  Rmarkdown template. Default \code{NULL}.
#' @param outputDir Specify the name of the output directory to save the 
#'  rendered HTML file. If \code{NULL} the file is stored to the current 
#'  working directory. Default \code{NULL}.
#' @param subtitle A character value specifying the subtitle to use in the 
#'  report. Default \code{NULL}.
#' @param authors A character value specifying the names of the authors to use 
#'  in the report. Default \code{NULL}.
#' @param showSession 	A logical value indicating if session information 
#'  should be displayed or not. Default is \code{FALSE}.
#' @param pdf A logical value indicating if a pdf should also be generated for 
#'  each figure in the report. Default is \code{FALSE}.
#'
#' @return A \code{\link[SingleCellExperiment]{SingleCellExperiment}} object
#'  with computations stored.
#' @export
reportSeuratMarkerSelection <- function(inSCE,
                                        biological.group = NULL,
                                        phenotype.groups = NULL,
                                        selected.markers = NULL,
                                        runMarkerSelection = TRUE,
                                        plotMarkerSelection = TRUE,
                                        countFeatures = 10,
                                        outputFile = NULL,
                                        outputDir = NULL,
                                        subtitle = NULL,
                                        authors =  NULL,
                                        showSession = FALSE,
                                        pdf = FALSE){
  
  if(is.null(biological.group)){
    stop("Must specify atleast one biological.group that is present in the colData of input object.")
  }
  
  if(!biological.group %in% names(colData(inSCE))){
    stop(biological.group, " not found in the colData of input object.")
  }
  
  if(!is.null(phenotype.groups)){
    if(!all(phenotype.groups %in% names(colData(inSCE)))){
      stop(phenotype.groups, " not found in the colData of input object.")
    }
  }
  
  if(is.null(outputDir)){
    outputDir <- getwd()
    message("No output directory defined, using current working directory ", outputDir, " instead.")
  }
  
  data <- inSCE
  
  rmarkdown::render(system.file("rmarkdown/seurat/reportSeuratMarkerSelection.Rmd",
                                package="singleCellTK"),
                    params = list(
                      subtitle = subtitle,
                      authors = authors,
                      sce = data,
                      biological.group = biological.group,
                      phenotype.groups = phenotype.groups,
                      selected.markers = selected.markers,
                      outputPath = outputDir,
                      showSession = showSession,
                      pdf = pdf,
                      runMarkerSelection = runMarkerSelection,
                      plotMarkerSelection = plotMarkerSelection,
                      countFeatures = countFeatures
                    ),
                    output_file = outputFile,
                    output_dir = outputDir,
                    intermediates_dir = getwd(),
                    knit_root_dir = getwd())
  
  path <- paste0(outputDir, "SCE_SeuratResults", "-", gsub(" ", "_", Sys.Date()), ".rds")
  saveRDS(data, path)
  message("Output SCE object stored as ", paste0("SCE_SeuratResults", "-", gsub(" ", "_", Sys.Date()), ".rds"), " in ", outputDir, ".")
  message("Output HTML file stored as ", outputFile, " in ", outputDir, ".")
  
  return(data)
}

#' Generates an HTML report for the complete Seurat workflow and returns the 
#'  SCE object with the results computed and stored inside the object.
#' @param inSCE Input \code{\link[SingleCellExperiment]{SingleCellExperiment}}
#'  object.
#' @param biological.group A character value that specifies the name of the 
#'  \code{colData()} column to use as the main biological group in the Seurat 
#'  report for marker selection and grouping.
#' @param phenotype.groups A character vector that specifies the names of the 
#'  \code{colData()} columns to use for differential expression in addition to 
#'  the \code{biological.group} parameter.
#' @param selected.markers A character vector containing the user-specified 
#'  gene symbols or feature names of marker genes that be used to generate 
#'  gene plots in addition to the gene markers computed from 
#'  differential expression.
#' @param clustering.resolution A numeric value indicating the user-specified 
#'  final resolution to use with clustering. Default is \code{0.8}.
#' @param variable.features A numeric value indicating the number of top 
#'  variable features to identify. Default \code{2000}.
#' @param pc.count A numeric value indicating the number of principal components
#'  to use in the analysis workflow. Default is \code{50}.
#' @param outputFile Specify the name of the generated output HTML file. 
#'  If \code{NULL} then the output file name will be based on the name of the 
#'  Rmarkdown template. Default \code{NULL}.
#' @param outputDir Specify the name of the output directory to save the 
#'  rendered HTML file. If \code{NULL} the file is stored to the current 
#'  working directory. Default \code{NULL}.
#' @param subtitle A character value specifying the subtitle to use in the 
#'  report. Default \code{NULL}.
#' @param authors A character value specifying the names of the authors to use 
#'  in the report. Default \code{NULL}.
#' @param showSession 	A logical value indicating if session information 
#'  should be displayed or not. Default is \code{FALSE}.
#' @param pdf A logical value indicating if a pdf should also be generated for 
#'  each figure in the report. Default is \code{FALSE}.
#' @param runHVG A logical value indicating if the feature selection
#'  computation should be run or not. Default is \code{TRUE}.
#' @param plotHVG A logical value indicating if the plot for the top most
#'  variable genes should be visualized in a mean-to-variance plot.
#'  Default is \code{TRUE}.
#' @param runDimRed A logical value indicating if PCA should be computed.
#'  Default is \code{TRUE}.
#' @param plotJackStraw A logical value indicating if JackStraw plot be
#'  visualized for the principal components. Default is \code{FALSE}.
#' @param plotElbowPlot A logical value indicating if the ElbowPlot be
#'  visualized for the principal components. Default is \code{TRUE}.
#' @param plotHeatmaps A logical value indicating if heatmaps should be plotted
#'  for the principal components. Default is \code{TRUE}.
#' @param runClustering A logical value indicating if clustering section should
#'  be run in the report. Default is \code{TRUE}.
#' @param plotTSNE A logical value indicating if TSNE plots should be visualized
#'  for clustering results. Default is \code{TRUE}.
#' @param plotUMAP A logical value indicating if the UMAP plots should be
#'  visualized for the clustering results. Default is \code{TRUE}.
#' @param minResolution A numeric value indicating the minimum resolution to
#'  use for clustering. Default is \code{0.3}.
#' @param maxResolution A numeric value indicating the maximum resolution to use
#'  for clustering. Default is \code{1.5}.
#' @param runMSClusters A logical value indicating if marker selection should
#'  be run between clusters. Default is \code{TRUE}.
#' @param runMSBioGroup A logical value indicating if marker selection should
#'  be run between the \code{biological.group} parameter. 
#'  Default is \code{TRUE}.
#' @param forceRun A logical value indicating if all algorithms should be
#'  re-run regardless if they have been computed previously in the input object.
#'  Default is \code{FALSE}. 
#'
#' @return A \code{\link[SingleCellExperiment]{SingleCellExperiment}} object
#'  with computations stored.
#' @export
reportSeurat <- function(
  inSCE,
  biological.group = NULL,
  phenotype.groups = NULL,
  selected.markers = NULL,
  clustering.resolution = 0.8,
  variable.features = 2000,
  pc.count = 50,
  outputFile = NULL,
  outputDir = NULL,
  subtitle = NULL,
  authors =  NULL,
  showSession = FALSE,
  pdf = FALSE,
  runHVG = TRUE,
  plotHVG = TRUE,
  runDimRed = TRUE,
  plotJackStraw = FALSE,
  plotElbowPlot = TRUE,
  plotHeatmaps = TRUE,
  runClustering = TRUE,
  plotTSNE = TRUE,
  plotUMAP = TRUE,
  minResolution = 0.3,
  maxResolution = 1.5,
  runMSClusters = TRUE,
  runMSBioGroup = TRUE,
  forceRun = FALSE){
  
  if(is.null(biological.group)){
    stop("Must specify atleast one biological.group that is present in the colData of input object.")
  }
  
  if(!biological.group %in% names(colData(inSCE))){
    stop(biological.group, " not found in the colData of input object.")
  }
  
  if(!is.null(phenotype.groups)){
    if(!all(phenotype.groups %in% names(colData(inSCE)))){
      stop(phenotype.groups, " not found in the colData of input object.")
    }
  }
  
  if(is.null(outputDir)){
    outputDir <- getwd()
    message("No output directory defined, using current working directory ", outputDir, " instead.")
  }
  
  data <- inSCE
  
  rmarkdown::render(system.file("rmarkdown/seurat/reportSeurat.Rmd",
                                package="singleCellTK"),
                    params = list(
                      subtitle = subtitle,
                      authors = authors,
                      sce = data,
                      biological.group = biological.group,
                      phenotype.groups = phenotype.groups,
                      selected.markers = selected.markers,
                      clustering.resolution = clustering.resolution,
                      variable.features = variable.features,
                      pc.count = pc.count,
                      outputPath = outputDir,
                      showSession = showSession,
                      pdf = pdf,
                      runHVG = runHVG,
                      plotHVG = plotHVG,
                      runDimRed = runDimRed,
                      plotJackStraw = plotJackStraw,
                      plotElbowPlot = plotElbowPlot,
                      plotHeatmaps = plotHeatmaps,
                      runClustering = runClustering,
                      plotTSNE = plotTSNE,
                      plotUMAP = plotUMAP,
                      minResolution = minResolution,
                      maxResolution = maxResolution,
                      runMSClusters = runMSClusters,
                      runMSBioGroup = runMSBioGroup,
                      forceRun = forceRun
                    ),
                    output_file = outputFile,
                    output_dir = outputDir,
                    intermediates_dir = getwd(),
                    knit_root_dir = getwd())
  
  path <- paste0(outputDir, "SCE_SeuratReport", "-", gsub(" ", "_", Sys.Date()), ".rds")
  saveRDS(data, path)
  message("Output SCE object stored as ", paste0("SCE_SeuratReport", "-", gsub(" ", "_", Sys.Date()), ".rds"), " in ", outputDir, ".")
  message("Output HTML file stored as ", outputFile, " in ", outputDir, ".")
  
  return(data)
}

#' @title Get diffAbundanceFET .html report
#' @description A function to generate .html Rmarkdown report containing the visualizations of the diffAbundanceFET function output
#' @param inSCE A \code{\link[SingleCellExperiment]{SingleCellExperiment}}
#' object.
#' @param cluster A single \code{character}, specifying the name to store the
#' cluster label in \code{\link{colData}}.
#' @param variable A single \code{character}, specifying the name to store the
#' phenotype labels in \code{\link{colData}}.
#' @param control \code{character}. Specifying one or more categories that can
#' be found in the vector specified by \code{variable}.
#' @param case \code{character}. Specifying one or more categories that can
#' be found in the vector specified by \code{variable}.
#' @param analysisName A single \code{character}. Will be used for naming the
#' result table, which will be saved in metadata slot.
#' @param output_dir name of the output directory to save the rendered file. If
#' \code{NULL} the file is stored to the current working directory.
#' Default \code{NULL}.
#' @param output_file name of the generated file. If \code{NULL} then the output
#' file name will be based on the name of the Rmarkdown template. Default
#' \code{NULL}.
#' @param pdf A \code{logical} value indicating if a pdf should also be
#'  generated for each figure in the report. Default is \code{TRUE}.
#' @param showSession A \code{logical} value indicating if session information
#'  should be displayed or not. Default is \code{TRUE}.
#' @return An HTML file of the report will be generated at the path specified
#' in the arguments.
#' @export
reportDiffAbundanceFET <-
    function(inSCE,
             cluster,
             variable,
             control,
             case,
             analysisName,
             output_dir = ".",
             output_file = "DifferentialAbundanceFET_Report",
             pdf = FALSE,
             showSession = TRUE) {
        inSCE <- diffAbundanceFET(inSCE, cluster, variable, control, 
                                  case, analysisName)
        rmarkdown::render(
            system.file("rmarkdown/DifferentialAbundanceFET_Report.Rmd",
                        package="singleCellTK"),
            params = list(
               sce = inSCE,
                analysisName = analysisName,
                pdf = isTRUE(pdf),
                showSession = isTRUE(showSession)
            ),
           output_file = output_file,
           output_dir = output_dir
        )
    }

#' @title Get plotClusterAbundance .html report
#' @description A function to generate .html Rmarkdown report containing the 
#' visualizations of the plotClusterAbundance function output
#' @param inSCE A \code{\link[SingleCellExperiment]{SingleCellExperiment}}
#' object.
#' @param cluster A single \code{character}, specifying the name to store the
#' cluster label in \code{\link{colData}}.
#' @param variable A single \code{character}, specifying the name to store the
#' phenotype labels in \code{\link{colData}}.
#' @param output_dir name of the output directory to save the rendered file. If
#' \code{NULL} the file is stored to the current working directory.
#' Default \code{NULL}.
#' @param output_file name of the generated file. If \code{NULL} then the output
#' file name will be based on the name of the Rmarkdown template. Default
#' \code{NULL}.
#' @param pdf A \code{logical} value indicating if a pdf should also be
#'  generated for each figure in the report. Default is \code{TRUE}.
#' @param showSession A \code{logical} value indicating if session information
#'  should be displayed or not. Default is \code{TRUE}.
#' @return An HTML file of the report will be generated at the path specified
#' in the arguments.
#' @export
reportClusterAbundance <- function(inSCE,
                                   cluster,
                                   variable,
                                   output_dir = ".",
                                   output_file = "plotClusterAbundance_Report",
                                   pdf = FALSE,
                                   showSession = TRUE) {
    rmarkdown::render(
        system.file("rmarkdown/PlotClusterAbundance_Report.Rmd",
        package="singleCellTK"),
        params = list(
            sce = inSCE,
            cluster = cluster,
            variable = variable,
            pdf = isTRUE(pdf),
            showSession = isTRUE(showSession)
        ),
        output_file = output_file,
        output_dir = output_dir
    )
}
