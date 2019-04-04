library(shiny)
library(shinyjs)
library(ComplexHeatmap)
library(limma)
library(ggplot2)
library(plotly)
library(data.table)
library(MAST)
library(colourpicker)
library(gridExtra)
library(cluster)
library(ggtree)
library(ape)
library(GSVA)
library(GSVAdata)
library(shinyalert)
library(enrichR)
library(matrixStats)
library(Biobase)
library(base)
library(SingleCellExperiment)
library(singleCellTK)
library(celda)
library(shinycssloaders)
library(shinythemes)
library(umap)

source("helpers.R")
source("colourGroupInput.R")
data("c2BroadSets")

#test internet connection for enrichR connectivity
internetConnection <- suppressWarnings(Biobase::testBioCConnection())

clusterChoice <- ""
sampleChoice <- ""
featureChoice <- ""
geneChoice <- ""
alertText <- ""
pcComponents <- ""
numClusters <- ""
currassays <- ""
currreddim <- ""
if (internetConnection){
  enrichedDB <- enrichR::listEnrichrDbs()$libraryName
} else {
  enrichedDB <- ""
}
numSamples <- 30
pcComponentsSelectedY <- NULL
if (!is.null(getShinyOption("inputSCEset"))){
  numSamples <- ncol(getShinyOption("inputSCEset"))
  clusterChoice <- colnames(colData(getShinyOption("inputSCEset")))
  geneChoice <- rownames(getShinyOption("inputSCEset"))
  sampleChoice <- colnames(getShinyOption("inputSCEset"))
  featureChoice <- colnames(rowData(getShinyOption("inputSCEset")))
  pcComponents <- paste("PC", 1:numSamples, sep = "")
  pcComponentsSelectedY <- pcComponents[2]
  numClusters <- 1:numSamples
  currassays <- names(assays(getShinyOption("inputSCEset")))
  currreddim <- names(reducedDims(getShinyOption("inputSCEset")))
  alertText <- HTML("<div class='alert alert-success alert-dismissible'>\
                    <span class='glyphicon glyphicon-ok' aria-hidden='true'>\
                    </span> Successfully Uploaded from Command Line! <button \
                    type='button' class='close' data-dismiss='alert'>&times;\
                    </button></div>")
}

if (is.null(getShinyOption("theme"))){
  shinyTheme <- "flatly"
} else {
  shinyTheme <- getShinyOption("theme")
}

source("ui_01_upload.R", local = TRUE) #creates shinyPanelUpload variable
source("ui_02_filter.R", local = TRUE) #creates shinyPanelFilter variable
source("ui_03_1_genewise_vis.R", local = TRUE) #creates shinyPanelCluster variable
source("ui_03_2_samplewise_vis.R", local = TRUE) #creates shinyPanelCluster variable
source("ui_03_3_celda.R", local = TRUE) #creates shinyPanelCelda variable
source("ui_04_batchcorrect.R", local = TRUE) #creates shinyPanelBatchcorrect variable
source("ui_05_1_diffex.R", local = TRUE) #creates shinyPanelDiffex variable
source("ui_05_2_mast.R", local = TRUE) #creates shinyPanelMAST variable
source("ui_06_1_pathway.R", local = TRUE) #creates shinyPanelPathway variable
source("ui_06_2_enrichR.R", local = TRUE) #creates shinyPanelEnrichR variable
source("ui_07_subsample.R", local = TRUE) #creates shinyPanelSubsample variable

if (is.null(getShinyOption("includeVersion"))){
  tooltitle <- paste("Single Cell Toolkit v",
                     packageVersion("singleCellTK"), sep = "")
} else {
  if (getShinyOption("includeVersion")){
    tooltitle <- paste("Single Cell Toolkit v",
                       packageVersion("singleCellTK"), sep = "")
  } else {
    tooltitle <- "Single Cell Toolkit"
  }
}

# Define UI for application that draws a histogram
shinyUI(
  navbarPage(
    tooltitle,
    theme = shinytheme(shinyTheme),
    #Upload Tab
    tabPanel("Upload", shinyPanelUpload),
    tabPanel("Data Summary & Filtering", shinyPanelFilter),
    navbarMenu(
      "Visualization & Clustering",
      tabPanel("Genewise Visualization", shinyPanelVis),
      tabPanel("Samplewise Vis & Clustering", shinyPanelCluster),
      tabPanel("Celda", shinyPanelCelda)
    ),
    tabPanel("Batch Correction", shinyPanelBatchcorrect),
    navbarMenu(
      "Differential Expression",
      tabPanel("Differential Expression", shinyPanelDiffex),
      tabPanel("MAST", shinyPanelMAST)
    ),
    navbarMenu(
      "Enrichment Analysis",
      tabPanel("GSVA", shinyPanelPathway),
      tabPanel("EnrichR", shinyPanelEnrichR)
    ),
    tabPanel("Sample Size", shinyPanelSubsample),
    footer = includeHTML("www/footer.html")
  )
)
