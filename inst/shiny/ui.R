library(shiny)
library(shinyjs)
library(ComplexHeatmap)
library(biomaRt)
library(circlize)
library(limma)
library(d3heatmap)
library(ggplot2)
library(plotly)
library(DESeq)
library(GGally)
library(data.table)
library(MAST)
library(rsvd)
library(pcaMethods)
library(colourpicker)
library(gridExtra)
library(cluster)
library(ggtree)
library(ape)
library(GSVA)
library(GSVAdata)
library(SingleCellExperiment)
library(singleCellTK)

source("helpers.R")
source("colourGroupInput.R")
data("c2BroadSets")

clusterChoice <- ""
sampleChoice <- ""
featureChoice <- ""
geneChoice <- ""
alertText <- ""
pcComponents <- ""
numClusters <- ""
currassays <- ""
numSamples <- 30
pcComponents_selectedY <- NULL
if (!is.null(getShinyOption("inputSCEset"))){
  numSamples <- ncol(getShinyOption("inputSCEset"))
  clusterChoice <- colnames(colData(getShinyOption("inputSCEset")))
  geneChoice <- rownames(getShinyOption("inputSCEset"))
  sampleChoice <- colnames(getShinyOption("inputSCEset"))
  featureChoice <- colnames(rowData(getShinyOption("inputSCEset")))
  pcComponents <- paste("PC", 1:numSamples, sep = "")
  pcComponents_selectedY <- pcComponents[2]
  numClusters <- 1:numSamples
  currassays <- names(assays(getShinyOption("inputSCEset")))
  alertText <- HTML("<div class='alert alert-success alert-dismissible'>\
                    <span class='glyphicon glyphicon-ok' aria-hidden='true'>\
                    </span> Successfully Uploaded from Command Line! <button \
                    type='button' class='close' data-dismiss='alert'>&times;\
                    </button></div>")
}

source("ui_01_upload.R", local = TRUE) #creates shiny_panel_upload variable
source("ui_02_filter.R", local = TRUE) #creates shiny_panel_filter variable
source("ui_03_cluster.R", local = TRUE) #creates shiny_panel_cluster variable
source("ui_04_diffex.R", local = TRUE) #creates shiny_panel_diffex variable
source("ui_05_subsample.R", local = TRUE) #creates shiny_panel_subsample variable
source("ui_06_batchcorrect.R", local = TRUE) #creates shiny_panel_batchcorrect variable
source("ui_07_pathway.R", local = TRUE) #creates shiny_panel_pathway variable
source("ui_08_mast.R", local = TRUE) #creates shiny_panel_mast variable

# Define UI for application that draws a histogram
shinyUI(
  navbarPage(
    paste("Single Cell Toolkit v", packageVersion("singleCellTK"), sep = ""),
    #bootstrap theme
    theme = "bootstrap.min.css",
    #Upload Tab
    tabPanel("Upload", shiny_panel_upload),
    tabPanel("Data Summary and Filtering", shiny_panel_filter),
    tabPanel("DR & Clustering", shiny_panel_cluster),
    tabPanel("Batch Correction", shiny_panel_batchcorrect),
    navbarMenu(
      "Differential Expression",
      tabPanel("Differential Expression", shiny_panel_diffex),
      tabPanel("MAST", shiny_panel_mast)
    ),
    tabPanel("Pathway Activity Analysis", shiny_panel_pathway),
    tabPanel("Sample Size", shiny_panel_subsample),
    footer = includeHTML("www/footer.html")
  )
)
