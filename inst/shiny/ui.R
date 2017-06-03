library(singleCellTK)
library(shiny)
library(shinyjs)
library(plotly)
library(d3heatmap)
library(ape)

source("helpers.R")
source("colourGroupInput.R")

clusterChoice <- ''
sampleChoice <- ''
featureChoice <- ''
geneChoice <- ''
alertText <- ''
pcComponents <- ''
numClusters <- ''
if(!is.null(getShinyOption("inputSCEset"))){
  clusterChoice <- colnames(pData(getShinyOption("inputSCEset")))
  geneChoice <- rownames(exprs(getShinyOption("inputSCEset"))[1:100])
  sampleChoice <- rownames(pData(getShinyOption("inputSCEset")))
  featureChoice <- colnames(fData(getShinyOption("inputSCEset")))
  pcComponents <- paste("PC",1:nrow(pData(getShinyOption("inputSCEset"))),sep="")
  numClusters <- 1:nrow(pData(getShinyOption("inputSCEset")))
  alertText <- HTML("<div class='alert alert-success alert-dismissible'>\
                    <span class='glyphicon glyphicon-ok' aria-hidden='true'>\
                    </span> Successfully Uploaded from Command Line! <button \
                    type='button' class='close' data-dismiss='alert'>&times;\
                    </button></div>")
}

source("ui_01_upload.R", local=TRUE) #creates shiny_panel_upload variable
source("ui_02_filter.R", local=TRUE) #creates shiny_panel_filter variable
source("ui_03_cluster.R", local=TRUE) #creates shiny_panel_cluster variable
source("ui_04_diffex.R", local=TRUE) #creates shiny_panel_diffex variable
source("ui_05_subsample.R", local=TRUE) #creates shiny_panel_subsample variable
source("ui_06_batchcorrect.R", local=TRUE) #creates shiny_panel_batchcorrect variable
source("ui_07_pathway.R", local=TRUE) #creates shiny_panel_pathway variable
source("ui_08_mast.R", local=TRUE) #creates shiny_panel_mast variable

# Define UI for application that draws a histogram
shinyUI(
  navbarPage(
    paste("Single Cell Toolkit v",packageVersion("singleCellTK"), sep=""),
    #bootstrap theme
    theme = "bootstrap.min.css",
    #Upload Tab
    tabPanel("Upload", shiny_panel_upload),
    tabPanel("Data Summary and Filtering", shiny_panel_filter),
    tabPanel("DR & Clustering", shiny_panel_cluster),
    tabPanel("Differential Expression", shiny_panel_diffex),
    tabPanel("Subsampling", shiny_panel_subsample),
    tabPanel("MAST", shiny_panel_mast),
    navbarMenu(
      "More",
      tabPanel("Batch Correction", shiny_panel_batchcorrect),      
      tabPanel("Pathway Activity Analysis", shiny_panel_pathway)
    )
  )
)
