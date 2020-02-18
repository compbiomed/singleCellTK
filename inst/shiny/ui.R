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
source("ui_03_2_samplewise_vis.R", local = TRUE) #creates shinyPanelCluster variable
source("ui_celda.R", local = TRUE) #creates shinyPanelCelda variable
source("ui_04_batchcorrect.R", local = TRUE) #creates shinyPanelBatchcorrect variable
source("ui_04_fs_dimred.R", local = TRUE) #creates shinyPanelFS_DimRed variable
source("ui_05_1_diffex.R", local = TRUE) #creates shinyPanelDiffex variable
source("ui_05_2_mast.R", local = TRUE) #creates shinyPanelMAST variable
source("ui_06_1_pathway.R", local = TRUE) #creates shinyPanelPathway variable
source("ui_06_2_enrichR.R", local = TRUE) #creates shinyPanelEnrichR variable
source("ui_07_subsample.R", local = TRUE) #creates shinyPanelSubsample variable
source("ui_08_viewers.R", local = TRUE) #creates shinyPanelViewers variable
source("ui_09_curatedworkflows.R", local = TRUE) #creates shinyPanelCuratedWorkflows variable



jsCode <- "

shinyjs.disableTabs = function() {
  let tabs = $('.nav li a').not('a[data-value=\"Upload\"]');
  tabs.bind('click', function(e) {
    e.preventDefault();
    return false;
  });
  
  tabs.addClass('disabled');
}

shinyjs.enableTabs = function() {
  let tabs = $('.nav li a');
  tabs.unbind('click');
  tabs.removeClass('disabled');
}
"

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

shinyUI(
    navbarPage(
      tooltitle,
      theme = shinytheme(shinyTheme),
      #Upload Tab
      tabPanel("Upload", shinyPanelUpload),
      navbarMenu("QC & Filtering", 
                 tabPanel("Filtering"), shinyPanelFilter),
      # tabPanel(title="QC & Filtering", shinyPanelFilter),
      tabPanel("Normalization & Batch Correction", shinyPanelBatchcorrect),
      tabPanel("Feature Selection & Dimensionality Reduction", shinyPanelFS_DimRed),
      tabPanel("Clustering", shinyPanelCluster),
      navbarMenu(
        "Differential Expression & Marker Selection",
        tabPanel("Differential Expression", shinyPanelDiffex),
        tabPanel("MAST", shinyPanelMAST)
      ),
      navbarMenu(
        "Cell Annotation & Pathway Analysis",
        tabPanel("GSVA", shinyPanelPathway),
        tabPanel("EnrichR", shinyPanelEnrichR)
      ),
      tabPanel("Sample Size Calculator", shinyPanelSubsample),
      navbarMenu(
        "Curated Workflows",
        tabPanel("CELDA", shinyPanelCelda),
        tabPanel("Seurat", h1("Seurat")),
        tabPanel("Bioconductor/OSCA", h1("Bioconductor/OSCA"))
      ),
      # tabPanel("Curated Workflows", shinyPanelCuratedWorkflows),
      navbarMenu("Viewers", 
                 tabPanel("Gene Visualization", shinyPanelViewers)),
      footer = includeHTML("www/footer.html"),
      # fluidRow(
      #   column(12, id = "consoleDiv",
      #          actionButton(inputId="consoleToggle", label = "Show/Hide Console Log"),
      #          verbatimTextOutput(outputId="console"),
      #          tags$head(tags$style("#console {height: 150px; margin-bottom: 0}")),
      #          tags$head(tags$style("#consoleDiv {position: fixed; bottom: 0; z-index: 3; padding: 0px"))
      #   )
      # )
      useShinyjs(),
      extendShinyjs(text = jsCode)
    )
)
