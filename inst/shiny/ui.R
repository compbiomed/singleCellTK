library(shiny)
library(shinyjs)
library(shinyFiles)
library(ComplexHeatmap)
library(limma)
library(ggplot2)
library(plotly)
library(data.table)
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
library(shinyWidgets);
library(shinyBS);
library(shinyjqui);
library(Seurat);
library(ggplotify);
library(ggplot2);
library(cowplot);
library(tidyverse)
library(dplyr)
library(readxl)
library(broom)
library(RColorBrewer)
library(grDevices)
library(shinyWidgets)
library(stringr)
library(Hmisc)


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
curraltExps <- ""
#from SCE
cell_list <- ""
gene_list <- ""
#from assays
method_list <- ""
#from reduced
approach_list <- ""
#from colData
annotation_list <- ""
#from RColorBrewer
colorbrewer_list <- rownames(RColorBrewer::brewer.pal.info)
color_table <- RColorBrewer::brewer.pal.info %>% data.frame()
color_seqdiv <- rownames(color_table[which(color_table$category == "div"
                                           |color_table$category == "seq"),])
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
  curraltExps <- names(altExp(getShinyOption("inputSCEset")))
  ###############################################################
  #from sce
  cell_list <- BiocGenerics::colnames(getShinyOption("inputSCEset"))
  gene_list <- BiocGenerics::rownames(getShinyOption("inputSCEset"))
  #from assays
  method_list <- names(assays(getShinyOption("inputSCEset")))
  #from reduced
  approach_list <- names(reducedDims(getShinyOption("inputSCEset")))
  #from colData
  annotation_list <- names(colData(getShinyOption("inputSCEset")))
  #from colorbrewer
  colorbrewer_list <- rownames(RColorBrewer::brewer.pal.info)
  color_table <- RColorBrewer::brewer.pal.info %>% data.frame()
  color_seqdiv <- rownames(color_table[which(color_table$category == "div"|color_table$category == "seq"),])
  ###############################################################
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

source("ui_01_import.R", local = TRUE) #creates shinyPanelImport variable
source("ui_01_gene_sets.R", local = TRUE) #creates shinyPanelGeneSets variable
source("ui_01_columnAnnotation.R", local = TRUE) #creates shinyPanelColumnAnnotation variable
source("ui_01_rowAnnotation.R", local = TRUE) #creates shinyPanelRowAnnotation variable
source("ui_export.R", local = TRUE) #creates shinyPanelExport variable
source("ui_02_qc_filter.R", local = TRUE) #creates shinyPanelQCFilter variable
source("ui_03_2_cluster.R", local = TRUE) #creates shinyPanelCluster variable
source("ui_09_3_celdaWorkflow.R", local = TRUE) #creates shinyPanelCelda variable
source("ui_04_batchcorrect.R", local = TRUE) #creates shinyPanelBatchcorrect variable
source("ui_04_fs_dimred.R", local = TRUE) #creates shinyPanelFS_DimRed variable
source("ui_05_1_diffex.R", local = TRUE) #creates shinyPanelDiffex variable
source("ui_05_2_findMarker.R", local = TRUE) #creates shinyPanelfindMarker variable
source("ui_06_1_pathway.R", local = TRUE) #creates shinyPanelPathway variable
source("ui_06_2_enrichR.R", local = TRUE) #creates shinyPanelEnrichR variable
source("ui_07_subsample.R", local = TRUE) #creates shinyPanelSubsample variable
source("ui_08_viewers.R", local = TRUE) #creates shinyPanelViewers variable
source("ui_08_2_cellviewer_v2.R", local = TRUE) #creates shinyPanelCellViewer variable
source("ui_08_3_heatmap.R", local = TRUE) #creates shinyPanelHeatmap variable
source("ui_09_curatedworkflows.R", local = TRUE) #creates shinyPanelCuratedWorkflows variable
source("ui_09_2_seuratWorkflow.R", local = TRUE) #creates shinyPanelSeurat variable
source("ui_export.R", local = TRUE) #creates shinyPanelExport variable

jsCode <- "

shinyjs.disableTabs = function() {
  let tabs = $('.nav li a').not('a[data-value=\"Data\"], a[data-value=\"Import\"]');
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
      id = "navbar",
      # selected="CellViewer",
      # theme = shinytheme(shinyTheme),
      theme = shinytheme("yeti"),
      navbarMenu(
        "Data",
        tabPanel("Import Single Cell Data", shinyPanelImport),
        tabPanel("Import Gene Sets", shinyPanelGeneSets),
        tabPanel("Column Annotation", shinyPanelColumnAnnotation),
        tabPanel("Row Annotation", shinyPanelRowAnnotation),
        tabPanel("Export Single Cell Data", shinyPanelExport)
      ),
      tabPanel("QC & Filtering", shinyPanelQCFilter),
      tabPanel("Normalization & Batch Correction", shinyPanelBatchcorrect),
      tabPanel("Feature Selection & Dimensionality Reduction", shinyPanelFS_DimRed),
      tabPanel("Clustering", shinyPanelCluster),
      navbarMenu(
        "Differential Expression & Marker Selection",
        tabPanel("Differential Expression", shinyPanelDiffex),
        tabPanel("Find Marker", shinyPanelfindMarker)
      ),
      navbarMenu(
        "Cell Annotation & Pathway Analysis",
        tabPanel("GSVA", shinyPanelPathway),
        tabPanel("EnrichR", shinyPanelEnrichR)
      ),
      tabPanel("Sample Size Calculator", shinyPanelSubsample),
      navbarMenu(
        "Curated Workflows",
        tabPanel("Celda", shinyPanelCelda),
        tabPanel("Seurat", shinyPanelSeurat),
        tabPanel("Bioconductor/OSCA", h1("Bioconductor/OSCA"))
      ),
      # tabPanel("Curated Workflows", shinyPanelCuratedWorkflows),
      navbarMenu("Viewers",
                 tabPanel("Gene Visualization", shinyPanelViewers),
                 tabPanel("Cell Viewer", value="CellViewer", shinyPanelCellViewer),
                 tabPanel("Heatmap", shinyPanelHeatmap)),
      footer = includeHTML("www/footer.html"),
      fluidRow(
        column(12, id = "consoleDiv",
               actionButton(inputId="consoleToggle", label = "Console Log"),
               hidden(verbatimTextOutput(outputId="console")),
        )
      ),
      useShinyjs(),
      extendShinyjs(text = jsCode, functions = c("enableTabs", "disableTabs"))
    )
)

