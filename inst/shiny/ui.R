# Check if CRAN packages are installed, otherwise prompt user to install them.
requiredPackages <- c("shinyjqui", "shinyWidgets", "shinythemes", "shinyFiles",
                      "shinyBS", "shinybusy", "tidyverse")
if(!all(requiredPackages %in% installed.packages())){
  missingPackages <- requiredPackages[which(requiredPackages %in% installed.packages() == FALSE)]
  message("Installing missing packages: ")
  message(paste0(missingPackages, collapse = " "))
  install.packages(missingPackages)
}

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
library(tibble)
# library(pushbar)
# library(spsComps)


source("helpers.R")
source("colourGroupInput.R")
data("c2BroadSets")

#source modules
source("module_nonLinearWorkflow.R")
source("module_filterTable.R")
source("module_renameCluster.R")

docs.base <- paste0("https://www.camplab.net/sctk/v",
                    package.version("singleCellTK"), "/")
docs.artPath <- paste0(docs.base, "articles/")

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
currGS <- ""
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
  enrichedDB <- listEnrichrDbs()$libraryName
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
  currGS <- sctkListGeneSetCollections(getShinyOption("inputSCEset"))
  ###############################################################
  #from sce
  cell_list <- colnames(getShinyOption("inputSCEset"))
  gene_list <- rownames(getShinyOption("inputSCEset"))
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

source("ui_00_helloWorld.R")
source("ui_01_import.R", local = TRUE) #creates shinyPanelImport variable
source("ui_01_gene_sets.R", local = TRUE) #creates shinyPanelGeneSets variable
source("ui_01_columnAnnotation.R", local = TRUE) #creates shinyPanelColumnAnnotation variable
source("ui_01_rowAnnotation.R", local = TRUE) #creates shinyPanelRowAnnotation variable
source("ui_01_removeData.R", local = TRUE) #creates shinyPanelRemove variable
source("ui_export.R", local = TRUE) #creates shinyPanelExport variable
source("ui_02_qc_filter.R", local = TRUE) #creates shinyPanelQCFilter variable
source("ui_03_2_cluster.R", local = TRUE) #creates shinyPanelCluster variable
source("ui_09_3_celdaWorkflow.R", local = TRUE) #creates shinyPanelCelda variable
source("ui_04_batchcorrect.R", local = TRUE) #creates shinyPanelBatchcorrect variable
source("ui_04_fs_dimred.R", local = TRUE) #creates shinyPanelFS_DimRed variable
source("ui_05_1_diffex.R", local = TRUE) #creates shinyPanelDiffex variable
source("ui_05_2_findMarker.R", local = TRUE) #creates shinyPanelfindMarker variable
source("ui_05_3_cellTypeLabel.R", local = TRUE) # creates shinyPanelLabelCellType variable
#source("ui_06_1_pathway.R", local = TRUE) #creates shinyPanelPathway variable
source("ui_06_2_enrichR.R", local = TRUE) #creates shinyPanelEnrichR variable
source("ui_06_1_pathwayAnalysis.R", local = TRUE) #creates shinyPanelvam variable
source("ui_10_1_TSCAN.R", local = TRUE) #creates shinyPanelTSCAN variable
source("ui_07_subsample.R", local = TRUE) #creates shinyPanelSubsample variable
source("ui_08_2_cellviewer.R", local = TRUE) #creates shinyPanelCellViewer variable
source("ui_08_3_heatmap.R", local = TRUE) #creates shinyPanelHeatmap variable
source("ui_08_4_bubbleplot.R", local = TRUE) #creates shinyPanelBubbleplot variable
#source("ui_09_curatedworkflows.R", local = TRUE) #creates shinyPanelCuratedWorkflows variable
source("ui_09_2_seuratWorkflow.R", local = TRUE) #creates shinyPanelSeurat variable
source("ui_09_4_scanpyWorkflow.R", local = TRUE) #creates shinyPanelSeurat variable
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

jsScriptAutoScrollConsole <- "
function mutate(mutations) {
  mutations.forEach(function(mutation) {
    alert(mutation.type);
  });
}

function startAutoScroll() {
                var $panel = $('#consolePanel');
    $panel.animate({scrollTop: $panel.prop('scrollHeight')});
}

var target = document.querySelector('#consoleText')
var observer = new MutationObserver( mutate );
var config = { characterData: false, attributes: false, childList: true, subtree: false };
observer.observe(target, config);
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
        tabPanel("Cell Annotation", shinyPanelColumnAnnotation),
        tabPanel("Feature Annotation", shinyPanelRowAnnotation),
        tabPanel("Export Single Cell Data", shinyPanelExport),
        tabPanel("Delete Single Cell Data", shinyPanelRemove)
      ),
      tabPanel("QC & Filtering", shinyPanelQCFilter),
      tabPanel("Normalization & Batch Correction", shinyPanelBatchcorrect),
      tabPanel("Feature Selection & Dimensionality Reduction", shinyPanelFS_DimRed),
      tabPanel("Clustering", shinyPanelCluster),
      navbarMenu(
        "Differential Expression & Cell Type Labeling",
        tabPanel("Find Marker", shinyPanelfindMarker),
        tabPanel("Differential Expression", shinyPanelDiffex),
        tabPanel("Cell Type Labeling", shinyPanelLabelCellType)

      ),
      navbarMenu(
        "Enrichment & Pathway Analysis",
        tabPanel("EnrichR", shinyPanelEnrichR),
        tabPanel("Pathway Activity", shinyPanelvam)


      ),
      navbarMenu(
        "Trajectory Analysis",
        tabPanel("TSCAN", value = "TSCANWorkflow", shinyPanelTSCAN)
      ),

      tabPanel("Sample Size Calculator", shinyPanelSubsample),
      navbarMenu(
        "Curated Workflows",
        tabPanel("Celda", value = "CeldaWorkflow", shinyPanelCelda),
        tabPanel("Seurat", shinyPanelSeurat),
        tabPanel("Scanpy", shinyPanelScanpy)
      ),
      # tabPanel("Curated Workflows", shinyPanelCuratedWorkflows),
      navbarMenu("Viewers",
                 tabPanel("Cell Viewer", value="CellViewer", shinyPanelCellViewer),
                 tabPanel("Heatmap", shinyPanelHeatmap),
                 tabPanel("Bubbleplot", shinyPanelBubbleplot)
                 ),
      footer = includeHTML("www/logo.html"),
      fluidRow(
        column(12, id = "consoleDiv",
               actionButton(inputId="consoleToggle", label = "Console Log"),
               tags$head(
                 tags$script(HTML(jsScriptAutoScrollConsole))
               ),
               hidden(div(id = "consolePanel", style = "overflow-y:scroll;
                          max-height: 220px; width: 100%; background-color: white;
                          position: relative; bottom: 0; align: centre; padding: 0px;",
                          verbatimTextOutput(outputId="consoleText", placeholder = TRUE)
               ))
        )
      ),
      useShinyjs(),
      extendShinyjs(text = jsCode, functions = c("enableTabs", "disableTabs")),

      # Following lines of code add a loading spinner when toolkit launches and
      # loads several ui elements/plots etc.
      includeCSS("busy-load-piccard21.css"),
      tags$script(src = "initialLoading.js"),
      tags$script(src = "busy-load-piccard21.js"),
      
      # Add ability to track usage with Google Analytics. Requires a 
      # link like:
      # https://www.googletagmanager.com/gtag/js?id=G-XXXXXXXXXX
      # with a code for the Google analytics project. It also requires a .js
      # file. See the following page on the wiki for more info:
      # https://github.com/compbiomed/singleCellTK/wiki/Google-Analytics
      tags$head(
        shiny::tags$script(
          src = "https://www.googletagmanager.com/gtag/js?id=G-NP0B0KLYE2",
          async = ""
        ),
        shiny::tags$script(
          src = "gtag.js"
        )
      )
    )
)
