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
source("ui_celda.R", local = TRUE) #creates shinyPanelCelda variable
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
                 tabPanel("Cell Viewer", value="CellViewer", fluidPage(tags$div(
                   class = "container",
                   h1("Cell Viewer"),
                   radioGroupButtons(
                     "viewertabs",
                     choices = c("Scatter Plot", "Bar Plot", "Violin/Box Plot"),
                     selected = NULL
                   ),
                   fluidRow(column(
                     3,
                     wellPanel(
                       # Section 1 - Assay Settings
                       actionButton("cv_button1", h4(strong("Select Coordinates"))),
                       # open by default
                       tags$div(
                         id = "cv_collapse1",
                         selectInput(
                           inputId = "QuickAccess",
                           label = NULL,
                           choices = c("", approach_list, "Custom")
                         ),
                         #-+-+-+-+-+-X-Axis###################################
                         conditionalPanel(
                           condition = sprintf("input['%s'] == 'Custom'", "QuickAccess"),
                           h5(strong("X-Axis")),
                           selectInput(
                             "TypeSelect_Xaxis",
                             h5("Data"),
                             choices = c("Reduced Dimensions", "Expression Assays", "Cell Annotation")
                           ),
                           #Reduced Dimensions condition
                           conditionalPanel(
                             condition = sprintf("input['%s'] == 'Reduced Dimensions'", "TypeSelect_Xaxis"),
                             selectizeInput(
                               "ApproachSelect_Xaxis",
                               label = h5("Approach"),
                               choices = c(approach_list)
                             ),
                             selectInput("ColumnSelect_Xaxis", h5("Dimension"), choices = NULL)
                           ),
                           
                           #Expression Assays condition
                           conditionalPanel(
                             condition = sprintf("input['%s'] == 'Expression Assays'", "TypeSelect_Xaxis"),
                             selectizeInput(
                               "AdvancedMethodSelect_Xaxis",
                               label = h5("Advanced Method"),
                               choices = c(method_list)
                             ),
                             selectizeInput(
                               "GeneSelect_Assays_Xaxis",
                               label = h5("Feature"),
                               choices = c(gene_list)
                             )
                           ),
                           
                           #Cell Annotation condition
                           conditionalPanel(
                             condition = sprintf("input['%s'] == 'Cell Annotation'", "TypeSelect_Xaxis"),
                             selectizeInput(
                               "AnnotationSelect_Xaxis",
                               label = h5("Annotation"),
                               choices = c(annotation_list)
                             )
                           ),
                           
                           #-+-+-+-+-+-Y-Axis###################################
                           h5(strong("Y-Axis")),
                           selectInput(
                             "TypeSelect_Yaxis",
                             h5("Data"),
                             choices = c("Reduced Dimensions", "Expression Assays", "Cell Annotation")
                           ),
                           #Reduced Dimensions condition
                           conditionalPanel(
                             condition = sprintf("input['%s'] == 'Reduced Dimensions'", "TypeSelect_Yaxis"),
                             selectizeInput(
                               "ApproachSelect_Yaxis",
                               label = h5("Approach"),
                               choices = c(approach_list)
                             ),
                             selectInput("ColumnSelect_Yaxis", h5("Dimension"), choices = NULL)
                           ),
                           
                           #Expression Assays condition
                           conditionalPanel(
                             condition = sprintf("input['%s'] == 'Expression Assays'", "TypeSelect_Yaxis"),
                             selectizeInput(
                               "AdvancedMethodSelect_Yaxis",
                               label = h5("Advanced Method"),
                               choices = c(method_list)
                             ),
                             selectizeInput(
                               "GeneSelect_Assays_Yaxis",
                               label = h5("Feature"),
                               choices = c(gene_list)
                             )
                           ),
                           
                           #Cell Annotation condition
                           conditionalPanel(
                             condition = sprintf("input['%s'] == 'Cell Annotation'", "TypeSelect_Yaxis"),
                             selectizeInput(
                               "AnnotationSelect_Yaxis",
                               label = h5("Annotation"),
                               choices = c(annotation_list)
                             )
                           )
                           
                         )
                       ),
                       
                       #-+-+-+-+-+-colorby part1###################################
                       tags$hr(),
                       #Select Color by Data
                       # Section 1 - Assay Settings
                       actionButton("cv_button2", h4(strong("Color"))),
                       
                       # open by default
                       tags$div(
                         id = "cv_collapse2",
                         selectInput(
                           'TypeSelect_Colorby', h5(strong('Color By')), choices = c(
                             "Single Color",
                             "Reduced Dimensions",
                             "Expression Assays",
                             "Cell Annotation"
                           ),
                         ),
                         # Single Color condition
                         conditionalPanel(
                           condition = sprintf("input['%s'] == 'Single Color'", "TypeSelect_Colorby"),
                           colourInput("Col", "", "purple", palette = 'limited')
                         ),
                         #Reduced Dimensions condition
                         conditionalPanel(
                           condition = sprintf("input['%s'] == 'Reduced Dimensions'", "TypeSelect_Colorby"),
                           selectizeInput(
                             "ApproachSelect_Colorby",
                             label = h5("Approach"),
                             choices = c(approach_list)
                           ),
                           selectInput("ColumnSelect_Colorby", h5("Dimension"), choices = NULL)
                         ),
                         
                         #Expression Assays condition
                         conditionalPanel(
                           condition = sprintf("input['%s'] == 'Expression Assays'", "TypeSelect_Colorby"),
                           selectizeInput(
                             "AdvancedMethodSelect_Colorby",
                             label = h5("Advanced Method"),
                             choices = c(method_list)
                           ),
                           selectizeInput(
                             "GeneSelect_Assays_Colorby",
                             label = h5("Feature"),
                             choices = c(gene_list)
                           )
                         ),
                         
                         #Cell Annotation condition
                         conditionalPanel(
                           condition = sprintf("input['%s'] == 'Cell Annotation'", "TypeSelect_Colorby"),
                           selectizeInput(
                             "AnnotationSelect_Colorby",
                             label = h5("Annotation"),
                             choices = c(annotation_list)
                           )
                         ),
                         
                         #-+-+-+-+-+-colorby part2###################################
                         conditionalPanel(
                           # condition = sprintf("input['%s'] != 'Single Color' && output.hide_typebtns == 'show'", "TypeSelect_Colorby"),
                           condition = sprintf("input['%s'] != 'Single Color'", "TypeSelect_Colorby"),
                           radioButtons(
                             "SelectColorType",
                             label = NULL,
                             choices = c("Categorical", "Continuous")
                           ),
                           conditionalPanel(
                             id = "continuousColorConditional",
                             condition = sprintf("input['%s'] == 'Continuous'", "SelectColorType"),
                             colourInput("highColor", "High Color", "blue", "background", "limited"),
                             colourInput("midColor", "Middle Color", "#666666", "background", "limited"),
                             colourInput("lowColor", "Low Color", "white", "background", "limited")
                           ), 
                           conditionalPanel(
                             id = "categoricalColorConditional",
                             condition = sprintf("input['%s'] == 'Categorical'", "SelectColorType"),
                             uiOutput("categoricalColorUI")
                           ),
                           conditionalPanel(
                             id="binningConditional",
                             condition = sprintf("input['%s'] == 'Continuous'", "SelectColorType"),
                             tags$hr(),
                             h5(style="display: inline-block; margin-top: 0px; margin-bottom: 20px","Perform Binning"),
                             switchInput(
                               inputId = "checkColorbinning",
                               onLabel = "Yes",
                               offLabel = "No",
                               value=FALSE,
                               size="mini",
                               inline = TRUE
                             )
                           ),
                           
                           conditionalPanel(
                             condition =  "input.checkColorbinning == 1",
                             numericInput(
                               "adjustColorbinning",
                               h5("Number of Bins"),
                               value = 2,
                               min = 2
                             )
                           )
                           #,
                           
                           
                           #selectizeInput("adjustbrewer", h5(strong("Color Palettes")), choices = NULL)
                         )
                       ),
                       #-+-+-+-+-+-group by###################################
                       tags$hr(),
                       shinyjs::useShinyjs(),
                       # Section 1 - Assay Settings
                       actionButton("cv_button3", h4(strong("Group By"))),
                       # open by default
                       tags$div(
                         id = "cv_collapse3",
                         selectizeInput(
                           inputId = "adjustgroupby",
                           label = NULL,
                           choices = c("None", annotation_list)
                         )
                         #,
                         #       conditionalPanel(condition = sprintf("input['%s'] != 'None'", "adjustgroupby"),
                         #                       radioButtons("SelectValueType",label = NULL,choices = c("Categorical", "Continuous")),
                         #                       conditionalPanel(condition = sprintf("input['%s'] == 'Continuous'", "SelectValueType"),
                         #                                         checkboxInput("checkbinning",h5("Perform Binning:"), value = FALSE)),
                         #                       conditionalPanel(condition = "input.checkbinning == 1",
                         #                                         numericInput("adjustbinning", h5("Number of Bins:"),value = 2, min =2))
                         #       )
                       ),
                       tags$hr(),
                       conditionalPanel(
                         id="violinConditional",
                         condition = sprintf("input['%s'] == 'Violin/Box Plot'", "viewertabs"),
                         # checkboxInput("vlnboxcheck", "Violin plot", value = FALSE),
                         h5(style="display: inline-block; margin-top: 0px; margin-bottom: 20px","Use Violin Plot"),
                         switchInput(
                           inputId = "vlnboxcheck",
                           onLabel = "Yes",
                           offLabel = "No",
                           value=FALSE,
                           size="mini",
                           inline = TRUE
                         )
                       ),
                       actionButton("runCellViewer", "Plot")
                     )
                   ), #sidebarPanel_end
                   #-+-+-+-+-+-mainPanel#################################
                   column(
                     9,
                     wellPanel(
                       plotlyOutput("scatter", height = "600px") %>% withSpinner(size = 3, color = "#0dc5c1", type = 8),
                       tags$br(),
                       # conditionalPanel("$('#scatter').hasClass('recalculating')",
                       #                  tags$div('Your plot is loading, due to large manipulation.
                       #                           This message will disappear once the plot is generated.')),
                       tags$hr(),
                       fluidRow(
                         column(6, textInput("adjusttitle", h5(strong(
                           "Title:"
                         )))),
                         column(6, textInput("adjustlegendtitle", h5(
                           strong("Legend title:")
                         ))),
                         column(6, sliderInput(
                           "adjustlegendtitlesize",
                           h5(strong("Legend title size:")),
                           min = 1,
                           max = 20,
                           value = 12
                         )),
                         column(6, sliderInput(
                           "adjustlegendsize",
                           h5(strong("Legend size:")),
                           min = 1,
                           max = 20,
                           value = 10
                         )),
                         column(6, sliderInput(
                           "adjustalpha",
                           h5(strong("Opacity:")),
                           min = 0,
                           max = 1,
                           value = 1
                         )),
                         column(6, sliderInput(
                           "adjustsize",
                           h5(strong("Dot size:")),
                           min = 0.1,
                           max = 0.8,
                           value = 0.45
                         )),
                         column(6, textInput("adjustxlab", h5(
                           strong("X-axis label:")
                         ))),
                         column(6, textInput("adjustylab", h5(
                           strong("Y-axis label:")
                         ))),
                         column(6, sliderInput(
                           "adjustaxissize",
                           h5(strong("Axis size:")),
                           min = 1,
                           max = 20,
                           value = 10
                         )),
                         column(6, sliderInput(
                           "adjustaxislabelsize",
                           h5(strong("Axis label size:")),
                           min = 1,
                           max = 20,
                           value = 10
                         ))
                       )
                     )
                   ))
                 ))),
                 tabPanel("Heatmap", shinyPanelHeatmap)),
      footer = includeHTML("www/footer.html"),
      fluidRow(
        column(12, id = "consoleDiv",
               actionButton(inputId="consoleToggle", label = "Console Log"),
               hidden(verbatimTextOutput(outputId="console")),
        )
      ),
      useShinyjs(),
      extendShinyjs(text = jsCode)
    )
)
