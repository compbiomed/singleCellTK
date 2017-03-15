library(singleCellTK)
library(shiny)
library(shinyjs)
library(plotly)
library(d3heatmap)

source("helpers.R")

clusterChoice <- ''
sampleChoice <- ''
alertText <- ''
pcComponents <- ''
if(!is.null(getShinyOption("inputSCEset"))){
  clusterChoice <- colnames(pData(getShinyOption("inputSCEset")))
  sampleChoice <- rownames(pData(getShinyOption("inputSCEset")))
  pcComponents <- paste("PC",1:nrow(pData(getShinyOption("inputSCEset"))),sep="")
  alertText <- HTML("<div class='alert alert-success alert-dismissible'>\
                    <span class='glyphicon glyphicon-ok' aria-hidden='true'>\
                    </span> Successfully Uploaded from Command Line! <button \
                    type='button' class='close' data-dismiss='alert'>&times;\
                    </button></div>")
}

# Define UI for application that draws a histogram
shinyUI(
  navbarPage(
    paste("Single Cell Toolkit v",packageVersion("singleCellTK"), sep=""),
    #bootstrap theme
    theme = "bootstrap.min.css",
    #Upload Tab
    tabPanel(
      "Upload",
      useShinyjs(),
      tags$style(appCSS),
      tags$div(
        class="jumbotron",
        tags$div(
          class="container",
          h1("Single Cell Toolkit"),
          p("Filter, cluster, and analyze single cell RNA-Seq data")
        )
      ),
      tags$div(
        class="container",
        tags$div(id="uploadAlert", alertText),
        fileInput('countsfile', 'Upload a matrix of counts here',
                  accept = c(
                    'text/csv',
                    'text/comma-separated-values',
                    'text/tab-separated-values',
                    'text/plain',
                    '.csv',
                    '.tsv'
                  )
        ),
        fileInput('annotfile', 'Optional: Upload a matrix of annotations here',
                  accept = c(
                    'text/csv',
                    'text/comma-separated-values',
                    'text/tab-separated-values',
                    'text/plain',
                    '.csv',
                    '.tsv'
                  )
        ),
        fileInput('featurefile', 'Optional: Upload a matrix of feature annotations here',
                  accept = c(
                    'text/csv',
                    'text/comma-separated-values',
                    'text/tab-separated-values',
                    'text/plain',
                    '.csv',
                    '.tsv'
                  )
        ),
        withBusyIndicatorUI(
          actionButton("uploadData", "Upload")
        ),
        tags$div(
          class="container",
          p("")
        )
      ),
      includeHTML('www/footer.html')
    ),#tab "upload" end
    tabPanel(
      "Data Summary and Filtering",
      tags$div(
        class="container",
        h1("Data Summary"),
        fluidPage(
          fluidRow(
            column(8, tableOutput('summarycontents')),
            column(
              4,
              wellPanel(
                checkboxInput("removeNoexpress", "Remove genes with 0 expression across all samples (Recommended)", value=TRUE),
                numericInput('minDetectGenect', label = 'Minimum Detected Genes per Sample.', value=1700, min = 1, max = 100000),
                numericInput("LowExpression", "% Low Gene Expression to Filter",value=40, min = 0, max = 100),
                h2("Delete Outliers"),
                selectInput("deletesamplelist","Select Samples:",
                            sampleChoice,
                            multiple = TRUE),
                actionButton("filterData", "Filter Data"),
                actionButton("resetData", "Reset")
              )
            )
          ),
          fluidRow(
            dataTableOutput('contents')
          )
        )
      ),
      includeHTML('www/footer.html')
    ),#tab "Data Summary and Filtering" end
    
    tabPanel(
      "DR & Clustering",
      tabsetPanel(
        tabPanel(
          "PCA",
          fluidRow(
            column(
              4,
              wellPanel(
                selectInput(inputId ="pcaAlgorithm",label = "Algorithm",choices = c("regular PCA","randomized PCA","robust PCA","RR PCA")),
                checkboxGroupInput(inputId = "pcaCheckbox",label = "Involving variables", choices = c("PC1","PC2","PC3"),selected = c("PC1","PC2","PC3")),
                selectInput(inputId = "selectAdditionalVariables",label = "Additional variables:",choices = clusterChoice,multiple = TRUE),
                selectInput("plotTypeId","Plot Type",c("Paired Plot","Single Plot"),"Paired Plot"),
                selectInput("colorClusters","Color Clusters By",clusterChoice),
                withBusyIndicatorUI(actionButton("plotPCA", "Plot PCA Data"))
              )#wellPanel
            ),
            column(
              8,
              plotOutput("pcaPlot")
            )
          ),#fluidRow
          mainPanel(
            h1("Instructions"),
            p(""), strong("Principle Component Analysis:"), 
            p("1. Choose algorithm (one of the PCA algorithm)"),
            p("2. Choose components to be involved in paired plot"),
            p("3. Choose 2 components in Checkbox to produce single plot(if more than 2 components are selected, only the first 2 will be used)"),
            p("3. Choose feature to color data by"),
            p("4. Visualize your data")
          )
        ),#subtab "Principle Component Analysis" end
        
        tabPanel("Dimensionality Reduction(tSNE)", 
                 fluidRow(
                   column(4,
                          wellPanel(
                          selectInput("selectDimRed","Algorithm",c("PCA","tSNE")),
                          selectInput("pcX", "X axis:", pcComponents),
                          selectInput("pcY", "Y axis:", pcComponents, selected = "PC2"),
                          selectInput("colorClusters","Color Clusters By",clusterChoice),
                          withBusyIndicatorUI(actionButton("plotData", "Plot Data"))
                          )),
                   column(8,
                          plotlyOutput("dimredPlot"))
  
                 ),
                 mainPanel(
                   h1("Instructions"),
                   p(""), strong("tSNE:"), 
                   p("1. Choose algorithm (PCA or tSNE)"),
                   p("2. Choose which components to use for the x and y axes"),
                   p("3. Choose feature to color data by"),
                   p("4. Visualize your data")
                 )
        ),#tSNE end
        tabPanel("Clustering",
                 fluidRow(
                   column(4,
                          wellPanel("Select clustering algorithm:",
                            selectInput("selectCluster", "Algorithm", c("K-Means")),
                            selectInput("selectK", "Cluster Centers", c(1, 2, 3, 4, 5)),
                            p("Select data to cluster:"),
                            selectInput("selectDataC", "Data", c("Raw Data", "PCA Components", "tSNE Components")),
                            selectInput("pcX", "X axis:", c("1"="PC1","2"="PC2","3"="PC3")),
                            selectInput("pcY", "Y axis:", c("1"="PC1","2"="PC2","3"="PC3")),
                            selectInput("colorClusters","Color Clusters By",c("Tissue")),
                            actionButton("plotClusters", "Plot Clusters")
                          )
                   ),
                   column(8,
                          plotlyOutput("clusterPlot"))
                 ),
                 mainPanel(
                   h1("Instructions"),
                   p(""), strong("Clustering:"),
                   p("1. Choose clustering algorithm (K-Means, ...)"),
                   p("2. Choose which data to use (raw data, principal component values, tSNE values)"), 
                   p("3. Choose which components to use for the x and y axes"),
                   p("4. Choose feature to color data by"),
                   p("5. Visualize your data and clusters")
                 )
      )
    ),#Clustering end
    includeHTML("www/footer.html")
    ),#tab "DR & Clustering" end
    tabPanel(
      "Differential Expression",
      tags$div(
        class="container",
        h1("Differential Expression"),
        fluidPage(
          fluidRow(
            column(4,
                   wellPanel(
                     selectInput("selectDiffex","Differential Expression",c("DESeq", "DESeq2", "limma")),
                     selectInput("selectDiffex_condition","Select Condition",clusterChoice),
                     sliderInput("selectNGenes", "Display Top N Genes:", 5, 500, 500, 5),
                     checkboxInput("applyCutoff", "Apply p-value Cutoff"),
                     checkboxInput("clusterRows", "Cluster Heatmap Rows", value=TRUE),
                     checkboxInput("clusterColumns", "Cluster Heatmap Columns", value=TRUE),
                     sliderInput("selectPval", "p-value cutoff:", 0.01, 0.2, 0.05),
                     selectInput("selectCorrection","Correction Type",c("FDR")),
                     withBusyIndicatorUI(actionButton("runDiffex", "Run Differential Expression")),
                     downloadButton("downloadGeneList","Download Results")
                   )),
            column(8,
                   tabsetPanel(
                     id = 'dataset',
                     tabPanel('Heatmap', 
                              fluidPage(
                                #fluidRow(
                                  column(4,
                                      wellPanel(
                                        checkboxInput("displayHeatmapRowLabels", "Display Row Labels", value=TRUE),
                                        checkboxInput("displayHeatmapColumnLabels", "Display Column Labels", value=TRUE),
                                        checkboxInput("displayHeatmapColumnDendrograms", "Display Column Dendrograms", value=TRUE),
                                        checkboxInput("displayHeatmapRowDendrograms", "Display Row Dendrograms", value=TRUE)
                                      )),
                                  column(8,
                                         plotOutput("diffPlot") )
                                #)
                              )),
                     tabPanel('Results Table', dataTableOutput('diffextable')),
                     tabPanel('Interactive Heatmap', d3heatmapOutput("interactivediffPlot"))
                   )
            )
          )
        )
      ),
      includeHTML('www/footer.html')
    ),#tab "Diff Expre" end
    tabPanel(
      "Subsampling",
      tags$div(
        class="container",
        h1("Subsampling"),
        fluidPage(
          fluidRow(
            column(4,
                   wellPanel(
                     selectInput("subCovariate", "Covariate for differential expression", clusterChoice),
                     selectInput("selectDiffMethod","Differential Expression Method",c("tpm.t","DESeq")),
                     numericInput('minSim', label = 'Minimum subsample Size.', value=1000, min = 1, max = 1000000),
                     numericInput('maxSim', label = 'Maximum subsample Size.', value=10000, min = 100, max = 100000000),
                     numericInput('iterations', label = 'Number of bootstrap iterations.', value=10, min = 2, max = 1000),
                     actionButton("runSubsample", "Run subsampler"),
                     actionButton("runDifferentialPower", "Run differential power analysis")
                   )),
            column(8,
                   plotOutput("downDone"))
          ),
          fluidRow(
            plotOutput("powerBoxPlot")
          )
        )
      ),
      includeHTML('www/footer.html')
    ),
    navbarMenu(
      "More",
      tabPanel(
        "DE Heatmap",
        class="container",
        h1("DE Heatmap"),
        fluidPage(
          fluidRow(
            column(4,
                   wellPanel(
                     selectInput("selectHeatmap","Select Heatmap",c("Standard","Complex","Interactive")),
                     actionButton("makeHeatmap", "Generate Heatmap")
                   )),
            column(8,
                   wellPanel(
                     d3heatmapOutput("heatmapPlot"))
            )
          )
        )
      ),
      tabPanel(
        "Batch Correction",
        tags$div(
          class="container",
          h1("Batch Correction")
        ),
        includeHTML('www/footer.html')
      ),      
      tabPanel(
        "Pathway Activity Analysis",
        tags$div(
          class="container",
          h1("Pathway Activity Analysis")
        ),
        includeHTML('www/footer.html')
      )
    )
  )
)
