# UI
nonLinearWorkflowUI <- function(id)
{
  ns <- NS(id)
  fluidPage(
    uiOutput(ns("ui"))
  )
}

# Server
nonLinearWorkflow <- function(input, output, session, parent, 
                              de = FALSE,
                              cv = FALSE, 
                              pa = FALSE,
                              qcf = FALSE,
                              nbc = FALSE,
                              cw = FALSE,
                              dr = FALSE,
                              fs = FALSE,
                              cl = FALSE)
{
  ns <- session$ns
  
  output$ui <- renderUI({
    bsCollapsePanel(tagList(icon("chevron-circle-down"), "Downstream Analysis"), 
                    uiOutput(ns("de")),
                    uiOutput(ns("pa")),
                    uiOutput(ns("qcf")),
                    uiOutput(ns("nbc")),
                    uiOutput(ns("cw")),
                    uiOutput(ns("dr")),
                    uiOutput(ns("fs")),
                    uiOutput(ns("cl")),
                    uiOutput(ns("cv")),
                    style = "success")
  })
  
  if(de){
    output$de <- renderUI({
                      panel(
                        heading = "Differential Expression & Marker Selection",
                        h5("Discover quantitative changes between experimental conditions using one of the many integrated statistical frameworks for differential expression:"),
                        actionButton(inputId = ns("goDE"), label = "Go to Differential Expression!"),
                        hr(),
                        h5("Compute and visualize marker genes from cluster classes using 'MAST', 'Limma' or 'DESeq2' methods:"),
                        actionButton(inputId = ns("goMS"), label = "Go to Marker Selection!"),
                      )
    })
  }
  
  if(cv){
    output$cv <- renderUI({
      panel(
        heading = "Cell Viewer",
        h5("Visualize data using scatter plot, violin plot, bar or box plots."),
        actionButton(inputId = ns("goCV"), label = "Go to Cell Viewer!")
      )
    })
  }
  
  if(pa){
    output$pa <- renderUI({
      panel(
        heading = "Pathway Analysis",
        h5("Explore biological activity or functions with pathway analysis using either 'GSVA' or 'EnrichR' statistical frameworks:"),
        actionButton(inputId = ns("goPathwayAnalysis"), label = "Go to Pathway Analysis!"),
        actionButton(inputId = ns("goEnrichR"), label = "Go to EnrichR!"),
      )
    })
  }
  
  if(qcf){
    output$qcf <- renderUI({
      panel(
        heading = "Quality Control & Filtering",
        h5("Perform quality control checks to get a better sense of the quality of the input data using a number of available methods. Additionally, perform filtering on the data to remove samples that do not pass a quality criteria."),
        actionButton(inputId = ns("goQC"), label = "Go to QC"),
        actionButton(inputId = ns("goFilter"), label = "Go to Filtering")
      )
    })
  }
  
  if(nbc){
    output$nbc <- renderUI({
      panel(
        heading = "Normalization & Batch-Correction",
        h5("Normalize the raw/filtered input data to remove biasness of technical variation from the data or additionally adjust for batch-effect if the data is processed in multiple batches."),
        actionButton(inputId = ns("goNBC"), label = "Go to Normalization/Batch-Correction!")
      )
    })
  }
  
  if(cw){
    output$cw <- renderUI({
      panel(
        heading = "Curated Workflows",
        h5("Interactively analyze raw or filtered data in a step-by-step curated workflow using either Seurat or Celda pipeline."),
        actionButton(inputId = ns("goSeurat"), label = "Go to Seurat"),
        actionButton(inputId = ns("goCelda"), label = "Go to Celda")
      )
    })
  }
  
  if(dr){
    output$dr <- renderUI({
      panel(
        heading = "Dimensionality Reduction",
        h5("Reduce high-dimensional normalized data into low-dimension using PCA/ICA or tSNE/UMAP which may further be used with some clustering algorithms as the recommended input data type."),
        actionButton(inputId = ns("goDR"), label = "Go to Dimensionality Reduction")
      )
    })
  }
  
  if(fs){
    output$fs <- renderUI({
      panel(
        heading = "Feature Selection",
        h5("Find and select a subset of features/genes that show the most heterogeneity among the data samples and use this subset in the downstream analysis."),
        actionButton(inputId = ns("goFS"), label = "Go to Feature Selection")
      )
    })
  }
  
  if(cl){
    output$cl <- renderUI({
      panel(
        heading = "Clustering",
        h5("Identify cluster labels for each of the sample of the input data by grouping together similar samples through a number of available clustering algorithms."),
        actionButton(inputId = ns("goCL"), label = "Go to Clustering")
      )
    })
  }
  
  observeEvent(input$goCV,{
    showTab(inputId = "navbar",
            target = "CellViewer",
            select = TRUE,
            session = parent)
  })
  
    observeEvent(input$goDE,{
      showTab(inputId = "navbar",
              target = "Differential Expression",
              select = TRUE,
              session = parent)
    })

    observeEvent(input$goMS,{
      showTab(inputId = "navbar",
              target = "Find Marker",
              select = TRUE,
              session = parent)
    })
    
    observeEvent(input$goPathwayAnalysis,{
      showTab(inputId = "navbar",
              target = "Pathway Activity",
              select = TRUE,
              session = parent)
    })
    
    observeEvent(input$goEnrichR,{
      showTab(inputId = "navbar",
              target = "EnrichR",
              select = TRUE,
              session = parent)
    })
    
    observeEvent(input$goQC,{
      showTab(inputId = "navbar",
              target = "QC & Filtering",
              select = TRUE,
              session = parent)
      updateTabsetPanel(session = parent, 
                        inputId = "QCFilterTabsetPanel", 
                        selected = "QC")
    })
    
    observeEvent(input$goFilter,{
      showTab(inputId = "navbar",
              target = "QC & Filtering",
              select = TRUE,
              session = parent)
      updateTabsetPanel(session = parent, 
                        inputId = "QCFilterTabsetPanel", 
                        selected = "Filtering")
    })
    
    observeEvent(input$goNBC,{
      showTab(inputId = "navbar",
              target = "Normalization & Batch Correction",
              select = TRUE,
              session = parent)
    })
    
    observeEvent(input$goCelda,{
      showTab(inputId = "navbar",
              target = "Celda",
              select = TRUE,
              session = parent)
    })
    
    observeEvent(input$goSeurat,{
      showTab(inputId = "navbar",
              target = "Seurat",
              select = TRUE,
              session = parent)
    })
    
    observeEvent(input$goDR,{
      showTab(inputId = "navbar",
              target = "Feature Selection & Dimensionality Reduction",
              select = TRUE,
              session = parent)
      updateTabsetPanel(session = parent, 
                        inputId = "FSDimRedTabsetPanel", 
                        selected = "Dimensionality Reduction")
    })
    
    observeEvent(input$goFS,{
      showTab(inputId = "navbar",
              target = "Feature Selection & Dimensionality Reduction",
              select = TRUE,
              session = parent)
      updateTabsetPanel(session = parent, 
                        inputId = "FSDimRedTabsetPanel", 
                        selected = "Feature Selection")
    })
    
    observeEvent(input$goCL,{
      showTab(inputId = "navbar",
              target = "Clustering",
              select = TRUE,
              session = parent)
    })
}
