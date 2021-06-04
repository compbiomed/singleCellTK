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
                              cl = FALSE,
                              curated = FALSE)
{
  ns <- session$ns
  
  output$ui <- renderUI({
    bsCollapsePanel("Downstream Analysis", 
                    uiOutput(ns("de")),
                    uiOutput(ns("cv")),
                    uiOutput(ns("pa")),
                    uiOutput(ns("qcf")),
                    uiOutput(ns("nbc")),
                    uiOutput(ns("cw")),
                    uiOutput(ns("dr")),
                    uiOutput(ns("fs")),
                    uiOutput(ns("cl")),
                    uiOutput(ns("curated")),
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
        h5("Discover quantitative changes between experimental conditions using one of the many integrated statistical frameworks:"),
        actionButton(inputId = ns("goCV"), label = "Go to Cell Viewer!")
      )
    })
  }
  
  if(pa){
    output$pa <- renderUI({
      panel(
        heading = "Pathway Analysis",
        h5("Explore biological activity or functions with pathway analysis using either 'GSVA' or 'EnrichR' statistical frameworks:"),
        actionButton(inputId = ns("goGSVA"), label = "Go to GSVA!"),
        actionButton(inputId = ns("goEnrichR"), label = "Go to EnrichR!"),
      )
    })
  }
  
  if(qcf){
    output$qcf <- renderUI({
      panel(
        heading = "Quality Control & Filtering",
        h5("Perform quality control checks / filtering."),
        actionButton(inputId = ns("goQC"), label = "Go to Quality Control!")
      )
    })
  }
  
  if(nbc){
    output$nbc <- renderUI({
      panel(
        heading = "Normalization & Batch-Correction",
        h5("go nbc"),
        actionButton(inputId = ns("goNBC"), label = "Go to Normalization/Batch-Correction!")
      )
    })
  }
  
  if(cw){
    output$cw <- renderUI({
      panel(
        heading = "Curated Workflows",
        h5("go cw"),
        actionButton(inputId = ns("goCelda"), label = "Go to Celda!"),
        actionButton(inputId = ns("goSeurat"), label = "Go to Seurat!")
      )
    })
  }
  
  if(dr){
    output$dr <- renderUI({
      panel(
        heading = "Dimensionality Reduction",
        h5("go dr"),
        actionButton(inputId = ns("goDR"), label = "Go to Dimensionality Reduction!")
      )
    })
  }
  
  if(fs){
    output$fs <- renderUI({
      panel(
        heading = "Feature Selection",
        h5("go fs"),
        actionButton(inputId = ns("goFS"), label = "Go to Feature Selection!")
      )
    })
  }
  
  if(cl){
    output$cl <- renderUI({
      panel(
        heading = "Clustering",
        h5("go cl"),
        actionButton(inputId = ns("goCL"), label = "Go to Clustering!")
      )
    })
  }
  
  if(curated){
    output$curated <- renderUI({
      panel(
        heading = "Curated Workflows",
        h5("go curated"),
        actionButton(inputId = ns("goSeurat"), label = "Go to Seurat!"),
        actionButton(inputId = ns("goSeurat"), label = "Go to Celda!")
      )
    })
  }
  
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
    
    observeEvent(input$goGSVA,{
      showTab(inputId = "navbar",
              target = "GSVA",
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
              target = "DR",
              select = TRUE,
              session = parent)
    })
    
    observeEvent(input$goFS,{
      showTab(inputId = "navbar",
              target = "FS",
              select = TRUE,
              session = parent)
    })
    
    observeEvent(input$goCL,{
      showTab(inputId = "navbar",
              target = "CL",
              select = TRUE,
              session = parent)
    })
}