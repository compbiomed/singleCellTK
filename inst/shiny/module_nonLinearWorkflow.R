# #UI
# nonLinearWorkflowUI <- function(id, 
#                                 de = FALSE, 
#                                 cv = FALSE, 
#                                 pa = FALSE, 
#                                 qcf = FALSE, 
#                                 nbc = FALSE, 
#                                 cw = FALSE, 
#                                 dr = FALSE, 
#                                 fs = FALSE,
#                                 cl = FALSE)
# {
#   ns <- NS(id)
#   fluidPage(
#     hidden(textInput(inputId = ns("de"),
#                      label = "",
#                      value = de),
#            textInput(inputId = ns("cv"),
#                      label = "",
#                      value = cv),
#            textInput(inputId = ns("pa"),
#                      label = "",
#                      value = pa),
#            textInput(inputId = ns("qcf"),
#                      label = "",
#                      value = qcf),
#            textInput(inputId = ns("nbc"),
#                      label = "",
#                      value = nbc),
#            textInput(inputId = ns("cw"),
#                      label = "",
#                      value = cw),
#            textInput(inputId = ns("dr"),
#                      label = "",
#                      value = dr),
#            textInput(inputId = ns("fs"),
#                      label = "",
#                      value = fs),
#            textInput(inputId = ns("cl"),
#                      label = "",
#                      value = cl)
#     ),
#       conditionalPanel(
#         condition = "output.displayDE",
#         ns = ns,
#         panel(
#           heading = "Differential Expression & Marker Selection",
#           h5("Discover quantitative changes between experimental conditions using one of the many integrated statistical frameworks for differential expression:"),
#           actionButton(inputId = ns("goDE"), label = "Go to Differential Expression!"),
#           hr(),
#           h5("Compute and visualize marker genes from cluster classes using 'MAST', 'Limma' or 'DESeq2' methods:"),
#           actionButton(inputId = ns("goMS"), label = "Go to Marker Selection!"),
#         ), br()
#       ),
#     conditionalPanel(
#       condition = "output.displayPA",
#       ns = ns,
#       panel(
#         heading = "Pathway Analysis",
#         h5("Explore biological activity or functions with pathway analysis using either 'GSVA' or 'EnrichR' statistical frameworks:"),
#         actionButton(inputId = ns("goGSVA"), label = "Go to GSVA!"),
#         actionButton(inputId = ns("goEnrichR"), label = "Go to EnrichR!"),
#       ), br()
#     ),
#     conditionalPanel(
#       condition = "output.displayCV",
#       ns = ns,
#       panel(
#         heading = "Cell Viewer",
#         h5("Discover quantitative changes between experimental conditions using one of the many integrated statistical frameworks:"),
#         actionButton(inputId = ns("goCV"), label = "Go to Cell Viewer!")
#       )
#     ),
#     conditionalPanel(
#       condition = "output.displayQCF",
#       ns = ns,
#       panel(
#         heading = "Quality Control & Filtering",
#         h5("Perform quality control checks / filtering."),
#         actionButton(inputId = ns("goQC"), label = "Go to Quality Control!")
#       )
#     ),
#     conditionalPanel(
#       condition = "output.displayNBC",
#       ns = ns,
#       panel(
#         heading = "Normalization & Batch-Correction",
#         h5("go nbc"),
#         actionButton(inputId = ns("goNBC"), label = "Go to Normalization/Batch-Correction!")
#       )
#     ),
#     conditionalPanel(
#       condition = "output.displayCW",
#       ns = ns,
#       panel(
#         heading = "Curated Workflows",
#         h5("go cw"),
#         actionButton(inputId = ns("goCelda"), label = "Go to Celda!"),
#         actionButton(inputId = ns("goSeurat"), label = "Go to Seurat!")
#       )
#     ),
#     conditionalPanel(
#       condition = "output.displayDR",
#       ns = ns,
#       panel(
#         heading = "Dimensionality Reduction",
#         h5("go dr"),
#         actionButton(inputId = ns("goDR"), label = "Go to Dimensionality Reduction!")
#       )
#     ),
#     conditionalPanel(
#       condition = "output.displayFS",
#       ns = ns,
#       panel(
#         heading = "Feature Selection",
#         h5("go fs"),
#         actionButton(inputId = ns("goFS"), label = "Go to Feature Selection!")
#       )
#     ),
#     conditionalPanel(
#       condition = "output.displayCL",
#       ns = ns,
#       panel(
#         heading = "Clustering",
#         h5("go cl"),
#         actionButton(inputId = ns("goCL"), label = "Go to Clustering!")
#       )
#     )
#   )
# }
# 
# #server
# nonLinearWorkflow <- function(input, output, session, parent)
# {
#   ns <- session$ns
# 
#   observeEvent(input$goDE,{
#     showTab(inputId = "navbar",
#             target = "Differential Expression",
#             select = TRUE,
#             session = parent)
#   })
# 
#   observeEvent(input$goMS,{
#     showTab(inputId = "navbar",
#             target = "Find Marker",
#             select = TRUE,
#             session = parent)
#   })
# 
#   observeEvent(input$goGSVA,{
#     showTab(inputId = "navbar",
#             target = "GSVA",
#             select = TRUE,
#             session = parent)
#   })
# 
#   observeEvent(input$goEnrichR,{
#     showTab(inputId = "navbar",
#             target = "EnrichR",
#             select = TRUE,
#             session = parent)
#   })
#   
#   observeEvent(input$goQC,{
#     showTab(inputId = "navbar",
#             target = "QC & Filtering",
#             select = TRUE,
#             session = parent)
#   })
#   
#   observeEvent(input$goNBC,{
#     showTab(inputId = "navbar",
#             target = "Normalization & Batch Correction",
#             select = TRUE,
#             session = parent)
#   })
#   
#   observeEvent(input$goCelda,{
#     showTab(inputId = "navbar",
#             target = "Celda",
#             select = TRUE,
#             session = parent)
#   })
#   
#   observeEvent(input$goSeurat,{
#     showTab(inputId = "navbar",
#             target = "Seurat",
#             select = TRUE,
#             session = parent)
#   })
#   
#   observeEvent(input$goDR,{
#     showTab(inputId = "navbar",
#             target = "DR",
#             select = TRUE,
#             session = parent)
#   })
#   
#   observeEvent(input$goFS,{
#     showTab(inputId = "navbar",
#             target = "FS",
#             select = TRUE,
#             session = parent)
#   })
#   
#   observeEvent(input$goCL,{
#     showTab(inputId = "navbar",
#             target = "CL",
#             select = TRUE,
#             session = parent)
#   })
# 
#     output$displayDE <- reactive({
#       returnedValue = input$de
#       return(as.logical(returnedValue))
#     })
# 
#     output$displayPA <- reactive({
#       returnedValue = input$pa
#       return(as.logical(returnedValue))
#     })
# 
#     output$displayCV <- reactive({
#       returnedValue = input$cv
#       return(as.logical(returnedValue))
#     })
#     
#     output$displayQCF <- reactive({
#       returnedValue = input$qcf
#       return(as.logical(returnedValue))
#     })
#     
#     output$displayNBC <- reactive({
#       returnedValue = input$nbc
#       return(as.logical(returnedValue))
#     })
#     
#     output$displayCW <- reactive({
#       returnedValue = input$cw
#       return(as.logical(returnedValue))
#     })
#     
#     output$displayDR <- reactive({
#       returnedValue = input$dr
#       return(as.logical(returnedValue))
#     })
#     
#     output$displayFS <- reactive({
#       returnedValue = input$fs
#       return(as.logical(returnedValue))
#     })
#     
#     output$displayCL <- reactive({
#       returnedValue = input$cl
#       return(as.logical(returnedValue))
#     })
# 
#   outputOptions(output, "displayDE", suspendWhenHidden = FALSE)
#   outputOptions(output, "displayPA", suspendWhenHidden = FALSE)
#   outputOptions(output, "displayCV", suspendWhenHidden = FALSE)
#   outputOptions(output, "displayQCF", suspendWhenHidden = FALSE)
#   outputOptions(output, "displayNBC", suspendWhenHidden = FALSE)
#   outputOptions(output, "displayCW", suspendWhenHidden = FALSE)
#   outputOptions(output, "displayDR", suspendWhenHidden = FALSE)
#   outputOptions(output, "displayFS", suspendWhenHidden = FALSE)
#   outputOptions(output, "displayCL", suspendWhenHidden = FALSE)
# }


##################################################
##################################################

#UI
nonLinearWorkflowUI <- function(id)
{
  ns <- NS(id)
  fluidPage(
    uiOutput(ns("ui"))
  )
}

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