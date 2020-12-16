#UI
nonLinearWorkflowUI <- function(id, de = FALSE, cv = FALSE, pa = FALSE)
{
  ns <- NS(id)
  fluidPage(
    hidden(textInput(inputId = ns("de"),
                     label = "",
                     value = de),
           textInput(inputId = ns("cv"),
                     label = "",
                     value = cv),
           textInput(inputId = ns("pa"),
                     label = "",
                     value = pa)
    ),
      conditionalPanel(
        condition = "output.displayDE",
        ns = ns,
        panel(
          heading = "Differential Expression & Marker Selection",
          h5("Discover quantitative changes between experimental conditions using one of the many integrated statistical frameworks for differential expression:"),
          actionButton(inputId = ns("goDE"), label = "Go to Differential Expression!"),
          hr(),
          h5("Compute and visualize marker genes from cluster classes using 'MAST', 'Limma' or 'DESeq2' methods:"),
          actionButton(inputId = ns("goMS"), label = "Go to Marker Selection!"),
        ), br()
      ),
    conditionalPanel(
      condition = "output.displayPA",
      ns = ns,
      panel(
        heading = "Pathway Analysis",
        h5("Explore biological activity or functions with pathway analysis using either 'GSVA' or 'EnrichR' statistical frameworks:"),
        actionButton(inputId = ns("goGSVA"), label = "Go to GSVA!"),
        actionButton(inputId = ns("goEnrichR"), label = "Go to EnrichR!"),
      ), br()
    ),
    conditionalPanel(
      condition = "output.displayCV",
      ns = ns,
      panel(
        heading = "Cell Viewer",
        h5("Discover quantitative changes between experimental conditions using one of the many integrated statistical frameworks:"),
        actionButton(inputId = ns("goCV"), label = "Go to Cell Viewer!")
      )
    )
  )
}

#server
nonLinearWorkflow <- function(input, output, session, parent)
{
  ns <- session$ns
  
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
  
    output$displayDE <- reactive({
      returnedValue = input$de
      return(as.logical(returnedValue))
    })
    
    output$displayPA <- reactive({
      returnedValue = input$pa
      return(as.logical(returnedValue))
    })
    
    output$displayCV <- reactive({
      returnedValue = input$cv
      return(as.logical(returnedValue))
    })
    
  outputOptions(output, "displayDE", suspendWhenHidden = FALSE)
  outputOptions(output, "displayPA", suspendWhenHidden = FALSE)
  outputOptions(output, "displayCV", suspendWhenHidden = FALSE)
}