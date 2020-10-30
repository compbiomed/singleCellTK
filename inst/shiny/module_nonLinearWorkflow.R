#UI
nonLinearWorkflowUI <- function(id, de = FALSE, ms = FALSE, cv = FALSE)
{
  ns <- NS(id)
  fluidPage(
    hidden(textInput(inputId = ns("de"),
                     label = "",
                     value = de),
           textInput(inputId = ns("ms"),
                     label = "",
                     value = ms),
           textInput(inputId = ns("cv"),
                     label = "",
                     value = cv)
    ),
      conditionalPanel(
        condition = "output.displayDE",
        ns = ns,
        panel(
          heading = "Differential Expression",
          h5("Discover quantitative changes between experimental conditions using one of the many integrated statistical frameworks:"),
          actionButton(inputId = ns("goDE"), label = "Go to Differential Expression!")
        ), br()
      ),
    conditionalPanel(
      condition = "output.displayMS",
      ns = ns,
      panel(
        heading = "Marker Selection",
        h5("Discover quantitative changes between experimental conditions using one of the many integrated statistical frameworks:"),
        actionButton(inputId = ns("goMS"), label = "Go to Marker Selection!")
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
  
    output$displayDE <- reactive({
      returnedValue = input$de
      return(as.logical(returnedValue))
    })
    
    output$displayMS <- reactive({
      returnedValue = input$ms
      return(as.logical(returnedValue))
    })
    
    output$displayCV <- reactive({
      returnedValue = input$cv
      return(as.logical(returnedValue))
    })
    
  outputOptions(output, "displayDE", suspendWhenHidden = FALSE)
  outputOptions(output, "displayMS", suspendWhenHidden = FALSE)
  outputOptions(output, "displayCV", suspendWhenHidden = FALSE)
}