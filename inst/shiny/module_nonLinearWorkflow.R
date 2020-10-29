#UI
nonLinearWorkflowUI <- function(id)
{
  ns <- NS(id)
  fluidPage(
    fluidRow(
      column(1),
      column(5,
             fluidRow(
               panel(
                 heading = "Differential Expression",
                 h5("Discover quantitative changes between experimental conditions using one of the many integrated statistical frameworks:"),
                 actionButton(inputId = ns("SeuratDE"), label = "Go DE!")
               )
             )
      ),
      column(5,
             fluidRow(
               panel(
                 heading = "Marker Selection",
                 h5("Differential Expression can be used:"),
                 actionButton(inputId = ns("da2"), label = "Go MS!")
               )
             )),
      column(1)
    )
  )
}

#server
nonLinearWorkflow <- function(input, output, session, parent)
{
  ns <- session$ns
  
  observeEvent(input$SeuratDE,{
    showTab(inputId = "navbar",
            target = "Differential Expression",
            select = TRUE,
            session = parent)
  })
}