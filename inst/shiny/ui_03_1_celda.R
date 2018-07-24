shinyPanelCelda <- fluidPage(
  tags$div(
    class = "container",
    h1("Celda"),
    h5(tags$a(href = "https://github.com/compbiomed/celda",
              "(help)", target = "_blank")),
    mainPanel(
      actionButton(inputId = "runCelda", label = "Run celda"),
      plotOutput("celdaPlot", height = "600px")
    )
  )
)
