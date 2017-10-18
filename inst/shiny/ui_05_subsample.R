shiny_panel_subsample <- fluidPage(
  tags$div(
    class = "container",
    h1("Subsampling"),
    fluidPage(
      fluidRow(
        column(4,
               wellPanel(
                 selectInput("subCovariate",
                             "Covariate for differential expression",
                             clusterChoice),
                 selectInput("selectDiffMethod",
                             "Differential Expression Method",
                             c("tpm.t", "DESeq")),
                 numericInput("minSim", label = "Minimum subsample Size.",
                              value = 1000, min = 1, max = 1000000),
                 numericInput("maxSim", label = "Maximum subsample Size.",
                              value = 10000, min = 100, max = 100000000),
                 numericInput("iterations",
                              label = "Number of bootstrap iterations.",
                              value = 10, min = 2, max = 1000),
                 actionButton("runSubsample", "Run subsampler"),
                 actionButton("runDifferentialPower",
                              "Run differential power analysis")
               )),
        column(8,
               plotOutput("downDone"))
      ),
      fluidRow(
        plotOutput("powerBoxPlot")
      )
    )
  )
)
