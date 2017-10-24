shiny_panel_subsample <- fluidPage(
  tabsetPanel(
    tabPanel(
      "Sequencing Depth",
      fluidRow(
        column(
          4,
          wellPanel(
            numericInput("minCount", "Minimum readcount to detect gene",
                         value = 10, min=1, max=10000),
            numericInput("minCells", "Minimum number of cells with nonzero expression to detect gene",
                         value=3, min=1, max=100000),
            #TO DO: Add a 'track a particular gene' option. Show counts for that gene.
            numericInput("iterations",
                         label = "Number of bootstrap iterations.",
                         value = 10, min = 2, max = 10000),
            sliderInput("maxDepth", "Maximum log10(number of simulated reads)", 3, 12, 5, 0.5),
            sliderInput("depthResolution", "how many values to simulate", 5, 100, 10, 5),
            actionButton("runSubsampleDepth", "Run subsampler")
          )
          tabPanel("Genes Detected", plotOutput("DepthDone")),
          tabPanel("Minimum Detectable Effect Size", plotOutput("MinEffectDone"))
          )
        )
      )
    )
    
    tabPanel(
      "# "
    )
  )
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
                             c("tpm.t", "DESeq2")),
                 numericInput("selectTotReads", label = "Total simulated reads (all cells)",
                              value = 1000000, min = 100, max = 1000000000),
                 numericInput("selectNCells", label = "Number of simulated cells",
                              value = 100, min = 10, max = 100000),
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
  ),
  includeHTML("www/footer.html")
)
