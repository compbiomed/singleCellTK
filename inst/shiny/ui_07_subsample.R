shinyPanelSubsample <- fluidPage(
  tags$h1("Sample Size Calculator"),
  h5(tags$a(href = "https://compbiomed.github.io/sctk_docs/articles/v09-tab07_Sample-Size.html",
            "(help)", target = "_blank")),
  tabsetPanel(
    tabPanel(
      "Sequencing Depth",
      fluidRow(
        column(
          4,
          wellPanel(
            selectInput("depthAssay", "Select Assay:", currassays),
            numericInput("minCount", "Minimum readcount to detect gene",
                         value = 10, min = 1, max = 10000),
            numericInput("minCells", "Minimum number of cells with nonzero expression to detect gene",
                         value = 3, min = 1, max = 100000),
            #TO DO: Add a 'track a particular gene' option. Show counts for that gene.
            numericInput("iterations",
                         label = "Number of bootstrap iterations.",
                         value = 10, min = 2, max = 10000),
            sliderInput("maxDepth", "Maximum log10(number of simulated reads)", 3, 12, 5, 0.5),
            sliderInput("depthResolution", "Number of values to simulate", 5, 100, 10, 5),
            selectInput("selectReadDepthCondition", "Condition for diffex", c("Random", clusterChoice), selected = "Random"),
            withBusyIndicatorUI(actionButton("runSubsampleDepth", "Run subsampler"))
          )
        ),
        column(8,
          tabsetPanel(
            tabPanel("Genes Detected", plotOutput("depthDone")),
            tabPanel("Minimum Detectable Effect Size", plotOutput("minEffectDone")),
            tabPanel("Number of Diffex Genes", plotOutput("sigNumDone"))
          )
        )
      )
    ),
    tabPanel(
      "Number of cells",
      fluidRow(
        column(
          4,
          wellPanel(
            selectInput("cellsAssay", "Select Assay:", currassays),
            numericInput("minCellNum", "Minimum number of cells to simulate",
                         value = 10, min = 1, max = 10000),
            numericInput("maxCellNum", "Maximum number of cells to simulate",
                         value = 100, min = 10, max = 100000),
            numericInput("iterations",
                         label = "Number of bootstrap iterations per cellcount.",
                         value = 10, min = 2, max = 10000),
            numericInput("totalReads", "Estimated number of aligned reads",
                         value = 1000000, min = 1000, max = 1000000000),
            checkboxInput("useReadCount", "Use the same number of reads as in original dataset"),
            selectInput("selectCellNumCondition", "Condition for diffex", clusterChoice),
            numericInput("minCount", "Minimum readcount to detect gene",
                         value = 10, min = 1, max = 10000),
            numericInput("minCells", "Minimum number of cells with nonzero expression to detect gene",
                         value = 3, min = 1, max = 100000),
            numericInput("depthResolution", "Number of dataset sizes to simulate",
                         value = 10, min = 1, max = 100),
            withBusyIndicatorUI(actionButton("runSubsampleCells", "Run resampler"))
          ),
          tabPanel("Genes Detected", plotOutput("CellCountDone")),
          tabPanel("Minimum Detectable Effect Size", plotOutput("MinEffectCellsDone")),
          tabPanel("Number of Significant DiffEx Genes", plotOutput("cellNumSigDone"))
        ),
        column(8,
          tabsetPanel(
            tabPanel("Genes Detected", plotOutput("cellsDone")),
            tabPanel("Minimum Detectable Effect Size", plotOutput("minEffectCells")),
            tabPanel("Number of Diffex Genes", plotOutput("sigNumCells"))
          )
        )
      )
    ),
    tabPanel(
      "Snapshot",
      fluidRow(
        column(
          4,
          wellPanel(
            selectInput("snapshotAssay", "Select Assay:", currassays),
            numericInput("numCellsSnap", "How many simulated cells?",
                         value = 100, min = 2, max = 10000),
            numericInput("numReadsSnap", "How many total reads?",
                         value = 1000000, min = 1000, max = 1000000000),
            selectInput("selectSnapshotCondition", "Condition for diffex", clusterChoice),
            numericInput("iterationsSnap", "Number of bootstrap iterations",
                         value = 10, min = 2, max = 1000),
            withBusyIndicatorUI(actionButton("runSnapshot", "Run resampling snapshot"))
          )
        ),
        column(8,
          plotOutput("Snaplot")
        )
      )
    )
  )
)

