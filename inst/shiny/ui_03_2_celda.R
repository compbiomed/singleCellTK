shinyPanelCelda <- fluidPage(
  tags$div(
    class = "container",
    h1("celda: CEllular Latent Dirichlet Allocation"),
    h5(tags$a(href = "https://github.com/compbiomed/celda",
              "(help)", target = "_blank")),
    fluidRow(
      sidebarLayout(
        sidebarPanel(
          # SHINYJS ACCORDION --------------------------
          # Section 1 - Basic Settings
          actionButton("celdaBasicSet", "Basic Settings"),
          # open by default
          tags$div(id = "celdaCollapse1",
            wellPanel(
              selectInput("celdaAssay", "Select Assay:", currassays),
              selectInput("celdaModel",
                "Select Celda Model:",
                c("celda_C", "celda_G", "celda_CG"),
                selected = "celda_CG"),
              # c("celda_C"),
              # selected = "celda_C"),
              conditionalPanel(
                condition = sprintf("input['%s'] == 'celda_C'", "celdaModel"),
                numericInput("cellClusterC",
                  label = "Number of Cell Clusters (K):",
                  value = 2,
                  min = 1,
                  max = 100000,
                  step = 1)
              ),
              conditionalPanel(
                condition = sprintf("input['%s'] == 'celda_G'", "celdaModel"),
                numericInput("geneModuleG",
                  label = "Number of Gene Modules (L):",
                  value = 2,
                  min = 1,
                  max = 100000,
                  step = 1)
              ),
              conditionalPanel(
                condition = sprintf("input['%s'] == 'celda_CG'", "celdaModel"),
                numericInput("cellClusterCG",
                  label = "Number of Cell Clusters (K):",
                  value = 2,
                  min = 1,
                  max = 100000,
                  step = 1),
                numericInput("geneModuleCG",
                  label = "Number of Gene Modules (L):",
                  value = 2,
                  min = 1,
                  max = 100000,
                  step = 1))
            )
          ),
          # Section 2 - Advanced Settings
          actionButton("celdaAdvSet", "Advanced Settings"),
          shinyjs::hidden(
            tags$div(id = "celdaCollapse2",
              wellPanel(
                conditionalPanel(
                  condition = sprintf("input['%s'] == 'celda_C' || input['%s'] ==
                    'celda_CG'", "celdaModel", "celdaModel"),
                  numericInput("celdaAlpha",
                    label = "Alpha",
                    value = 1,
                    min = 0.00000001,
                    max = 100000)
                  ),
                numericInput("celdaBeta",
                  label = "Beta",
                  value = 1,
                  min = 0.00000001,
                  max = 100000),
                conditionalPanel(
                  condition = sprintf("input['%s'] == 'celda_G' || input['%s'] ==
                    'celda_CG'", "celdaModel", "celdaModel"),
                  numericInput("celdaDelta",
                    label = "Delta:",
                    value = 1,
                    min = 0.00000001,
                    max = 100000),
                  numericInput("celdaGamma",
                    label = "Gamma:",
                    value = 1,
                    min = 0.00000001,
                    max = 100000)
                  ),
                numericInput("celdaMaxIter",
                  label = "Maximum Number of Iterations:",
                  value = 200,
                  min = 1,
                  max = 100000,
                  step = 1),
                numericInput("celdaStopIter",
                  label = "Number of Converging Iterations for Gibbs sampler to
                  stop:",
                  value = 10,
                  min = 1,
                  max = 100000,
                  step = 1),
                numericInput("celdaSplitIter",
                  label = "Split on Every This Number of Iteration:",
                  value = 10,
                  min = 1,
                  max = 100000,
                  step = 1),
                numericInput("celdaNChains",
                  label = "Number of Gibbs Sampling chains to run for every K/L
                  combination:",
                  value = 1,
                  min = 1,
                  max = 100000,
                  step = 1),
                numericInput("celdaCores",
                  label = "Number of Cores to use for parallel Gibbs Sampling:",
                  value = 1,
                  min = 1,
                  max = 100000,
                  step = 1),
                numericInput("celdaSeed",
                  label = "Base Seed For Random Number Generation:",
                  value = 12345,
                  min = 1,
                  max = 100000,
                  step = 1)
                )
              )
            ),
          tags$hr(),
          withBusyIndicatorUI(actionButton(inputId = "runCelda",
            label = "Run celda")),
          tags$hr(),
          downloadButton("downloadSCECelda", "Download SCtkExperiment")
        ),
        mainPanel(
          conditionalPanel(
            condition = sprintf("input['%s'] == 'celda_C' || input['%s'] ==
              'celda_G' || input['%s'] == 'celda_CG'",
              "celdaModel",
              "celdaModel",
              "celdaModel"),
            plotOutput("celdaPlot", height = "600px")
            )
          )
        )
      )
    )
  )
