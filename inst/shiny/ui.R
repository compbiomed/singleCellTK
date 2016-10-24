library(shiny)

# Define UI for application that draws a histogram
shinyUI(
  navbarPage(
    "Single Cell Toolkit (alpha)",
    #bootstrap theme
    theme = "bootstrap.min.css",

    #Upload Tab
    tabPanel(
      "Upload",
      tags$div(
        class="jumbotron",
        tags$div(
          class="container",
          h1("Single Cell Toolkit"),
          p("Filter, cluster, and analyze single cell RNA-Seq data")
        )
      ),
      tags$div(
        class="container",
        fileInput('countsfile', 'Choose file to upload',
                  accept = c(
                    'text/csv',
                    'text/comma-separated-values',
                    'text/tab-separated-values',
                    'text/plain',
                    '.csv',
                    '.tsv'
                  )
        )
      ),
      includeHTML('www/footer.html')
    ),
    tabPanel(
      "Data Summary",
      tags$div(
        class="container",
        h1("Data Summary"),
        fluidPage(
          fluidRow(
            column(8, tableOutput('summarycontents')),
            column(
              4,
              wellPanel(
                numericInput('minDetectGenect', label = 'Minimum Detected Genes per Sample.', value=1700, min = 1, max = 100000),
                numericInput("LowExpression", "% Low Gene Expression to Filter",value=40, min = 0, max = 100),
                actionButton("filterData", "Filter Data"),
                actionButton("filterData", "Reset")
              )
            )
          ),
          fluidRow(
            dataTableOutput('contents')
          )
        )
      ),
      includeHTML('www/footer.html')
    ),
    tabPanel(
      "Clustering",
      tags$div(
        class="container",
        h1("Clustering"),
        fluidPage(
          fluidRow(
            column(4,
                   wellPanel(
                     selectInput("selectCustering","Clustering Algorithm",c("PCA","tSNE")),
                     actionButton("clusterData", "Cluster Data")
                   )),
            column(8,
                   plotOutput("clusterPlot"))
          )
        )
      ),
      includeHTML('www/footer.html')
    ),
    tabPanel(
      "Batch Correction",
      tags$div(
        class="container",
        h1("Batch Correction")
      ),
      includeHTML('www/footer.html')
    ),
    tabPanel(
      "Pathway",
      tags$div(
        class="container",
        h1("Pathway Profiling")
      ),
      includeHTML('www/footer.html')
    ),

    #Data Summary Tab
    tabPanel(
      "Example Tab",
      tags$div(
        class="container",
        fluidPage(
          titlePanel("My Shiny App"),
          sidebarLayout(
            sidebarPanel(
              h2("Installation"),
              p("Shiny is available on CRAN, so you can install it in the usual way from your R console:"),
              code('install.packages("shiny")'),
              br(),
              br(),
              br(),
              br(),
              "shiny is a product of ",
              span("RStudio", style = "color:blue"),
              br(),
              sliderInput("bins", "Number of bins:", min = 1, max = 50, value = 30)
            ),
            mainPanel(
              h1("Introducing Shiny"),
              p("Shiny is a new package from RStudio that makes it ",
                em("incredibly easy"),
                " to build interactive web applications with R."),
              br(),
              p("For an introduction and live examples, visit the ", a("Shiny homepage.", href = "http://www.rstudio.com/shiny")),
              br(),
              h2("Features"),
              p("* Build useful web applications with only a few lines of code—no JavaScript required."),
              p("* Shiny applications are automatically “live” in the same way that ",
                strong("spreadsheets"),
                " are live. Outputs change instantly as users modify inputs, without requiring a reload of the browser."),
              plotOutput("distPlot")
            )
          )
        )
      ),
      includeHTML('www/footer.html')
    ),
    navbarMenu(
      "More",
      tabPanel(
        "Sub-Component A",
        includeHTML('www/footer.html')),
      tabPanel(
        "Sub-Component B",
        includeHTML('www/footer.html'))
    )
  )
)
