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
        actionButton("primaryButton", "Upload Count Matrix", class="btn btn-primary btn-lg")
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
    tabPanel(
      "Data Summary",
      h1("Data Summary"),
      includeHTML('www/footer.html')
    ),
    tabPanel(
      "Clustering",
      h1("Clustering"),
      includeHTML('www/footer.html')
    ),
    tabPanel(
      "Batch Correction",
      h1("Batch Correction"),
      includeHTML('www/footer.html')
    ),
    tabPanel(
      "Pathway",
      h1("Pathway Profiling"),
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
