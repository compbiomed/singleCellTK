cxdsHelpModal <- function() {
  modalDialog(
    tags$style(HTML("
      div {
        word-wrap: break-word;
      }
      ")),
    tags$div(
      h3("CXDS Parameters"),
      fluidRow(
        column(4, tags$b("Parameter Name")),
        column(8, tags$b("Description"))
      ),
      tags$hr(),
      fluidRow(
        column(4, "ntop"),
        column(8, "Integer. Indimessageing number of top variance genes to consider. Default: 500")
      ),
      tags$hr(),
      fluidRow(
        column(4, "binThresh"),
        column(8, 'Integer. Minimum counts to consider a gene "present" in a cell. Default: 0')
      ),
      tags$hr(),
      fluidRow(
        column(4, "verb"),
        column(8, "Check off to print log messages to the console. Default: TRUE.")
      ),
      tags$hr(),
      fluidRow(
        column(4, "retRes"),
        column(8, "Check off to return gene pair scores & top-scoring gene pairs.")
      ),
      tags$hr(),
      fluidRow(
        column(4, "estNdbl"),
        column(8, "Check off to estimate  the numer of doublets be estimated from the data. Enables doublet calls.")
      ),
    )
  )
}

