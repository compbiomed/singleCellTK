library(shiny)

# Define server logic required to draw a histogram
shinyServer(function(input, output) {

  # Expression that generates a histogram. The expression is
  # wrapped in a call to renderPlot to indicate that:
  #
  #  1) It is "reactive" and therefore should re-execute automatically
  #     when inputs change
  #  2) Its output type is a plot

  output$distPlot <- renderPlot({
    x    <- faithful[, 2]  # Old Faithful Geyser data
    bins <- seq(min(x), max(x), length.out = input$bins + 1)

    # draw the histogram with the specified number of bins
    hist(x, breaks = bins, col = 'darkgray', border = 'white')
  })

  vals <- reactiveValues(
    counts = NULL
  )

  observeEvent(input$countsfile, {
    vals$counts <- read.table(input$countsfile$datapath, header = TRUE, sep = "\t")
  })

  summarizeTable <- function(indata){
    return(data.frame("Metric"=c("Number of Samples","Number of Genes", "Samples with <1700 detected genes"),
                      "Value"=c(ncol(indata)-1,
                                nrow(indata),
                                sum(apply(indata[,2:ncol(indata)], 2, function(x) sum(as.numeric(x)==0)) < 1700))))
  }

  output$contents <- renderDataTable({
    vals$counts
  }, options = list(scrollX = TRUE))

  output$summarycontents <- renderTable({
    summarizeTable(vals$counts)
  })

  observeEvent(input$filterData, {
    vals$counts <- vals$counts[1:100,]
  })

  #we need to filter 0s
  clusterDataframe <- eventReactive(input$clusterData, {
    prcomp(t(vals$counts[,2:ncol(vals$counts)]), center=T, scale. = T)
  })

  output$clusterPlot <- renderPlot({
    plot(clusterDataframe()$x[,1],clusterDataframe()$x[,2])
  })
})
