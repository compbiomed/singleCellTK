library(shiny)
library(scater)

#1GB max upload size
options(shiny.maxRequestSize=1000*1024^2)

#move this to the package itself?
summarizeTable <- function(indata){
  return(data.frame("Metric"=c("Number of Samples","Number of Genes", "Samples with <1700 detected genes"),
                    "Value"=c(ncol(indata),
                              nrow(indata),
                              sum(apply(indata, 2, function(x) sum(as.numeric(x)==0)) < 1700))))
}

createSCESet <- function(countfile, annotfile){
  countsin <- read.table(countfile, sep="\t", header=T, row.names=1)
  annotin <- read.table(annotfile, sep="\t", header=T, row.names=1)
  pd <- new("AnnotatedDataFrame", data = annotin)

  gene_df <- data.frame(Gene = rownames(countsin))
  rownames(gene_df) <- gene_df$Gene
  fd <- new("AnnotatedDataFrame", data = gene_df)
  return(newSCESet(countData = countsin, phenoData = pd,
                   featureData = fd))
}

# Define server logic required to draw a histogram
shinyServer(function(input, output, session) {

  vals <- reactiveValues(
    counts = getShinyOption("inputSCEset")
  )

  observeEvent(input$uploadData, {
    vals$counts <- createSCESet(input$countsfile$datapath, input$annotfile$datapath)
    updateSelectInput(session, "colorClusters", choices = colnames(pData(vals$counts)))
  })

  output$contents <- renderDataTable({
    exprs(vals$counts)
  }, options = list(scrollX = TRUE))

  output$summarycontents <- renderTable({
    summarizeTable(exprs(vals$counts))
  })

  observeEvent(input$filterData, {
    vals$counts <- vals$counts[1:100,]
  })

  clusterDataframe <- observeEvent(input$clusterData, {
    if(input$selectCustering == "PCA"){
      output$clusterPlot <- renderPlot({
        plotPCA(vals$counts, colour_by=input$colorClusters)
      })
    } else if(input$selectCustering == "tSNE"){
      output$clusterPlot <- renderPlot({
        plotTSNE(vals$counts, colour_by=input$colorClusters)
      })
    }
  })
})
