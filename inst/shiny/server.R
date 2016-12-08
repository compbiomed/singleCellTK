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
    updateSelectInput(session, "subCovariate", choices = colnames(pData(vals$counts)))
    insertUI(
      selector = '#uploadAlert',
      ## wrap element in a div with id for ease of removal
      ui = tags$div(class="alert alert-success", "Successfully Uploaded!")
      )
  })

  output$contents <- renderDataTable({
    if(!(is.null(vals$counts))){
      exprs(vals$counts)
    }
  }, options = list(scrollX = TRUE))

  output$summarycontents <- renderTable({
    if(!(is.null(vals$counts))){
      summarizeTable(exprs(vals$counts))
    }
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
  
  runDownsampler <- observeEvent(input$runSubsample, {
    subData <- reactiveValues(
      counts=Downsample(counts(vals$counts), newcounts=floor(2^seq.int(from=log2(input$minSim), to=log2(input$maxSim), length.out=10)), iterations=input$iterations)
    )
    output$downDone <- renderPlot({
      heatmap(as.matrix(subData$counts[order(apply(subData$counts[,,10,1],1,sum),decreasing=TRUE)[1:20],,10,1]))
    })
  })
  
  runDiffPower <- observeEvent(input$runDifferentialPower, {
#    if(exists('subData$counts')){
      output$powerBoxPlot <- renderPlot({
        subData <- Downsample(counts(vals$counts), newcounts=floor(2^seq.int(from=log2(input$minSim), to=log2(input$maxSim), length.out=10)), iterations=input$iterations)
        diffPower <- differentialPower(datamatrix=counts(vals$counts), downmatrix=subData, conditions=phenoData(vals$counts)[[input$subCovariate]], method=input$selectDiffMethod)
        boxplot(diffPower)#,names=floor(2^seq.int(from=log2(input$minSim), to=log2(input$maxSim), length.out=10)))
      })
#    }
#    else{
#      output$powerBoxPlot <- renderPlot({
#        plot(c(0,1),c(0,1),main="You need to run the subsampler first.")
#      })
#    }
  })
})
