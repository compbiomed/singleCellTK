library(shiny)
library(scater)
library(ComplexHeatmap)
library(biomaRt)
library(circlize)
library(limma)
library(d3heatmap)
library(ggplot2)
library(plotly)
library(DESeq)

#1GB max upload size
options(shiny.maxRequestSize=1000*1024^2)

# Define server logic required to draw a histogram
shinyServer(function(input, output, session) {

  vals <- reactiveValues(
    counts = getShinyOption("inputSCEset"),
    original = getShinyOption("inputSCEset")
  )

  observeEvent(input$uploadData, {
    vals$counts <- createSCESet(input$countsfile$datapath,
                                input$annotfile$datapath)
    updateSelectInput(session, "colorClusters",
                      choices = colnames(pData(vals$counts)))
    updateSelectInput(session, "selectDiffex_condition",
                      choices = colnames(pData(vals$counts)))
    updateSelectInput(session, "subCovariate",
                      choices = colnames(pData(vals$counts)))
    insertUI(
      selector = '#uploadAlert',
      ## wrap element in a div with id for ease of removal
      ui = tags$div(class="alert alert-success", "Successfully Uploaded!")
      )
    vals$original <- vals$counts
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
    nkeeprows <- ceiling((1-(0.01 * input$LowExpression)) * as.numeric(nrow(vals$counts)))
    tokeep <- order(rowSums(counts(vals$counts)), decreasing = TRUE)[1:nkeeprows]
    vals$counts <- vals$counts[tokeep,]
  })
  
  observeEvent(input$resetData, {
    vals$counts <- vals$original
  })

  clusterDataframe <- observeEvent(input$clusterData, {
    if(input$selectCustering == "PCA"){
      pca <- scater::plotPCA(vals$counts, return_SCESet=TRUE)
      g <- reducedDimension(pca)
      l <- data.frame(g)
      w <- input$colorClusters
      l$Treatment <- eval(parse(text = paste("pData(vals$counts)$",w,sep="")))
      g <- ggplot(l, aes(PC1, PC2, color=Treatment))+geom_point()
      output$clusterPlot <- renderPlotly({
        ggplotly(g)
        #plotPCA(vals$counts, colour_by=input$colorClusters)
      })
    } else if(input$selectCustering == "tSNE"){
      output$clusterPlot <- renderPlot({
        scater::plotTSNE(vals$counts, colour_by=input$colorClusters)
      })
    }
  })
  
  deHeatmapDataframe <- observeEvent(input$makeHeatmap, {
    if(input$selectHeatmap == "Standard") {
      output$heatmapPlot <- renderPlot({
        heatmap(counts(vals$counts)[1:50,], labCol = FALSE, labRow = FALSE)})
    } else if(input$selectHeatmap == "Complex") {
      # Do Something
    } else if(input$selectHeatmap == "Interactive") {
      output$heatmapPlot <- renderD3heatmap({
        plotHeatmap(vals$counts)})
    }
  })
  
  diffexDataframe <- observeEvent(input$runDiffex, {
    output$diffPlot <- renderPlot({
      scDiffEx(vals$counts, input$selectDiffex_condition, input$selectPval,
               input$selectNGenes, input$applyCutoff)
    }, height=600)
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
