library(shiny)
library(shinyjs)
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
    withBusyIndicatorServer("uploadData", {
      vals$counts <- createSCESet(countfile = input$countsfile$datapath,
                                  annotfile = input$annotfile$datapath,
                                  featurefile = input$featurefile$datapath)
      updateSelectInput(session, "colorClusters",
                        choices = colnames(pData(vals$counts)))
      updateSelectInput(session, "deletesamplelist",
                        choices = rownames(pData(vals$counts)))
      updateSelectInput(session, "selectDiffex_condition",
                        choices = colnames(pData(vals$counts)))
      updateSelectInput(session, "subCovariate",
                        choices = colnames(pData(vals$counts)))
      updateSelectInput(session, "pcX",
                        choices = paste("PC",1:nrow(pData(vals$counts)),sep=""),
                        selected = "PC1")
      updateSelectInput(session, "pcY",
                        choices = paste("PC",1:nrow(pData(vals$counts)),sep=""),
                        selected = "PC2")
      insertUI(
        selector = '#uploadAlert',
        ## wrap element in a div with id for ease of removal
        ui = tags$div(class="alert alert-success alert-dismissible", HTML("<span \
                      class='glyphicon glyphicon-ok' aria-hidden='true'></span> \
                      Successfully Uploaded! <button type='button' class='close' \
                      data-dismiss='alert'>&times;</button>"))
        )
      vals$original <- vals$counts
    })
  })

  output$contents <- renderDataTable({
    if(!(is.null(vals$counts)) && nrow(pData(vals$counts)) < 50){
      temptable <- cbind(rownames(fData(vals$counts)),exprs(vals$counts))
      colnames(temptable)[1] <- "Gene"
      temptable
    }
  }, options = list(scrollX = TRUE))

  output$summarycontents <- renderTable({
    if(!(is.null(vals$counts))){
      summarizeTable(vals$counts)
    }
  })

  observeEvent(input$filterData, {
    if(is.null(vals$original)){
      alert("Warning: Upload data first!")
    }
    else{
      vals$counts <- vals$original
      vals$counts <- filterSCData(vals$counts,
                                  deletesamples=input$deletesamplelist,
                                  remove_noexpress=input$removeNoexpress,
                                  remove_bottom=0.01 * input$LowExpression,
                                  minimum_detect_genes=input$minDetectGenect)
    }
  })

  observeEvent(input$resetData, {
    if(is.null(vals$original)){
      alert("Warning: Upload data first!")
    }
    else{
      vals$counts <- vals$original
      updateSelectInput(session, "deletesamplelist",
                        choices = rownames(pData(vals$counts)))
    }
  })

  drDataframe <- observeEvent(input$plotData, {
    withBusyIndicatorServer("plotData", {
      if(is.null(vals$counts)){
        alert("Warning: Upload data first!")
      }
      else{
        g <- runDimRed(input$selectDimRed, vals$counts, input$colorClusters, input$pcX, input$pcY)
        output$dimredPlot <- renderPlotly({
          ggplotly(g)
        })
      }
    })
  })

  # Below needs to be put into a function (partially runDimRed) and/or make a new function or redesign both functions (Emma - 2/16/17)
  clusterDataFrame <- observeEvent(input$plotClusters, {
    if(input$selectCluster == "K-Means" && input$selectDataC == "PCA Components") {
      pc1 <- input$pcX
      pc2 <- input$pcY
      k <- input$selectK
      pca <- scater::plotPCA(vals$counts, return_SCESet=TRUE)
      pca <- data.frame(reducedDimension(pca))
      pca <- setNames(cbind(rownames(pca), pca, row.names=NULL), c("Sample", colnames(pca)))
      w <- input$colorClusters
      d <- c("Treatment")
      pca$Treatment <- eval(parse(text = paste("pData(vals$counts)$",w,sep="")))
      pca$Sample <- rownames(pData(vals$counts))
      cl <- kmeans(pca[,(strtoi(strsplit(pc1, split = "PC")[[1]][2])+1):(strtoi(strsplit(pc2, split = "PC")[[1]][2])+1)], k)
      g <- ggplot(pca, aes(PC1, PC2, label=Sample, color=factor(cl$cluster), shape=Treatment)) + 
        geom_point() +
        theme(legend.title=element_blank())
    } else if(input$selectCluster == "K-Means" && input$selectDataC == "tSNE Components"){
      k <- input$selectK
      tsne <- scater::plotTSNE(vals$counts, return_SCESet=TRUE)
      tsne <- data.frame(reducedDimension(tsne))
      tsne <- setNames(cbind(rownames(tsne), tsne, row.names=NULL), c("Sample", colnames(tsne)))
      w <- input$colorClusters
      tsne$Treatment <- eval(parse(text = paste("pData(vals$counts)$",w,sep="")))
      tsne$Sample <- rownames(pData(vals$counts))
      cl <- kmeans(tsne[,2:3], k)
      g <- ggplot(tsne, aes(X1, X2, label=Sample, color=factor(cl$cluster), shape=Treatment)) + 
        geom_point() +
        theme(legend.title=element_blank())
    } else if(input$selectCluster == "K-Means" && input$selectDataC == "Raw Data"){
      pc1 <- input$pcX
      pc2 <- input$pcY
      k <- input$selectK
      # Need to Transpose vals$counts so samples are clustered
      cl <- kmeans(t(exprs(vals$counts)), k)
      pca <- scater::plotPCA(vals$counts, return_SCESet=TRUE)
      pca <- data.frame(reducedDimension(pca))
      pca <- setNames(cbind(rownames(pca), pca, row.names=NULL), c("Sample", colnames(pca)))
      w <- input$colorClusters
      d <- c("Treatment")
      pca$Treatment <- eval(parse(text = paste("pData(vals$counts)$",w,sep="")))
      pca$Sample <- rownames(pData(vals$counts))
      g <- ggplot(pca, aes(PC1, PC2, label=Sample, color=factor(cl$cluster), shape=Treatment)) + 
        geom_point() +
        theme(legend.title=element_blank())
    }
    output$clusterPlot <- renderPlotly({
      ggplotly(g)
    })

  })
  # END Emma's Note

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
    if(is.null(vals$counts)){
      alert("Warning: Upload data first!")
    }
    else{
      withBusyIndicatorServer("runDiffex", {
        #run diffex to get gene list and pvalues
        vals$diffexgenelist <- scDiffEx(vals$counts, input$selectDiffex_condition,
                                        input$selectPval, input$selectNGenes, input$applyCutoff,
                                        diffexmethod=input$selectDiffex,
                                        clusterRow=input$clusterRows,
                                        clusterCol=input$clusterColumns)
      })
    }
  })

  output$diffPlot <- renderPlot({
    if(!is.null(vals$diffexgenelist)){
      draw(plot_DiffEx(vals$counts, input$selectDiffex_condition,
                  rownames(vals$diffexgenelist), clusterRow=input$clusterRows,
                  clusterCol=input$clusterColumns,
                  displayRowLabels=input$displayHeatmapRowLabels,
                  displayColumnLabels=input$displayHeatmapColumnLabels,
                  displayRowDendrograms=input$displayHeatmapRowDendrograms,
                  displayColumnDendrograms=input$displayHeatmapColumnDendrograms))
    }
  }, height=600)

  output$interactivediffPlot <- renderD3heatmap({
    if(!is.null(vals$diffexgenelist)){
      plot_d3DiffEx(vals$counts, input$selectDiffex_condition,
                    rownames(vals$diffexgenelist), clusterRow=input$clusterRows,
                    clusterCol=input$clusterColumns)
    }
  })

  output$diffextable <- renderDataTable({
    if(!is.null(vals$diffexgenelist)){
      temptable <- cbind(rownames(vals$diffexgenelist),data.frame(vals$diffexgenelist))
      colnames(temptable)[1] <- "Gene"
      temptable
    }
  })

  # Need to modify the scDiffEx function to return gene list initially and then
  # returns
  output$downloadGeneList <- downloadHandler(
    filename = function() {
      paste("diffex_results-", Sys.Date(), ".csv", sep="")
    },
    content = function(file) {
      write.csv(vals$diffexgenelist, file)
    }
  )

  runDownsampler <- observeEvent(input$runSubsample, {
    if(is.null(vals$counts)){
      alert("Warning: Upload data first!")
    }
    else{
      subData <- reactiveValues(
        counts=Downsample(counts(vals$counts), newcounts=floor(2^seq.int(from=log2(input$minSim), to=log2(input$maxSim), length.out=10)), iterations=input$iterations)
      )
      output$downDone <- renderPlot({
        heatmap(as.matrix(subData$counts[order(apply(subData$counts[,,10,1],1,sum),decreasing=TRUE)[1:20],,10,1]))
      })
    }
  })

  runDiffPower <- observeEvent(input$runDifferentialPower, {

    if(is.null(vals$counts)){
      alert("Warning: Upload data first!")
    }
    else{
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
    }
  })
})
