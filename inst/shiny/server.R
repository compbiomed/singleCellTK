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
library(GGally)
library(data.table)
library(MAST)
library(rsvd)
library(pcaMethods)
library(colourpicker)
library(gridExtra)
library(cluster)
library(ggtree)


#1GB max upload size
options(shiny.maxRequestSize=1000*1024^2)

# Define server logic required to draw a histogram
shinyServer(function(input, output, session) {
  
  #-----------------------------------------------------------------------------
  # MISC - Used throughout app
  #-----------------------------------------------------------------------------
  
  #reactive values object
  vals <- reactiveValues(
    counts = getShinyOption("inputSCEset"),
    original = getShinyOption("inputSCEset"),
    PCA = NULL
  )
  
  #Update all of the columns that depend on pvals columns
  updateAllPdataInputs <- function(){
    updateSelectInput(session, "colorDims",
                      choices = colnames(pData(vals$counts)))
    updateSelectInput(session, "clusterChoice",
                      choices = colnames(pData(vals$counts)))
    updateSelectInput(session, "colorClusters",
                      choices = colnames(pData(vals$counts)))
    updateSelectInput(session, "shapeClusters",
                      choices = colnames(pData(vals$counts)))
    updateSelectInput(session, "colorClusters_Plot",
                      choices = c("Cluster Label", colnames(pData(vals$counts))))
    updateSelectInput(session, "shapeClusters_Plot",
                      choices = c("Cluster Label", colnames(pData(vals$counts))))
    updateSelectInput(session, "colorClusters_MAST",
                      choices = colnames(pData(vals$counts)))
    updateSelectInput(session, "selectDiffex_condition",
                      choices = colnames(pData(vals$counts)))
    updateSelectInput(session, "subCovariate",
                      choices = colnames(pData(vals$counts)))
    updateSelectInput(session, "selectAdditionalVariables",
                      choices = colnames(pData(vals$counts)))
    updateSelectInput(session, "deletepdatacolumn",
                      choices = colnames(pData(vals$counts)))
    updateSelectInput(session, "colorBy",
                      choices = c("No Color", "Gene Expression", colnames(pData(vals$counts))))
    updateSelectInput(session, "shapeBy",
                      choices = c("No Shape", colnames(pData(vals$counts))))
    updateSelectInput(session, "Knumber",
                      choices = 1:nrow(pData(vals$counts)))
  }
  
  # Close app on quit
  session$onSessionEnded(stopApp)
  
  #-----------------------------------------------------------------------------
  # Page 1: Upload
  #-----------------------------------------------------------------------------
  
  #Upload data through shiny app
  observeEvent(input$uploadData, {
    withBusyIndicatorServer("uploadData", {
      if(input$uploadChoice == "files"){
        vals$counts <- createSCESet(countfile = input$countsfile$datapath,
                                    annotfile = input$annotfile$datapath,
                                    featurefile = input$featurefile$datapath)
      } else {
        vals$counts <- createSCESet(countfile = eval(as.symbol(input$selectExampleData))$counts,
                                    annotfile = eval(as.symbol(input$selectExampleData))$annot,
                                    inputdataframes = TRUE)
      }
      updateAllPdataInputs()
      updateSelectInput(session, "deletesamplelist",
                        choices = rownames(pData(vals$counts)))
      updateSelectInput(session, "colorGenes",
                        choices = rownames(vals$counts))
      updateSelectInput(session, "pcX",
                        choices = paste("PC",1:nrow(pData(vals$counts)),sep=""),
                        selected = "PC1")
      updateSelectInput(session, "pcY",
                        choices = paste("PC",1:nrow(pData(vals$counts)),sep=""),
                        selected = "PC2")
      updateSelectInput(session, "pcX_Clustering_Data",
                        choices = paste("PC",1:nrow(pData(vals$counts)),sep=""),
                        selected = "PC1")
      updateSelectInput(session, "pcY_Clustering_Data",
                        choices = paste("PC",1:nrow(pData(vals$counts)),sep=""),
                        selected = "PC2")
      updateSelectInput(session, "pcX_Clustering_Plot",
                        choices = paste("PC",1:nrow(pData(vals$counts)),sep=""),
                        selected = "PC1")
      updateSelectInput(session, "pcY_Clustering_Plot",
                        choices = paste("PC",1:nrow(pData(vals$counts)),sep=""),
                        selected = "PC2")
      updateSelectInput(session, "numberKClusters",
                        choices = 1:nrow(pData(vals$counts)))
      updateSelectInput(session, "numberHClusters",
                        choices = 1:nrow(pData(vals$counts)))
      updateSelectInput(session, "Knumber",
                        choices = 1:nrow(pData(vals$counts)))
      insertUI(
        selector = '#uploadAlert',
        ## wrap element in a div with id for ease of removal
        ui = tags$div(class="alert alert-success alert-dismissible", HTML("<span \
                      class='glyphicon glyphicon-ok' aria-hidden='true'></span> \
                      Successfully Uploaded! <button type='button' class='close' \
                      data-dismiss='alert'>&times;</button>"))
        )
      vals$original <- vals$counts
      vals$PCA <- NULL
      vals$TSNE <- NULL
    })
  })
  
  #-----------------------------------------------------------------------------
  # Page 2: Data Summary
  #-----------------------------------------------------------------------------

  #Render data table if there are fewer than 50 samples
  output$contents <- renderDataTable({
    if(!(is.null(vals$counts)) && nrow(pData(vals$counts)) < 50){
      temptable <- cbind(rownames(fData(vals$counts)),exprs(vals$counts))
      colnames(temptable)[1] <- "Gene"
      temptable
    }
  }, options = list(scrollX = TRUE))

  #Render summary table
  output$summarycontents <- renderTable({
    if(!(is.null(vals$counts))){
      summarizeTable(vals$counts)
    }
  })

  #Filter the data based on the options
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
      #Refresh things for the clustering tab
      vals$PCA <- NULL
      vals$TSNE <- NULL
      updateAllPdataInputs()
    }
  })

  #Reset the data to the original uploaded dataset
  observeEvent(input$resetData, {
    if(is.null(vals$original)){
      alert("Warning: Upload data first!")
    }
    else{
      vals$counts <- vals$original
      updateSelectInput(session, "deletesamplelist",
                        choices = rownames(pData(vals$counts)))
      #Refresh things for the clustering tab
      vals$PCA <- NULL
      vals$TSNE <- NULL
      updateAllPdataInputs()
    }
  })
  
  #Delete a column from the pData annotations
  observeEvent(input$deletepDatabutton, {
    if(is.null(vals$original)){
      alert("Warning: Upload data first!")
    }
    else{
      pData(vals$counts) <- pData(vals$counts)[,!(colnames(pData(vals$counts)) %in% input$deletepdatacolumn), drop=F]
      updateAllPdataInputs()
    }
  })
  
  output$downloadSCESet <- downloadHandler(
    filename = function() {
      paste('SCESet-', Sys.Date(), '.rds', sep='')
    },
    content = function(file) {
      saveRDS(vals$counts, file)
  })
  
  #-----------------------------------------------------------------------------
  # Page 3: DR & Clustering
  #-----------------------------------------------------------------------------

  #Multi PCA plot by Lloyd
  #singlCellTK::runDimRed may want to change
  multipcaDataFrame <- observeEvent(input$plotPCA, {
   withBusyIndicatorServer("plotPCA", {
     if(is.null(vals$counts)){
       alert("Warning: Upload data first!")
     }
     else{
       g = runPCA(plot.type = input$plotTypeId,
                  method = input$pcaAlgorithm,
                  countm = exprs(vals$counts),
                  annotm = pData(vals$counts),
                  featurem = fData(vals$counts),
                  involving.variables = input$pcaCheckbox,
                  additional.variables = input$selectAdditionalVariables,
                  colorClusters = input$colorClusters_MAST)
       output$pcaPlot <- renderPlot({
         g
       })
     }
   }) 
  })
  #end Lloyd's Code
  
  output$clusterPlot <- renderPlotly({
    if(is.null(vals$counts)){
      alert("Warning: Upload data first!")
    } else{
      if(input$dimRedPlotMethod=="PCA"){
        if (is.null(vals$PCA)) {
          vals$PCA <- getPCA(vals$counts)
        }
        if(!is.null(vals$PCA)){
          if (input$colorBy != "Gene Expression") {
            g <- singleCellTK::plotPCA(vals$counts, vals$PCA, input$colorBy, input$shapeBy, input$pcX, input$pcY)
          } else if (input$colorGenes == ""){
            g <- singleCellTK::plotPCA(vals$counts, vals$PCA, "No Color", "No Shape", input$pcX, input$pcY)
          } 
          ggplotly(g)
        } else {
          ggplotly(ggplot() + geom_point())
        }
      } else if(input$dimRedPlotMethod=="tSNE"){
        if (is.null(vals$TSNE)) {
          vals$TSNE <- getTSNE(vals$counts)
        }
        if(!is.null(vals$TSNE)){
          if (input$colorBy != "Gene Expression") {
            g <- singleCellTK::plotTSNE(vals$counts, vals$TSNE, input$colorBy, input$shapeBy)
          } else if (input$colorGenes == ""){
            g <- singleCellTK::plotTSNE(vals$counts, vals$TSNE, "No Color", "No Shape")
          }
          ggplotly(g)
        } else {
        ggplotly(ggplot() + geom_point())
        }
      } else{
        ggplotly(ggplot() + theme_bw() + theme(plot.background = element_rect(fill='white')) + theme(panel.border = element_rect(colour = "white")))
      }
    }
  })
  
  
  output$pctable <- renderTable(
    if(is.null(vals$counts)){
      alert("Warning: Upload data first!")
    } else{
      if (input$dimRedPlotMethod == "PCA") {
        data.frame(PC = paste("PC",1:nrow(pData(vals$counts)),sep=""), Variances = attr(vals$PCA,"percentVar")*100)[1:10,]
      }
    })
  
  output$geneExpressionPlot <- renderPlot(
    if(is.null(vals$counts)){
      alert("Warning: Upload data first!")
    } else {
      if(input$dimRedPlotMethod=="PCA"){
        if (is.null(vals$PCA)) {
          vals$PCA <- getPCA(vals$counts)
        }
        if(!is.null(vals$PCA)){
          if(is.null(input$colorBy)) {
            return()}
          if (input$colorBy == "Gene Expression") {
            if (is.null(input$colorGenes)){
              ggplot() + theme_bw() + theme(plot.background = element_rect(fill='white')) + theme(panel.border = element_rect(colour = "white"))
            } else {
              g <- plotBiomarker(vals$counts, input$colorGenes, input$colorBinary, "PCA", input$shapeBy,vals$PCA, input$pcX, input$pcY)
              g
            }
          }
        }
      } else if(input$dimRedPlotMethod == "tSNE"){
        if (is.null(vals$TSNE)) {
          vals$TSNE <- getTSNE(vals$counts)
        } 
        if(!is.null(vals$TSNE)){
          if (input$colorBy == "Gene Expression") {
            if (is.null(input$colorGenes)){
              ggplot() + theme_bw() + theme(plot.background = element_rect(fill='white')) + theme(panel.border = element_rect(colour = "white"))
            } else {
              g <- plotBiomarker(vals$counts, input$colorGenes, input$colorBinary, "tSNE", input$shapeBy, vals$TSNE)
              g
            }
          }
        }
      }
    }
  )
  
  output$treePlot <- renderPlot(
    if(is.null(vals$counts)){
      alert("Warning: Upload data first!")
    } else {
      if(input$dimRedPlotMethod=="Dendrogram"){
        if (input$clusteringAlgorithmD == "Phylogenetic Tree") {
          data <- getClusterInputData(vals$counts, input$selectClusterInputData, vals)
          d <- dist(data)
          h <- hclust(d, "ward.D")
          g <- ggtree(as.phylo(h), layout = "circular", open.angle = 360) + geom_tiplab2(size=2)
          g
        } else if (input$clusteringAlgorithmD == "Hierarchical") {
          data <- getClusterInputData(vals$counts, input$selectClusterInputData, vals)
          d <- dist(data)
          h <- hclust(d, "ward.D")
          g <- ggtree(as.phylo(h)) + theme_tree2() + geom_tiplab(size=2)
          g
        }
      } 
    }
  )
  
  clusterDataFrame <- observeEvent(input$clusterData, {
    withBusyIndicatorServer("clusterData", {
      if(input$clusteringAlgorithm == "K-Means"){
        data <- getClusterInputData(vals$counts, input$selectClusterInputData, vals)
        koutput <- kmeans(data, input$Knumber)
        pData(vals$counts)$Kmeans <- factor(koutput$cluster)
        updateAllPdataInputs()
      } else if(input$clusteringAlgorithm == "Clara") {
        data <- getClusterInputData(vals$counts, input$selectClusterInputData, vals)
        coutput <- clara(data, input$Cnumber)
        pData(vals$counts)$Clara <- factor(coutput$clustering)
        updateAllPdataInputs()
      } 
    })
  })
  
  #-----------------------------------------------------------------------------
  # Page 4: Differential Expression
  #-----------------------------------------------------------------------------
  
  #For conditions with more than two factors, select the factor of interest
  output$selectDiffex_conditionofinterestUI <- renderUI({
    if(length(unique(pData(vals$counts)[,input$selectDiffex_condition])) > 2){
      selectInput("selectDiffex_conditionofinterest",
                  "Select Factor of Interest",
                  unique(sort(pData(vals$counts)[,input$selectDiffex_condition])))
    }
  })
  
  #Run differential expression
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
                                        clusterCol=input$clusterColumns,
                                        levelofinterest = input$selectDiffex_conditionofinterest)
      })
    }
  })
  
  #Plot the differential expression results
  output$diffPlot <- renderPlot({
    if(!is.null(vals$diffexgenelist)){
      if (input$displayHeatmapColorBar){
        names = names(input)
        conditions = unique(pData(vals$counts)[,input$selectDiffex_condition])
        condNames = paste0("hmColorBar_", conditions)
        names = names[names %in% condNames]
        if (length(names)==0){
          colors = sapply(1:length(conditions), function (i) palette()[(i %% length(palette()))+1])
          names(colors) = conditions
        } else {
          colors = rep("", length(names))
          names(colors) = names
          for (i in 1:length(colors)){
            n = names[i]
            colors[i]<-input[[as.character(n)]]
            names(colors)[i] = gsub("hmColorBar_", "", n)
          }
        }
      } else {
        colors=NULL
      }
      draw(plot_DiffEx(vals$counts, input$selectDiffex_condition,
                  rownames(vals$diffexgenelist), clusterRow=input$clusterRows,
                  clusterCol=input$clusterColumns,
                  displayRowLabels=input$displayHeatmapRowLabels,
                  displayColumnLabels=input$displayHeatmapColumnLabels,
                  displayRowDendrograms=input$displayHeatmapRowDendrograms,
                  displayColumnDendrograms=input$displayHeatmapColumnDendrograms,
                  annotationColors=colors,
                  columnTitle=input$heatmapColumnsTitle))
    }
  }, height=600)

  #Create the color bar options
  output$colorBarOptions <- renderUI({
    if (!is.null(vals$counts)){
      conditions = unique(pData(vals$counts)[,input$selectDiffex_condition])
      L = vector("list", length(conditions))
      for (i in 1:length(L)){
        id=paste0("hmColorBar_", conditions[i])
        if (is.null(input[[id]])){
          color = palette()[(i %% length(palette()))+1]
        } else {
          color = input[[id]]
        }
        L[[i]] = list(colourInput(
          id, 
          conditions[i], 
          color))
      }
      return(L)  
    }
  })
  
  #Render the interactive heatmap
  output$interactivediffPlot <- renderD3heatmap({
    if(!is.null(vals$diffexgenelist)){
      plot_d3DiffEx(vals$counts, input$selectDiffex_condition,
                    rownames(vals$diffexgenelist), clusterRow=input$clusterRows,
                    clusterCol=input$clusterColumns)
    }
  })
  
  #Create the differential expression results table
  output$diffextable <- renderDataTable({
    if(!is.null(vals$diffexgenelist)){
      temptable <- cbind(rownames(vals$diffexgenelist),data.frame(vals$diffexgenelist))
      colnames(temptable)[1] <- "Gene"
      temptable
    }
  })

  # Download the differential expression results table
  output$downloadGeneList <- downloadHandler(
    filename = function() {
      paste("diffex_results-", Sys.Date(), ".csv", sep="")
    },
    content = function(file) {
      write.csv(vals$diffexgenelist, file)
    }
  )
  
  observeEvent(input$saveBiomarker, {
    if(is.null(input$biomarkerName)){
      alert("Warning: Specify biomarker name!")
    } else{
      fData(vals$counts)[,input$biomarkerName] <- ifelse(rownames(fData(vals$counts)) %in% rownames(vals$diffexgenelist), 1,0)
    }
  })
  
  #-----------------------------------------------------------------------------
  # Page 5: Subsampling
  #-----------------------------------------------------------------------------

  #Run subsampling analysis
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

  #Run differential power analysis
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
