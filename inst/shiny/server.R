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
    original = getShinyOption("inputSCEset")
  )
  
  #Update all of the columns that depend on pvals columns
  updateAllPdataInputs <- function(){
    updateSelectInput(session, "colorDims",
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
  }
  
  #-----------------------------------------------------------------------------
  # Page 1: Upload
  #-----------------------------------------------------------------------------
  
  #Upload data through shiny app
  observeEvent(input$uploadData, {
    withBusyIndicatorServer("uploadData", {
      vals$counts <- createSCESet(countfile = input$countsfile$datapath,
                                  annotfile = input$annotfile$datapath,
                                  featurefile = input$featurefile$datapath)
      updateAllPdataInputs()
      updateSelectInput(session, "deletesamplelist",
                        choices = rownames(pData(vals$counts)))
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
  
  #Render histogram of read counts per cell
  output$countshist <- renderPlotly({
    if(!(is.null(vals$counts))){
      f <- list(family = "Courier New, monospace", size = 18, color = "#7f7f7f")
      x <- list(title = "Reads per cell", titlefont = f)
      y <- list(title = "Number of cells", titlefont = f)
      plot_ly(x = apply(scater::counts(vals$counts), 2, function(x) sum(x)), type="histogram") %>%
        layout(xaxis=x, yaxis=y)
    } else {
      plotly_empty(type="scatter") %>% add_trace(mode="lines")
      }
  })
  
  #Render histogram of genes detected per cell
  output$geneshist <- renderPlotly({
    if(!(is.null(vals$counts))){
      f <- list(family = "Courier New, monospace", size = 18, color = "#7f7f7f")
      x <- list(title = "Genes detected per cell", titlefont = f)
      y <- list(title = "Number of cells", titlefont = f)
      plot_ly(x = apply(scater::counts(vals$counts), 2, function(x) sum(x>0)), type="histogram") %>%
        layout(xaxis=x, yaxis=y)
    } else {
      plotly_empty(type="scatter") %>% add_trace(mode="lines")
      }
  })  
  
  #Render summary table
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

  #Reset the data to the original uploaded dataset
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
 
  #Create UI for filtering samples based on annotations
  output$FilterSamples <- renderUI({
    annot = ''
    if (!is.null(vals$counts)){
      annot = colnames(pData(vals$counts))
    }
    selectInput("filteredSample", "Select Annotation:", c("none", annot))
  })
  observeEvent(input$filteredSample, {
    output$filterSampleOptions <- renderUI({
      if (input$filteredSample != "none")({
        L = vector("list", 2)
        L[[1]] = list(checkboxGroupInput("filterSampleChoices", 
                                         label = "Select samples to keep",
                                         choices = unique(pData(vals$counts)[,input$filteredSample])))
        L[[2]] = list(actionButton("runFilterSample", "Filter"))
        return(L)
      }) else {
        L = list()
      }
    })
  })
  
  #Filter the selected samples
  observeEvent(input$runFilterSample, {
    filter = pData(vals$counts)[,input$filteredSample] %in% input$filterSampleChoices
    vals$counts = vals$counts[,filter]
  })
   
  #Create UI for filtering genes based on feature annotations
  observeEvent(input$filteredFeature, {
    output$filterFeatureOptions <- renderUI({
      if (input$filteredFeature != "none")({
        L = vector("list", 2)
        L[[1]] = list(checkboxGroupInput("filterFeatureChoices", 
                                         label = "Select features to keep",
                                         choices = unique(fData(vals$counts)[,input$filteredFeature])))
        L[[2]] = list(actionButton("runFilterFeature", "Filter"))
        return(L)
      }) else {
        L = list()
      }
    })
  })
  
  #Filter the selected features
  observeEvent(input$runFilterFeature, {
    filter = fData(vals$counts)[,input$filteredFeature] %in% input$filterFeatureChoices
    vals$counts = vals$counts[filter,]
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
  
  #Plot dimensionality reduction plot
  drDataframe <- observeEvent(input$plotData, {
    withBusyIndicatorServer("plotData", {
      if(is.null(vals$counts)){
        alert("Warning: Upload data first!")
      }
      #Sebastian's Code
      else{
        if(input$selectDimRed == "PCA"){
          vals$PCA <- getPCA(vals$counts)
          #updateVals('PCA')
          algo <- vals$PCA
          PC <- paste("PC",1:nrow(pData(vals$counts)),sep="")
          Variances <- attr(vals$PCA,"percentVar")*100
          var <- data.frame(PC,Variances)
          output$pctable <- renderTable(var[1:10,])
        }
        else{
          vals$TSNE <- getTSNE(vals$counts)
          algo <- vals$TSNE
          output$pctable <- renderTable(data.frame(NULL))
        }
        g <- plotDimRed(input$selectDimRed, algo, vals$counts, input$colorDims, input$pcX, input$pcY)
        output$dimredPlot <- renderPlotly({
          ggplotly(g)
        })
      }
    })
    #End Sebastian's Code
  })

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
  
  ### Start Emma's Code (Still being modulerized)
  # Run PCA or tSNE if selected as data set for clustering
  # This will eventually be replaced with a check if the PCA/tSNE values already exist (Sebastian)
  observeEvent(input$selectDataC, {
    if(input$selectDataC == "PCA Components") {
      pc1 <- input$pcX_Clustering_Data
      pc2 <- input$pcY_Clustering_Data
      pca <- scater::plotPCA(vals$counts, return_SCESet=TRUE)
      pca <- data.frame(reducedDimension(pca))
      pca <- setNames(cbind(rownames(pca), pca, row.names=NULL), c("Sample", colnames(pca)))
      pData(vals$counts)$PCA <- pca
    } else if(input$selectDataC == "tSNE Components") {
      tsne <- scater::plotTSNE(vals$counts, return_SCESet=TRUE)
      tsne <- data.frame(reducedDimension(tsne))
      tsne <- setNames(cbind(rownames(tsne), tsne, row.names=NULL), c("Sample", colnames(tsne)))
      pData(vals$counts)$tSNE <- tsne
    }
  })
  
  # If K-Means is Selected
  # This assumes the PCAs are stored in pData(vals$counts)$PCA (this will change based off
  # of Sebastian's implementation)
  observeEvent(input$numberKClusters, {
    k <- input$numberKClusters
    if(input$selectDataC == "Raw Data") {
      cl <- kmeans(t(exprs(vals$counts)), k) # Need to Transpose vals$counts so samples are clustered
      pData(vals$counts)$Kmeans <- cl$cluster
      updateAllPdataInputs()
    } else if(input$selectDataC == "PCA Components") {
      pc1 <- input$pcX_Clustering_Data
      pc2 <- input$pcY_Clustering_Data
      cl <- kmeans(pData(vals$counts)$PCA[,(strtoi(strsplit(pc1, split = "PC")[[1]][2])+1):(strtoi(strsplit(pc2, split = "PC")[[1]][2])+1)], k)
      pData(vals$counts)$Kmeans <- cl$cluster
      updateAllPdataInputs()
    } else if(input$selectDataC == "tSNE Components") {
      cl <- kmeans(tsne[,2:3], k)
      pData(vals$counts)$Kmeans <- cl$cluster
      updateAllPdataInputs()
    }
    
  })
  
  # If visualize with PCA is Selected
  # Run PCA or tSNE if selected as data set for clustering
  # This will eventually be replaced with a check if the PCA/tSNE values already exist (Sebastian)
  clusterDataFramePCA <- observeEvent(input$plotClustersPCA, {
    ### Run PCA
    pca <- scater::plotPCA(vals$counts, return_SCESet=TRUE)
    pca <- data.frame(reducedDimension(pca))
    pca <- setNames(cbind(rownames(pca), pca, row.names=NULL), c("Sample", colnames(pca)))
    
    ### Plotting Parameters
    PC1 <- input$pcX_Clustering_Plot
    PC2 <- input$pcY_Clustering_Plot
    # Color
    if (input$colorClusters_Plot == "Cluster Label"){
      pca$color <- factor(pData(vals$counts)$Kmeans)
    } else {
      pca$color <- eval(parse(text = paste("pData(vals$counts)$",input$colorClusters_Plot,sep="")))
    }
    # Shape
    if (input$shapeClusters_Plot == "Cluster Label"){
      pca$shape <- factor(pData(vals$counts)$Kmeans)
    } else {
      pca$shape <- eval(parse(text = paste("pData(vals$counts)$",input$shapeClusters_Plot,sep="")))
    }
    # Label
    pca$label <- rownames(pData(vals$counts))
    
    ### Create Plot Object
    g <- ggplot(pca, aes(PC1, PC2, label=label, color=color, shape=shape)) + 
      geom_point() +
      theme(legend.title=element_blank())
    
    output$clusterPlot <- renderPlotly({
      ggplotly(g)
    })
  }) 
  
  # If plot dendogram is selected
  # Sebastian: Change input to PCA/tSNE
  clusterDataFrameDendogram <- observeEvent(input$plotClustersDendogram, {
    if(input$selectDataC == "Raw Data") {
      e <- exprs(vals$counts)
    } else if(input$selectDataC == "PCA Components") {
      pc1 <- input$pcX_Clustering_Data
      pc2 <- input$pcY_Clustering_Data
      # Get PCA Values (in format such that rows are PCs, columns are samples)
      e <- data.frame(t(pData(vals$counts)$PCA[,(strtoi(strsplit(pc1, split = "PC")[[1]][2])+1):(strtoi(strsplit(pc2, split = "PC")[[1]][2])+1)]))
      colnames(e) <- colnames(exprs(vals$counts))
      e <- e[-1,]
    } else if(input$selectDataC == "tSNE Components") {
      e <- data.frame(t(pData(vals$counts)$tSNE))
      colnames(e) <- colnames(exprs(vals$counts))
      e <- e[-1,]
    }
    d <- dist(t(e))
    h <- hclust(d, "ward.D")
    k <- as.integer(input$numberHClusters)
    pData(vals$counts)$Hierarchical <- cutree(h, k=k)
    updateAllPdataInputs()
    
    output$dendoPlot <- renderPlot({
      plot(h, hang=-1, main=sprintf("%s Clusters", k))
      if (k > 1) {
        rect.hclust(h, k=k, border="red")
      }
    })
  })
  
  # If plot hiererchical clustering is selected
  # Sebastian: Change input to PCA/tSNE
  clusterDataFramePhylogenetic <- observeEvent(input$plotClustersPTree, {
    if(input$selectDataC == "Raw Data") {
      e <- exprs(vals$counts)
    } else if(input$selectDataC == "PCA Components") {
      pc1 <- input$pcX_Clustering_Data
      pc2 <- input$pcY_Clustering_Data
      # Get PCA Values (in format such that rows are PCs, columns are samples)
      e <- data.frame(t(pData(vals$counts)$PCA[,(strtoi(strsplit(pc1, split = "PC")[[1]][2])+1):(strtoi(strsplit(pc2, split = "PC")[[1]][2])+1)]))
      colnames(e) <- colnames(exprs(vals$counts))
      e <- e[-1,]
    } else if(input$selectDataC == "tSNE Components") {
      e <- data.frame(t(pData(vals$counts)$tSNE))
      colnames(e) <- colnames(exprs(vals$counts))
      e <- e[-1,]
    }
    d <- dist(t(e))
    h <- hclust(d, "ward.D")
    
    output$phyloPlot <- renderPlot({
      plot(as.phylo(h), type="unrooted")
    })
  })
  
  # clusterDataFrametSNE <- observeEvent(input$plotClustersTSNE, {
  #   ### Run tSNE
  #   tsne <- scater::plotTSNE(vals$counts, return_SCESet=TRUE)
  #   tsne <- data.frame(reducedDimension(tsne))
  #   tsne <- setNames(cbind(rownames(tsne), tsne, row.names=NULL), c("Sample", colnames(tsne)))
  #   
  #   ### Plotting Parameters
  #   # Color
  #   if (input$colorClusters == "Cluster Label"){
  #     tsne$color <- factor(phenoData(vals$counts)$Kmeans$cluster)
  #   } else {
  #     tsne$color <- eval(parse(text = paste("pData(vals$counts)$",input$colorClusters,sep="")))
  #   }
  #   # Shape
  #   if (input$shapeClusters == "Cluster Label"){
  #     tsne$shape <- factor(phenoData(vals$counts)$Kmeans)
  #   } else {
  #     tsne$shape <- eval(parse(text = paste("pData(vals$counts)$",input$shapeClusters,sep="")))
  #   }
  #   # Label
  #   tsne$label <- rownames(pData(vals$counts))
  #   
  #   ### Create Plot Object
  #   g <- ggplot(tsne, aes(X1, X2, label=label, color=color, shape=shape)) + 
  #     geom_point() +
  #     theme(legend.title=element_blank())
  #   
  #   output$clusterPlot <- renderPlotly({
  #     ggplotly(g)
  #   })
  # })
  
  #clusterDataFrameScatter <- observeEvent( {
    
  #})
  
  #clusterDataFrameTable <- observeEvent(input$makeCTable, {
   # k <- input$numberClusters
    #cl <- kmeans(t(exprs(vals$counts)), k)
    #output$clusterTable <- renderDataTable({
     # data.frame(cl$cluster)})
#  })
      
  ### END Emma's Note
  


  #-----------------------------------------------------------------------------
  # Page 4: Differential Expression
  #-----------------------------------------------------------------------------
  
  #For conditions with more than two factors, select the factor of interest
  output$selectDiffex_conditionofinterestUI <- renderUI({
    if(!is.null(vals$counts)){
      if(length(unique(pData(vals$counts)[,input$selectDiffex_condition])) > 2){
        selectInput("selectDiffex_conditionofinterest",
                    "Select Factor of Interest",
                    unique(sort(pData(vals$counts)[,input$selectDiffex_condition])))
      }
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
  
  
  output$colorBarCondition <- renderUI({
    if (is.null(vals$counts)){
      selectInput("colorBar_Condition", "Select Condition", c())
    } else {
      selectizeInput("colorBar_Condition", "Select Condition", 
                     names(pData(vals$counts)), multiple=TRUE)
    }
  })
  
  annotationColors <- reactiveValues(cols = list())
  
  output$HeatmapSampleAnnotations <- renderUI({
    if (!is.null(vals$counts)){
      h = input$colorBar_Condition
      L = lapply(1:length(h), function(i) colourGroupInput(paste0("colourGroup", i)))
      annotationColors$cols = lapply(1:length(h), function(i) callModule(colourGroup, paste0("colourGroup", i),
                                                                  heading = h[i], 
                                                                  options = unique(unlist(pData(vals$counts)[,h[i]]))))
      return(L)
    }
    
  })

  #Plot the differential expression results
  output$diffPlot <- renderPlot({
    if(!is.null(vals$diffexgenelist)){
      if (input$displayHeatmapColorBar){
        if (is.null(input$colorBar_Condition)){
          colors = NULL
        } else {
          colors = lapply(annotationColors$cols, function(col) col())
          names(colors) = input$colorBar_Condition
          if (is.null(colors[[length(colors)]][[1]])){colors=NULL}
        }
      } else {
        colors=NULL
      }
      
      draw(plot_DiffEx(vals$counts, input$colorBar_Condition,
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
      updateSelectInput(session, "filteredFeature", choices = c("none", colnames(fData(vals$counts))))
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
