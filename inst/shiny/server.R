#1GB max upload size
options(shiny.maxRequestSize = 1000 * 1024 ^ 2)

# Define server logic required to draw a histogram
shinyServer(function(input, output, session) {

  #-----------------------------------------------------------------------------
  # MISC - Used throughout app
  #-----------------------------------------------------------------------------

  #reactive values object
  vals <- reactiveValues(
    counts = getShinyOption("inputSCEset"),
    original = getShinyOption("inputSCEset"),
    selectgenes = rownames(getShinyOption("inputSCEset"))
  )

  #Update all of the columns that depend on pvals columns
  updateAllPdataInputs <- function(){
    pdata_options <- colnames(colData(vals$counts))
    updateSelectInput(session, "selectDiffex_condition",
                      choices = pdata_options)
    updateSelectInput(session, "subCovariate",
                      choices = pdata_options)
    updateSelectInput(session, "deleterowdatacolumn",
                      choices = pdata_options)
    updateSelectInput(session, "hurdlecondition",
                      choices = pdata_options)
    updateSelectInput(session, "colorBy",
                      choices = c("No Color", "Gene Expression", pdata_options))
    updateSelectInput(session, "shapeBy",
                      choices = c("No Shape", pdata_options))
    updateSelectInput(session, "Knumber",
                      choices = 1:ncol(vals$counts))
    selectthegenes <- rownames(vals$counts)
    updateSelectizeInput(session, "colorGenes",
                      choices = selectthegenes, server = TRUE)
  }

  # Close app on quit
  session$onSessionEnded(stopApp)

  #-----------------------------------------------------------------------------
  # Page 1: Upload
  #-----------------------------------------------------------------------------

  #Upload data through shiny app
  observeEvent(input$uploadData, {
    withBusyIndicatorServer("uploadData", {
      if (input$uploadChoice == "files"){
        vals$original <- createSCE(countfile = input$countsfile$datapath,
                                   annotfile = input$annotfile$datapath,
                                   featurefile = input$featurefile$datapath,
                                   create_logcounts = input$createLogcounts)
      } else {
        data(list = paste0(input$selectExampleData, "_sce"))
        vals$original <- base::eval(parse(text = paste0(input$selectExampleData, "_sce")))
      }
      vals$counts <- vals$original
      updateAllPdataInputs()
      updateSelectInput(session, "deletesamplelist",
                        choices = colnames(vals$counts))
      updateSelectInput(session, "dimRedAssaySelect",
                        choices = names(assays(vals$counts)))
      updateSelectizeInput(session, "colorGenes",
                        choices = rownames(vals$counts))
      updateSelectInput(session, "pcX",
                        choices = paste("PC", 1:ncol(vals$counts), sep = ""),
                        selected = "PC1")
      updateSelectInput(session, "pcY",
                        choices = paste("PC", 1:ncol(vals$counts), sep = ""),
                        selected = "PC2")
      updateSelectInput(session, "Knumber",
                        choices = 1:ncol(vals$counts))
      updateSelectInput(session, "Cnumber",
                        choices = 1:ncol(vals$counts))
      insertUI(
        selector = "#uploadAlert",
        ## wrap element in a div with id for ease of removal
        ui = tags$div(class = "alert alert-success alert-dismissible", HTML("<span \
                      class='glyphicon glyphicon-ok' aria-hidden='true'></span> \
                      Successfully Uploaded! <button type='button' class='close' \
                      data-dismiss='alert'>&times;</button>"))
      )
    })
  })

  #-----------------------------------------------------------------------------
  # Page 2: Data Summary
  #-----------------------------------------------------------------------------

  #Render data table if there are fewer than 50 samples
  output$contents <- renderDataTable({
    if (!(is.null(vals$counts)) && ncol(vals$counts) < 50){
      temptable <- cbind(rownames(vals$counts), log2(assay(vals$counts,
                                                           "counts") + 1))
      colnames(temptable)[1] <- "Gene"
      temptable
    }
  }, options = list(scrollX = TRUE))

  #Render histogram of read counts per cell
  output$countshist <- renderPlotly({
    if (!(is.null(vals$counts))){
      f <- list(family = "Courier New, monospace", size = 18, color = "#7f7f7f")
      x <- list(title = "Reads per cell", titlefont = f)
      y <- list(title = "Number of cells", titlefont = f)
      plot_ly(x = apply(assay(vals$counts, "counts"), 2, function(x) sum(x)),
              type = "histogram") %>%
        layout(xaxis = x, yaxis = y)
    } else {
      plotly_empty(type = "scatter") %>% add_trace(mode = "lines")
    }
  })

  #Render histogram of genes detected per cell
  output$geneshist <- renderPlotly({
    if (!(is.null(vals$counts))){
      f <- list(family = "Courier New, monospace", size = 18, color = "#7f7f7f")
      x <- list(title = "Genes detected per cell", titlefont = f)
      y <- list(title = "Number of cells", titlefont = f)
      plot_ly(x = apply(assay(vals$counts, "counts"), 2,
                        function(x) sum(x > 0)), type = "histogram") %>%
        layout(xaxis = x, yaxis = y)
    } else {
      plotly_empty(type = "scatter") %>% add_trace(mode = "lines")
      }
  })

  #Render summary table
  output$summarycontents <- renderTable({
    if (!(is.null(vals$counts))){
      summarizeTable(vals$counts)
    }
  })

  #Filter the data based on the options
  observeEvent(input$filterData, {
    if (is.null(vals$original)){
      alert("Warning: Upload data first!")
    }
    else{
      vals$counts <- vals$original
      deletesamples <- input$deletesamplelist
      vals$counts <- filterSCData(vals$counts,
                                  deletesamples = deletesamples,
                                  remove_noexpress = input$removeNoexpress,
                                  remove_bottom = 0.01 * input$LowExpression,
                                  minimum_detect_genes = input$minDetectGenect)
      #Refresh things for the clustering tab
      updateAllPdataInputs()
    }
  })

  #Reset the data to the original uploaded dataset
  observeEvent(input$resetData, {
    if (is.null(vals$original)){
      alert("Warning: Upload data first!")
    }
    else{
      vals$counts <- vals$original
      updateSelectInput(session, "deletesamplelist",
                        choices = colnames(vals$counts))
      #Refresh things for the clustering tab
      updateAllPdataInputs()
    }
  })

  #Delete a column from the colData annotations
  observeEvent(input$deleterowDatabutton, {
    if (is.null(vals$original)){
      alert("Warning: Upload data first!")
    }
    else{
      colData(vals$counts) <- colData(vals$counts)[, !(colnames(colData(vals$counts)) %in% input$deleterowdatacolumn), drop = F]
      updateAllPdataInputs()
    }
  })

  #Create UI for filtering samples based on annotations
  output$FilterSamples <- renderUI({
    annot <- ""
    if (!is.null(vals$counts)){
      annot <- colnames(colData(vals$counts))
    }
    selectInput("filteredSample", "Select Annotation:", c("none", annot))
  })

  observeEvent(input$filteredSample, {
    output$filterSampleOptions <- renderUI({
      if (input$filteredSample != "none")({
        if (length(unique(colData(vals$counts)[, input$filteredSample])) < 100){
          L <- vector("list", 3)
          L[[1]] <- renderText("Select samples to keep")
          L[[2]] <- wellPanel(style = "overflow-y:scroll; max-height: 100px",
                        list(checkboxGroupInput("filterSampleChoices",
                                            label = NULL,
                                            choices = unique(colData(vals$counts)[, input$filteredSample])))
                             )
          L[[3]] <- list(actionButton("runFilterSample", "Filter"))
          return(L)
        } else {
          L <- list(renderText("Annotation must have fewer than 100 options"))
          return(L)
        }
      }) else {
        L <- list()
      }
    })
  })

  #Filter the selected samples
  observeEvent(input$runFilterSample, {
    filter <- colData(vals$counts)[, input$filteredSample] %in% input$filterSampleChoices
    vals$counts <- vals$counts[, filter]
  })

  #Create UI for filtering genes based on feature annotations
  output$filterFeatures <- renderUI({
    fannot <- ""
    if (!is.null(vals$counts)){
      fannot <- colnames(rowData(vals$counts))
    }
    selectInput("filteredFeature", "Select Feature:", c("none", fannot))
  })

  observeEvent(input$filteredFeature, {
    output$filterFeatureOptions <- renderUI({
      if (input$filteredFeature != "none")({
        if (length(unique(rowData(vals$counts)[, input$filteredFeature])) < 100){
          L <- vector("list", 3)
          L[[1]] <- renderText("Select features to keep")
          L[[2]] <- wellPanel(style = "overflow-y:scroll; max-height: 100px",
                             list(checkboxGroupInput("filterFeatureChoices",
                                                     label = NULL,
                                                     choices = unique(rowData(vals$counts)[, input$filteredFeature]))))
          L[[3]] <- list(actionButton("runFilterFeature", "Filter"))
          return(L)
        } else {
          L <- list(renderText("Annotation must have fewer than 100 options"))
          return(L)
        }
      }) else {
        L <- list()
      }
    })
  })

  #Filter the selected features
  observeEvent(input$runFilterFeature, {
    filter <- rowData(vals$counts)[, input$filteredFeature] %in% input$filterFeatureChoices
    vals$counts <- vals$counts[filter, ]
  })

  output$downloadSCE <- downloadHandler(
    filename <- function() {
      paste("SCE-", Sys.Date(), ".rds", sep = "")
    },
    content <- function(file) {
      saveRDS(vals$counts, file)
  })

  #-----------------------------------------------------------------------------
  # Page 3: DR & Clustering
  #-----------------------------------------------------------------------------

  output$clusterPlot <- renderPlotly({
    if (is.null(vals$counts)){
      ggplotly(ggplot())
    } else{
      if (input$dimRedPlotMethod == "PCA"){
        pcadimname <- paste0("PCA", "_", input$dimRedAssaySelect)
        if (is.null(reducedDim(vals$counts, pcadimname))) {
          vals$counts <- getPCA(vals$counts,
                                use_assay = input$dimRedAssaySelect,
                                reducedDimName = pcadimname)
        }
        if (!is.null(reducedDim(vals$counts, pcadimname))){
          if (input$colorBy != "Gene Expression") {
            g <- singleCellTK::plotPCA(vals$counts, input$colorBy,
                                       input$shapeBy, input$pcX, input$pcY,
                                       use_assay = input$dimRedAssaySelect,
                                       reducedDimName = pcadimname)
          } else if (input$colorGenes == ""){
            g <- singleCellTK::plotPCA(vals$counts, "No Color", "No Shape",
                                       input$pcX, input$pcY,
                                       use_assay = input$dimRedAssaySelect,
                                       reducedDimName = pcadimname)
          }
          ggplotly(g)
        } else {
          ggplotly(ggplot() + geom_point())
        }
      } else if (input$dimRedPlotMethod == "tSNE"){
        tsnedimname <- paste0("TSNE", "_", input$dimRedAssaySelect)
        if (is.null(reducedDim(vals$counts, tsnedimname))) {
          vals$counts <- getTSNE(vals$counts,
                                 use_assay = input$dimRedAssaySelect,
                                 reducedDimName = tsnedimname)
        }
        if (!is.null(reducedDim(vals$counts, tsnedimname))){
          if (input$colorBy != "Gene Expression") {
            g <- singleCellTK::plotTSNE(vals$counts, input$colorBy, input$shapeBy,
                                        use_assay = input$dimRedAssaySelect,
                                        reducedDimName = tsnedimname)
          } else if (input$colorGenes == ""){
            g <- singleCellTK::plotTSNE(vals$counts, "No Color", "No Shape",
                                        use_assay = input$dimRedAssaySelect,
                                        reducedDimName = tsnedimname)
          }
          ggplotly(g)
        } else {
          ggplotly(ggplot() + geom_point())
        }
      } else{
        ggplotly(ggplot() + theme_bw() +
                   theme(plot.background = element_rect(fill = "white")) +
                   theme(panel.border = element_rect(colour = "white")))
      }
    }
  })

  #TODO: this doesn't work with multiple pca dims
  output$pctable <- renderTable({
    if (is.null(vals$counts)){
    } else{
      if (input$dimRedPlotMethod == "PCA") {
        if(nrow(pca_variances(vals$counts)) == ncol(vals$counts)){
          data.frame(PC = paste("PC", 1:ncol(vals$counts), sep = ""),
                     Variances = pca_variances(vals$counts)$percentVar * 100)[1:10, ]
        }
      }
    }
  })

  output$geneExpressionPlot <- renderPlot({
    if (is.null(vals$counts)){
    } else {
      if (input$dimRedPlotMethod == "PCA"){
        pcadimname <- paste0("PCA", "_", input$dimRedAssaySelect)
        if (is.null(reducedDim(vals$counts, pcadimname))) {
          vals$counts <- getPCA(vals$counts)
        }
        if (!is.null(reducedDim(vals$counts, pcadimname))){
          if (is.null(input$colorBy)) {
            return()
          }
          if (input$colorBy == "Gene Expression") {
            if (input$colorGeneBy == "Biomarker (from DE tab)"){
              if (input$colorGenesBiomarker == ""){
                ggplot() + theme_bw() + theme(plot.background = element_rect(fill = "white")) + theme(panel.border = element_rect(colour = "white"))
              } else {
                biomarkers <- data.frame(eval(parse(text = paste("rowData(vals$counts)[,'", input$colorGenesBiomarker, "']", sep = ""))))
                rownames(biomarkers) <- rowData(vals$counts)[, "Gene"]
                biomarkers <- rownames(subset(biomarkers, biomarkers[, 1] == 1))
                g <- plotBiomarker(vals$counts, biomarkers, input$colorBinary,
                                   "PCA", input$shapeBy, input$pcX, input$pcY,
                                   use_assay = input$dimRedAssaySelect,
                                   reducedDimName = pcadimname)
                g
              }
            } else if (input$colorGeneBy == "Manual Input") {
              if (is.null(input$colorGenes)){
                ggplot() + theme_bw() + theme(plot.background = element_rect(fill = "white")) + theme(panel.border = element_rect(colour = "white"))
              } else {
                g <- plotBiomarker(vals$counts, input$colorGenes,
                                   input$colorBinary, "PCA", input$shapeBy,
                                   input$pcX, input$pcY,
                                   use_assay = input$dimRedAssaySelect,
                                   reducedDimName = pcadimname)
                g
              }
            }
          }
        }
      } else if (input$dimRedPlotMethod == "tSNE"){
        tsnedimname <- paste0("TSNE", "_", input$dimRedAssaySelect)
        if (is.null(reducedDim(vals$counts, tsnedimname))){
          vals$counts <- getTSNE(vals$counts)
        }
        if (!is.null(reducedDim(vals$counts, tsnedimname))){
          if (input$colorBy == "Gene Expression") {
            if (input$colorGeneBy == "Biomarker (from DE tab)"){
              if (input$colorGenesBiomarker == ""){
                ggplot() + theme_bw() + theme(plot.background = element_rect(fill = "white")) + theme(panel.border = element_rect(colour = "white"))
              } else {
                biomarkers <- data.frame(eval(parse(text = paste("rowData(vals$counts)[,'", input$colorGenesBiomarker, "']", sep = ""))))
                rownames(biomarkers) <- rowData(vals$counts)[, "Gene"]
                biomarkers <- rownames(subset(biomarkers, biomarkers[, 1] == 1))
                g <- plotBiomarker(vals$counts, biomarkers, input$colorBinary,
                                   "tSNE", input$shapeBy,
                                   use_assay = input$dimRedAssaySelect,
                                   reducedDimName = tsnedimname)
                g
              }
            } else if (input$colorGeneBy == "Manual Input") {
              if (is.null(input$colorGenes)){
                ggplot() + theme_bw() + theme(plot.background = element_rect(fill = "white")) + theme(panel.border = element_rect(colour = "white"))
              } else {
                g <- plotBiomarker(vals$counts, input$colorGenes,
                                   input$colorBinary, "tSNE", input$shapeBy,
                                   use_assay = input$dimRedAssaySelect,
                                   reducedDimName = tsnedimname)
                g
              }
            }
          }
        }
      }
    }
  })

  output$treePlot <- renderPlot({
    if (is.null(vals$counts)){
      alert("Warning: Upload data first!")
    } else {
      if (input$dimRedPlotMethod == "Dendrogram"){
        if (input$clusteringAlgorithmD == "Phylogenetic Tree") {
          data <- getClusterInputData(vals$counts, input$selectClusterInputData)
          d <- dist(data)
          h <- hclust(d, input$dendroDistanceMetric)
          g <- ggtree(as.phylo(h), layout = "circular", open.angle = 360) + geom_tiplab2(size = 2)
          g
        } else if (input$clusteringAlgorithmD == "Hierarchical") {
          data <- getClusterInputData(vals$counts, input$selectClusterInputData)
          d <- dist(data)
          h <- hclust(d, input$dendroDistanceMetric)
          g <- ggtree(as.phylo(h)) + theme_tree2() + geom_tiplab(size = 2)
          g
        }
      }
    }
  }, height = 600)

  clusterDataFrame <- observeEvent(input$clusterData, {
    if (is.null(vals$counts)){
      alert("Warning: Upload data first!")
    } else {
      withBusyIndicatorServer("clusterData", {
        currdimname <- NULL
        if (input$selectClusterInputData == "PCA Components") {
          currdimname <- paste0("PCA", "_", input$dimRedAssaySelect)
        } else if (input$selectClusterInputData == "tSNE Components") {
          currdimname <- paste0("TSNE", "_", input$dimRedAssaySelect)
        }
        if (input$clusteringAlgorithm == "K-Means"){
          data <- getClusterInputData(vals$counts, input$selectClusterInputData,
                                      use_assay = input$dimRedAssaySelect,
                                      reducedDimName = currdimname)
          koutput <- kmeans(data, input$Knumber)
          name <- input$clusterName
          colData(vals$counts)[, paste(name)] <- factor(koutput$cluster)
          updateAllPdataInputs()
          updateSelectInput(session, "colorBy",
                            choices = c("No Color", "Gene Expression",
                                        colnames(colData(vals$counts))),
                            selected = input$clusterName)
        } else if (input$clusteringAlgorithm == "Clara") {
          data <- getClusterInputData(vals$counts, input$selectClusterInputData,
                                      use_assay = input$dimRedAssaySelect,
                                      reducedDimName = currdimname)
          coutput <- clara(data, input$Cnumber)
          colData(vals$counts)$Clara <- factor(coutput$clustering)
          updateAllPdataInputs()
          updateSelectInput(session, "colorBy",
                            choices = c("No Color", "Gene Expression",
                                        colnames(colData(vals$counts))),
                            selected = input$clusterName)
        }
      })
    }
  })

  observeEvent(input$reRunTSNE, {
    if (is.null(vals$original)){
      alert("Warning: Upload data first!")
    }
    else{
      withBusyIndicatorServer("reRunTSNE", {
        vals$counts <- getTSNE(vals$counts)
      })
    }
  })

  #-----------------------------------------------------------------------------
  # Page 4: Differential Expression
  #-----------------------------------------------------------------------------

  #For conditions with more than two factors, select the factor of interest
  output$selectDiffex_conditionofinterestUI <- renderUI({
    if (!is.null(vals$counts)){
      if (length(unique(colData(vals$counts)[, input$selectDiffex_condition])) > 2 & input$selectDiffex != "ANOVA"){
        selectInput("selectDiffex_conditionofinterest",
                    "Select Factor of Interest",
                    unique(sort(colData(vals$counts)[, input$selectDiffex_condition])))
      }
    }
  })

  #Run differential expression
  diffexDataframe <- observeEvent(input$runDiffex, {
    if (is.null(vals$counts)){
      alert("Warning: Upload data first!")
    }
    else{
      withBusyIndicatorServer("runDiffex", {
        #run diffex to get gene list and pvalues
        vals$diffexgenelist <- scDiffEx(vals$counts, input$selectDiffex_condition,
                                        input$selectPval, input$selectNGenes, input$applyCutoff,
                                        diffexmethod = input$selectDiffex,
                                        clusterRow = input$clusterRows,
                                        clusterCol = input$clusterColumns,
                                        levelofinterest = input$selectDiffex_conditionofinterest)
        updateSelectInput(session, "colorBar_Condition", selected = input$selectDiffex_condition)
      })
    }
  })

  output$colorBarCondition <- renderUI({
    if (is.null(vals$counts)){
      selectInput("colorBar_Condition", "Select Condition", c())
    } else {
      selectInput("colorBar_Condition", "Select Condition",
                  colnames(colData(vals$counts)), multiple = TRUE)
    }
  })

  annotationColors <- reactiveValues(cols = list())

  output$HeatmapSampleAnnotations <- renderUI({
    if (!is.null(vals$counts) && length(input$colorBar_Condition) > 0){
      h <- input$colorBar_Condition
      L <- lapply(1:length(h), function(i) colourGroupInput(paste0("colorGroup", i)))
      annotationColors$cols <- lapply(1:length(h),
                                      function(i) callModule(colourGroup, paste0("colorGroup", i),
                                                             heading = h[i],
                                                             options = unique(unlist(colData(vals$counts)[, h[i]]))))
      return(L)
    }
  })

  #Plot the differential expression results
  output$diffPlot <- renderPlot({
    if (!is.null(vals$diffexgenelist)){
      if (input$displayHeatmapColorBar){
        if (is.null(input$colorBar_Condition)){
          colors <- NULL
        } else {
          colors <- lapply(annotationColors$cols, function(col) col())
          names(colors) <- input$colorBar_Condition
          if (is.null(colors[[length(colors)]][[1]])){
            colors <- NULL
          }
        }
      } else {
        colors <- NULL
      }
      draw(plot_DiffEx(vals$counts, input$colorBar_Condition,
                       rownames(vals$diffexgenelist),
                       clusterRow = input$clusterRows,
                       clusterCol = input$clusterColumns,
                       displayRowLabels = input$displayHeatmapRowLabels,
                       displayColumnLabels = input$displayHeatmapColumnLabels,
                       displayRowDendrograms = input$displayHeatmapRowDendrograms,
                       displayColumnDendrograms = input$displayHeatmapColumnDendrograms,
                       annotationColors = colors,
                       columnTitle = input$heatmapColumnsTitle))
    }
  }, height = 600)

  #Render the interactive heatmap
  output$interactivediffPlot <- renderD3heatmap({
    if (!is.null(vals$diffexgenelist)){
      plot_d3DiffEx(vals$counts, input$selectDiffex_condition,
                    rownames(vals$diffexgenelist), clusterRow = input$clusterRows,
                    clusterCol = input$clusterColumns)
    }
  })

  #Create the differential expression results table
  output$diffextable <- renderDataTable({
    if (!is.null(vals$diffexgenelist)){
      temptable <- cbind(rownames(vals$diffexgenelist), data.frame(vals$diffexgenelist))
      colnames(temptable)[1] <- "Gene"
      temptable
    }
  })

  # Download the differential expression results table
  output$downloadGeneList <- downloadHandler(
    filename = function() {
      paste("diffex_results-", Sys.Date(), ".csv", sep = "")
    },
    content = function(file) {
      write.csv(vals$diffexgenelist, file)
    }
  )

  observeEvent(input$saveBiomarker, {
    if (is.null(input$biomarkerName)){
      alert("Warning: Specify biomarker name!")
    } else{
      rowData(vals$counts)[, input$biomarkerName] <- ifelse(rownames(rowData(vals$counts)) %in% rownames(vals$diffexgenelist), 1, 0)
      updateSelectInput(session, "filteredFeature", choices = c("none", colnames(rowData(vals$counts))))
    }
  })

  #-----------------------------------------------------------------------------
  # Page 5: Subsampling
  #-----------------------------------------------------------------------------

  #Run subsampling analysis
  runDownsampler <- observeEvent(input$runSubsample, {
    if (is.null(vals$counts)){
      alert("Warning: Upload data first!")
    }
    else{
      subData <- reactiveValues(
        counts = Downsample(assay(vals$counts, "counts"),
                            newcounts = floor(2 ^ seq.int(from = log2(input$minSim),
                                                          to = log2(input$maxSim),
                                                          length.out = 10)),
                            iterations = input$iterations)
      )
      output$downDone <- renderPlot({
        heatmap(as.matrix(subData$counts[order(apply(subData$counts[, , 10, 1], 1, sum), decreasing = TRUE)[1:20], , 10, 1]))
      })
    }
  })

  #Run differential power analysis
  runDiffPower <- observeEvent(input$runDifferentialPower, {
    if (is.null(vals$counts)){
      alert("Warning: Upload data first!")
    }
    else{
      #    if(exists('subData$counts')){
      output$powerBoxPlot <- renderPlot({
        subData <- Downsample(assay(vals$counts, "counts"), newcounts = floor(2 ^ seq.int(from = log2(input$minSim), to = log2(input$maxSim), length.out = 10)), iterations = input$iterations)
        diffPower <- differentialPower(datamatrix = assay(vals$counts, "counts"), downmatrix = subData, conditions = colData(vals$counts)[, input$subCovariate], method = input$selectDiffMethod)
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

  #-----------------------------------------------------------------------------
  # Page 6: MAST
  #-----------------------------------------------------------------------------

  #For conditions with more than two factors, select the factor of interest
  output$hurdleconditionofinterestUI <- renderUI({
    if (!is.null(vals$counts)){
      if (length(unique(colData(vals$counts)[, input$hurdlecondition])) > 2){
        selectInput("hurdleconditionofinterest",
                    "Select Factor of Interest",
                    unique(sort(colData(vals$counts)[, input$hurdlecondition])))
      }
    }
  })

  #Run MAST differential expression
  observeEvent(input$runDEhurdle, {
    if (is.null(vals$counts)){
      alert("Warning: Upload data first!")
    }
    else{
      withBusyIndicatorServer("runDEhurdle", {
        #run diffex to get gene list and pvalues
        vals$mastgenelist <- MAST(vals$counts,
                                  condition = input$hurdlecondition,
                                  interest.level = input$hurdleconditionofinterest,
                                  freq_expressed = input$hurdlethresh,
                                  fc_threshold = input$FCthreshold,
                                  p.value = input$hurdlepvalue,
                                  usethresh = input$useAdaptThresh)
      })
    }
  })

  observeEvent(input$runThreshPlot, {
    if (is.null(vals$counts)){
      alert("Warning: Upload data first!")
    }
    else{
      withBusyIndicatorServer("runThreshPlot", {
        output$threshplot <- renderPlot({
          vals$thres <- thresholdGenes(vals$counts)
          par(mfrow = c(5, 4))
          plot(vals$thres)
          par(mfrow = c(1, 1))
        }, height = 600)
      })
    }
  })

  #TODO: Choose assay here
  output$hurdleviolin <- renderPlot({
    if (!(is.null(vals$mastgenelist))){
      MASTviolin(SCEdata = vals$counts, use_assay = "logcounts",
                 fcHurdleSig = vals$mastgenelist,
                 variable = input$hurdlecondition)
    }
  }, height = 600)

  #TODO: Choose assay here
  output$hurdlelm <- renderPlot({
    if (!(is.null(vals$mastgenelist))){
      MASTregression(SCEdata = vals$counts, use_assay = "logcounts",
                     fcHurdleSig = vals$mastgenelist,
                     variable = input$hurdlecondition)
    }
  }, height = 600)

  output$hurdleHeatmap <- renderPlot({
    if (!(is.null(vals$mastgenelist))){
      draw(plot_DiffEx(vals$counts, input$hurdlecondition,
                       rownames(vals$mastgenelist),
                       annotationColors = c("red","green")))
    }
  }, height = 600)

  #Create the MAST results table
  output$mastresults <- renderDataTable({
    if (!is.null(vals$mastgenelist)){
      vals$mastgenelist
    }
  })

  #download mast results
})
