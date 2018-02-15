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
    combatstatus = "",
    diffexgenelist = NULL,
    gsva_res = NULL,
    gsva_limma = NULL
  )

  #Update all of the columns that depend on pvals columns
  updateColDataNames <- function(){
    pdata_options <- colnames(colData(vals$counts))
    updateSelectInput(session, "filteredSample",
                      choices = c("none", pdata_options))
    updateSelectInput(session, "deleterowdatacolumn",
                      choices = pdata_options)
    updateSelectInput(session, "colorBy",
                      choices = c("No Color", "Gene Expression", pdata_options))
    updateSelectInput(session, "shapeBy",
                      choices = c("No Shape", pdata_options))
    updateSelectInput(session, "selectDiffex_condition",
                      choices = pdata_options)
    updateSelectInput(session, "subCovariate",
                      choices = pdata_options)
    updateSelectInput(session, "combatBatchVar",
                      choices = pdata_options)
    updateSelectInput(session, "combatConditionVar",
                      choices = pdata_options)
    updateSelectInput(session, "hurdlecondition",
                      choices = pdata_options)
    updateSelectInput(session, "pathwayPlotVar",
                      choices = pdata_options)
    updateSelectInput(session, "select_ReadDepth_Condition",
                      choices = pdata_options)
    updateSelectInput(session, "select_CellNum_Condition",
                      choices = pdata_options)
    updateSelectInput(session, "select_Snapshot_Condition",
                      choices = pdata_options)
  }

  updateGeneNames <- function(){
    selectthegenes <- rownames(vals$counts)
    updateSelectizeInput(session, "colorGenes",
                         choices = selectthegenes, server = TRUE)
  }

  updateFeatureAnnots <- function(){
    updateSelectInput(session, "filteredFeature",
                      choices = c("none", colnames(rowData(vals$counts))))
  }

  updateNumSamples <- function(){
    numsamples <- ncol(vals$counts)
    updateSelectInput(session, "Knumber",
                      choices = 1:numsamples)
    updateSelectInput(session, "Cnumber",
                      choices = 1:numsamples)
    updateSelectInput(session, "pcX",
                      choices = paste("PC", 1:numsamples, sep = ""),
                      selected = "PC1")
    updateSelectInput(session, "pcY",
                      choices = paste("PC", 1:numsamples, sep = ""),
                      selected = "PC2")
    updateNumericInput(session, "downsampleNum", value = numsamples,
                       max = numsamples)
  }

  updateAssayInputs <- function(){
    currassays <- names(assays(vals$counts))
    updateSelectInput(session, "dimRedAssaySelect", choices = currassays)
    updateSelectInput(session, "combatAssay", choices = currassays)
    updateSelectInput(session, "diffexAssay", choices = currassays)
    updateSelectInput(session, "mastAssay", choices = currassays)
    updateSelectInput(session, "pathwayAssay", choices = currassays)
    updateSelectInput(session, "delAssayType", choices = currassays)
    updateSelectInput(session, "filterAssaySelect", choices = currassays)
  }

  updateReddimInputs <- function(){
    currreddim <- names(reducedDims(vals$counts))
    updateSelectInput(session, "delRedDimType", choices = currreddim)
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
      updateColDataNames()
      updateNumSamples()
      updateAssayInputs()
      updateGeneNames()
      updateReddimInputs()
      updateSelectInput(session, "deletesamplelist",
                        choices = colnames(vals$counts))
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
  # Page 2: Data Summary and Filtering
  #-----------------------------------------------------------------------------

  #Render data table if there are fewer than 50 samples
  output$contents <- renderDataTable({
    if (!is.null(getShinyOption("inputSCEset"))){
      updateGeneNames()
    }
    if (!(is.null(vals$counts)) && ncol(vals$counts) < 50){
      temptable <- cbind(rownames(vals$counts), assay(vals$counts, input$filterAssaySelect))
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
      plot_ly(x = apply(assay(vals$counts, input$filterAssaySelect), 2, function(x) sum(x)),
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
      plot_ly(x = apply(assay(vals$counts, input$filterAssaySelect), 2,
                        function(x) sum(x > 0)), type = "histogram") %>%
        layout(xaxis = x, yaxis = y)
    } else {
      plotly_empty(type = "scatter") %>% add_trace(mode = "lines")
      }
  })

  #random downsample of samples
  observeEvent(input$downsampleGo, {
    req(vals$counts)
    withBusyIndicatorServer("downsampleGo", {
      vals$counts <- vals$counts[, sample(ncol(vals$counts), input$downsampleNum)]
      updateNumSamples()
      vals$diffexgenelist <- NULL
      vals$gsva_res <- NULL
    })
  })

  #Render summary table
  output$summarycontents <- renderTable({
    req(vals$counts)
    summarizeTable(indata = vals$counts,
                   use_assay = input$filterAssaySelect,
                   expression_cutoff = input$minDetectGene)
  })

  #Filter the data based on the options
  observeEvent(input$filterData, {
    if (is.null(vals$original)){
      alert("Warning: Upload data first!")
    }
    else{
      withBusyIndicatorServer("filterData", {
        deletesamples <- input$deletesamplelist
        vals$counts <- filterSCData(insceset = vals$counts,
                                    use_assay = input$filterAssaySelect,
                                    deletesamples = deletesamples,
                                    remove_noexpress = input$removeNoexpress,
                                    remove_bottom = 0.01 * input$LowExpression,
                                    minimum_detect_genes = input$minDetectGene) #TODO: user decides to filter spikeins
        vals$diffexgenelist <- NULL
        vals$gsva_res <- NULL
        #Refresh things for the clustering tab
        updateGeneNames()
        if (!is.null(input$deletesamplelist)){
          updateSelectInput(session, "deletesamplelist",
                            choices = colnames(vals$counts))
        }
      })
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
      vals$diffexgenelist <- NULL
      vals$gsva_res <- NULL
      #Refresh things for the clustering tab
      updateColDataNames()
      updateNumSamples()
      updateAssayInputs()
      updateGeneNames()
    }
  })

  #Delete a column from the colData annotations
  observeEvent(input$deleterowDatabutton, {
    if (is.null(vals$original)){
      alert("Warning: Upload data first!")
    }
    else{
      colData(vals$counts) <- colData(vals$counts)[, !(colnames(colData(vals$counts)) %in% input$deleterowdatacolumn), drop = F]
      updateColDataNames()
    }
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
          L[[3]] <- list(withBusyIndicatorUI(actionButton("runFilterSample", "Filter")))
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
    withBusyIndicatorServer("runFilterSample", {
      filter <- colData(vals$counts)[, input$filteredSample] %in% input$filterSampleChoices
      vals$counts <- vals$counts[, filter]
      vals$diffexgenelist <- NULL
      vals$gsva_res <- NULL
      updateNumSamples()
    })
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

  observeEvent(input$orgOrganism, {
    library(input$orgOrganism, character.only = TRUE)
    indb <- get(paste(input$orgOrganism))
    output$orgConvertColumns <- renderUI({
      tagList(
        selectInput("orgFromCol", "Select From Annotation:", columns(indb)),
        selectInput("orgToCol", "Select To Annotation:", columns(indb))
      )
    })
  })

  observeEvent(input$convertGenes, {
    withBusyIndicatorServer("convertGenes", {
      vals$counts <- convert_gene_ids(inSCESet = vals$counts,
                                      in_symbol = input$orgFromCol,
                                      out_symbol = input$orgToCol,
                                      database = input$orgOrganism)
      updateGeneNames()
      vals$diffexgenelist <- NULL
      vals$gsva_res <- NULL
    })
  })

  #Filter the selected features
  observeEvent(input$runFilterFeature, {
    filter <- rowData(vals$counts)[, input$filteredFeature] %in% input$filterFeatureChoices
    vals$counts <- vals$counts[filter, ]
    updateGeneNames()
    vals$diffexgenelist <- NULL
    vals$gsva_res <- NULL
  })

  output$downloadSCE <- downloadHandler(
    filename <- function() {
      paste("SCE-", Sys.Date(), ".rds", sep = "")
    },
    content <- function(file) {
      saveRDS(vals$counts, file)
  })

  output$assayList <- renderTable(
    if (!is.null(vals$counts) && length(names(assays(vals$counts))) > 0){
      data.table(assays = names(assays(vals$counts)))
    }
  )

  output$reducedDimsList <- renderTable(
    if (!is.null(vals$counts) && length(names(reducedDims(vals$counts))) > 0){
      data.table("Reduced Dimension" = names(reducedDims(vals$counts)))
    }
  )

  observeEvent(input$addAssay, {
    req(vals$counts)
    if (input$addAssayType %in% names(assays(vals$counts))){
      alert("assay already exists!")
    } else {
      withBusyIndicatorServer("addAssay", {
        if (input$addAssayType == "logcounts"){
          if("counts" %in% names(assays(vals$counts))){
            assay(vals$counts, "logcounts") <- log2(assay(vals$counts, "counts") + 1)
          } else {
            alert("A matrix named counts is required to calculate logcounts")
          }
        } else if (input$addAssayType == "cpm") {
          if("counts" %in% names(assays(vals$counts))){
            assay(vals$counts, "cpm") <- apply(assay(vals$counts, "counts"), 2, function(x) { x / (sum(x) / 1000000) })
          } else {
            alert("Count matrix required for cpm calculation")
          }
        } else if (input$addAssayType == "logcpm") {
          #try to calculate from cpm
          if ("cpm" %in% names(assays(vals$counts))){
            assay(vals$counts, "logcpm") <- log2(assay(vals$counts, "cpm") + 1)
          } else if ("counts" %in% names(assays(vals$counts))){
          #calculate from counts
            assay(vals$counts, "logcpm") <- log2(apply(assay(vals$counts, "counts"), 2, function(x) { x / (sum(x) / 1000000) }) + 1)
          } else {
            alert("Count matrix or cpm required for logcpm calculation")
          }
        } else {
          stop("unsupported assay type")
        }
        updateAssayInputs()
      })
    }
  })

  observeEvent(input$delAssay, {
    req(vals$counts)
    if (!(input$delAssayType %in% names(assays(vals$counts)))){
      alert("assay does not exist!")
    } else {
      withBusyIndicatorServer("delAssay", {
        assay(vals$counts, input$delAssayType) <- NULL
        updateAssayInputs()
      })
    }
  })

  observeEvent(input$delRedDim, {
    req(vals$counts)
    if (!(input$delRedDimType %in% names(reducedDims(vals$counts)))){
      alert("reducedDim does not exist!")
    } else {
      withBusyIndicatorServer("delRedDim", {
        reducedDim(vals$counts, input$delRedDimType) <- NULL
        updateReddimInputs()
      })
    }
  })

  output$colDataDataFrame <- DT::renderDataTable({
    if (!is.null(vals$counts)){
      data.frame(colData(vals$counts))
    }
  }, options = list(scrollX = TRUE, pageLength = 30))

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
          vals$counts <- getPCA(count_data = vals$counts,
                                use_assay = input$dimRedAssaySelect,
                                reducedDimName = pcadimname)
          updateReddimInputs()
        }
        if (!is.null(reducedDim(vals$counts, pcadimname))){
          if (input$colorBy != "Gene Expression") {
            g <- singleCellTK::plotPCA(count_data = vals$counts,
                                       colorBy = input$colorBy,
                                       shape = input$shapeBy, pcX = input$pcX,
                                       pcY = input$pcY,
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
          vals$counts <- getTSNE(count_data = vals$counts,
                                 use_assay = input$dimRedAssaySelect,
                                 reducedDimName = tsnedimname)
          updateReddimInputs()
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
    if (is.null(vals$counts) | !(class(vals$counts) == "SCtkExperiment")){
    } else{
      if (input$dimRedPlotMethod == "PCA") {
        if (nrow(pca_variances(vals$counts)) == ncol(vals$counts)){
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
          vals$counts <- getPCA(count_data = vals$counts,
                                use_assay = input$dimRedAssaySelect,
                                reducedDimName = pcadimname)
          updateReddimInputs()
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
          vals$counts <- getTSNE(count_data = vals$counts,
                                 use_assay = input$dimRedAssaySelect,
                                 reducedDimName = tsnedimname)
          updateReddimInputs()
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
  }, height = 600)

  output$treePlot <- renderPlot({
    if (is.null(vals$counts)){
      alert("Warning: Upload data first!")
    } else {
      if (input$dimRedPlotMethod == "Dendrogram" && paste0("PCA", "_", input$dimRedAssaySelect) %in% names(reducedDims(vals$counts))){
        data <- getClusterInputData(count_data = vals$counts,
                                    inputData = "PCA Components",
                                    use_assay = input$dimRedAssaySelect,
                                    reducedDimName = paste0("PCA", "_", input$dimRedAssaySelect))
        d <- stats::dist(data)
        h <- stats::hclust(d, input$dendroDistanceMetric)
        if (input$clusteringAlgorithmD == "Phylogenetic Tree") {
          g <- ggtree(as.phylo(h), layout = "circular", open.angle = 360) + geom_tiplab2(size = 2)
        } else if (input$clusteringAlgorithmD == "Hierarchical") {
          g <- ggtree(as.phylo(h)) + theme_tree2() + geom_tiplab(size = 2)
        } else {
          stop("Input clustering algorithm not found ", input$clusteringAlgorithmD)
        }
        g
      }
    }
  }, height = 600)

  observeEvent(input$clusterData, {
    if (is.null(vals$counts)){
      alert("Warning: Upload data first!")
    } else if (input$clusterName == "") {
      alert("Cluster name required!")
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
          updateColDataNames()
          updateSelectInput(session, "colorBy",
                            choices = c("No Color", "Gene Expression",
                                        colnames(colData(vals$counts))),
                            selected = input$clusterName)
        } else if (input$clusteringAlgorithm == "Clara") {
          data <- getClusterInputData(count_data = vals$counts,
                                      inputData = input$selectClusterInputData,
                                      use_assay = input$dimRedAssaySelect,
                                      reducedDimName = currdimname)
          coutput <- clara(data, input$Cnumber)
          name <- input$clusterName
          colData(vals$counts)[, paste(name)] <- factor(coutput$clustering)
          updateColDataNames()
          updateSelectInput(session, "colorBy",
                            choices = c("No Color", "Gene Expression",
                                        colnames(colData(vals$counts))),
                            selected = input$clusterName)
        }
      })
    }
  })

  observe({
    if (!is.null(vals$original)){
      if(input$dimRedPlotMethod == "PCA"){
        pcadimname <- paste0("PCA", "_", input$dimRedAssaySelect)
        if (!is.null(reducedDim(vals$counts, pcadimname))) {
          curr_pcs <- colnames(reducedDim(vals$counts, "PCA_counts"))
          updateSelectInput(session, "pcX", choices = curr_pcs,
                            selected = curr_pcs[1])
          updateSelectInput(session, "pcY", choices = curr_pcs,
                            selected = curr_pcs[2])
        }
      }
    }
  })

  observeEvent(input$reRunTSNE, {
    if (is.null(vals$original)){
      alert("Warning: Upload data first!")
    }
    else{
      withBusyIndicatorServer("reRunTSNE", {
        vals$counts <- getTSNE(count_data = vals$counts,
                               use_assay = input$dimRedAssaySelect,
                               reducedDimName = paste0("TSNE", "_",
                                                       input$dimRedAssaySelect))
        updateReddimInputs()
      })
    }
  })

  observeEvent(input$reRunPCA, {
    if (is.null(vals$original)){
      alert("Warning: Upload data first!")
    }
    else{
      withBusyIndicatorServer("reRunPCA", {
        vals$counts <- getPCA(count_data = vals$counts,
                              use_assay = input$dimRedAssaySelect,
                              reducedDimName = paste0("PCA", "_",
                                                      input$dimRedAssaySelect))
        updateReddimInputs()
      })
    }
  })

  #-----------------------------------------------------------------------------
  # Page 4: Batch Correction
  #-----------------------------------------------------------------------------

  output$selectCombat_refbatchUI <- renderUI({
    if (!is.null(vals$counts)){
      if (input$combatRef){
        selectInput("combatRefBatch", "Choose Reference Batch:",
                    unique(sort(colData(vals$counts)[, input$combatBatchVar])))
      }
    }
  })

  observeEvent(input$combatRun, {
    if (is.null(vals$counts)){
      alert("Warning: Upload data first!")
    }
    else{
      withBusyIndicatorServer("combatRun", {
        if (input$batchMethod == "ComBat"){
          #check for zeros
          if (any(rowSums(assay(vals$counts, input$combatAssay)) == 0)){
            alert("Warning: Rows with a sum of zero found. Filter data to continue")
          } else {
            saveassayname <- gsub(" ", "_", input$combatSaveAssay)
            if (input$combatRef){
              assay(vals$counts, saveassayname) <-
                ComBat_SCE(SCEdata = vals$counts, batch = input$combatBatchVar,
                           use_assay = input$combatAssay,
                           par.prior = input$combatParametric,
                           covariates = input$combatConditionVar,
                           mean.only = input$combatMeanOnly,
                           ref.batch = input$combatRefBatch)
            } else {
              assay(vals$counts, saveassayname) <-
                ComBat_SCE(SCEdata = vals$counts, batch = input$combatBatchVar,
                           use_assay = input$combatAssay,
                           par.prior = input$combatParametric,
                           covariates = input$combatConditionVar,
                           mean.only = input$combatMeanOnly)
            }
            updateAssayInputs()
            vals$combatstatus <- "ComBat Complete"
          }
        } else {
          alert("Unsupported Batch Correction Method!")
        }
      })
    }
  })

  output$combatStatus <- renderUI({
    h2(vals$combatstatus)
  })

  output$combatBoxplot <- renderPlot({
    if (!is.null(vals$counts) &
        !is.null(input$combatBatchVar) &
        length(input$combatConditionVar) == 1){
      plotBatchVariance(inSCESet = vals$counts,
                        use_assay = input$combatAssay,
                        batch = input$combatBatchVar,
                        condition = input$combatConditionVar)
    }
  }, height = 600)

  #-----------------------------------------------------------------------------
  # Page 5.1: Differential Expression
  #-----------------------------------------------------------------------------

  output$selectDiffex_conditionUI <- renderUI({
    if (!is.null(vals$counts)){
      if (input$selectDiffex == "ANOVA") {
        tagList(
          selectInput("selectDiffex_condition", "Select Condition(s):",
                      colnames(colData(vals$counts)), multiple = TRUE)
        )
      } else {
        tagList(
          selectInput("selectDiffex_condition",
                      "Select Condition:",
                      colnames(colData(vals$counts))),
          selectInput("selectDiffex_covariates",
                      "Select Additional Covariates:",
                      colnames(colData(vals$counts)), multiple = TRUE)
        )
      }
    }
  })

  #For conditions with more than two factors, select the factor of interest
  output$selectDiffex_conditionlevelUI <- renderUI({
    if (!is.null(vals$counts) & length(colnames(colData(sce))) > 0){
      if (length(unique(colData(vals$counts)[, input$selectDiffex_condition])) > 2 & input$selectDiffex != "ANOVA"){
        tagList(
          conditionalPanel(
            condition = "input.selectDiffex == 'DESeq2'",
            radioButtons("selectDiffexConditionMethod", "Select Analysis Method:",
                         choiceNames = c("Biomarker (1 vs all)", "Factor of Interest vs. Control Factor"),
                         choiceValues = c("biomarker", "contrast"))
          ),
          selectInput("selectDiffex_conditionofinterest",
                      "Select Factor of Interest",
                      unique(sort(colData(vals$counts)[, input$selectDiffex_condition]))),
          conditionalPanel(
            condition = "input.selectDiffexConditionMethod == 'contrast' && input.selectDiffex == 'DESeq2'",
            selectInput("selectDiffex_controlcondition",
                        "Select Control Factor",
                        unique(sort(colData(vals$counts)[, input$selectDiffex_condition])))
          )
        )
      } else if (input$selectDiffex == "ANOVA") {
        tagList(
          selectInput("anovaCovariates", "Select Additional Covariates:",
                      colnames(colData(vals$counts)), multiple = TRUE)
        )
      }
    }
  })

  #Run differential expression
  observeEvent(input$runDiffex, {
    if (is.null(vals$counts)){
      alert("Warning: Upload data first!")
    }
    else{
      withBusyIndicatorServer("runDiffex", {
        #run diffex to get gene list and pvalues
        if(input$selectDiffex == "ANOVA"){
          use_covariates <- input$anovaCovariates
        } else {
          use_covariates <- input$selectDiffex_covariates
        }
        vals$diffexgenelist <- scDiffEx(inSCESet = vals$counts,
                                        use_assay = input$diffexAssay,
                                        condition = input$selectDiffex_condition,
                                        covariates = use_covariates,
                                        significance = input$selectPval,
                                        ntop = input$selectNGenes,
                                        usesig = input$applyCutoff,
                                        diffexmethod = input$selectDiffex,
                                        clusterRow = input$clusterRows,
                                        clusterCol = input$clusterColumns,
                                        levelofinterest = input$selectDiffex_conditionofinterest,
                                        analysis_type = input$selectDiffexConditionMethod,
                                        controlLevel = input$selectDiffex_controlcondition)
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
      draw(plot_DiffEx(inSCESet = vals$counts,
                       use_assay = input$diffexAssay,
                       condition = input$colorBar_Condition,
                       geneList = rownames(vals$diffexgenelist),
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
    if (input$biomarkerName == ""){
      alert("Warning: Specify biomarker name!")
    } else{
      withBusyIndicatorServer("saveBiomarker", {
        biomarker_name <- gsub(" ", "_", input$biomarkerName)
        rowData(vals$counts)[, biomarker_name] <- ifelse(rownames(vals$counts) %in% rownames(vals$diffexgenelist), 1, 0)
        updateFeatureAnnots()
      })
    }
  })

  #-----------------------------------------------------------------------------
  # Page 5.2: MAST
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
    } else {
      withBusyIndicatorServer("runDEhurdle", {
        #run diffex to get gene list and pvalues
        vals$mastgenelist <- MAST(SCEdata = vals$counts,
                                  use_assay = input$mastAssay,
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
          vals$thres <- thresholdGenes(SCEdata = vals$counts,
                                       use_assay = input$mastAssay)
          par(mfrow = c(5, 4))
          plot(vals$thres)
          par(mfrow = c(1, 1))
        }, height = 600)
      })
    }
  })

  output$hurdleviolin <- renderPlot({
    if (!(is.null(vals$mastgenelist))){
      MASTviolin(SCEdata = vals$counts, use_assay = input$mastAssay,
                 fcHurdleSig = vals$mastgenelist,
                 variable = input$hurdlecondition,
                 threshP = input$useAdaptThresh)
    }
  }, height = 600)

  output$hurdlelm <- renderPlot({
    if (!(is.null(vals$mastgenelist))){
      MASTregression(SCEdata = vals$counts, use_assay = input$mastAssay,
                     fcHurdleSig = vals$mastgenelist,
                     variable = input$hurdlecondition,
                     threshP = input$useAdaptThresh)
    }
  }, height = 600)

  output$hurdleHeatmap <- renderPlot({
    if (!(is.null(vals$mastgenelist))){
      draw(plot_DiffEx(vals$counts, use_assay = input$mastAssay,
                       condition = input$hurdlecondition,
                       geneList = vals$mastgenelist$Gene,
                       annotationColors = "auto",
                       columnTitle = "MAST"))
    }
  }, height = 600)

  #Create the MAST results table
  output$mastresults <- renderDataTable({
    if (!is.null(vals$mastgenelist)){
      vals$mastgenelist
    }
  })

  #download mast results
  output$downloadHurdleResult <- downloadHandler(
    filename = function() {
      paste("mast_results-", Sys.Date(), ".csv", sep = "")
    },
    content = function(file) {
      write.csv(vals$mastgenelist, file)
    }
  )

  #-----------------------------------------------------------------------------
  # Page 6: Pathway Activity Analysis
  #-----------------------------------------------------------------------------

  output$selectPathwayGeneLists <- renderUI({
    if (input$genelistSource == "Manual Input"){
      if (!is.null(vals$counts)){
        selectizeInput("pathwayGeneLists", "Select Gene List(s):",
                       colnames(rowData(vals$counts)), multiple = T)
      } else {
        h4("Note: upload data.")
      }
    } else {
      selectInput("pathwayGeneLists", "Select Gene List(s):",
                  c("ALL", names(c2BroadSets)), multiple = T)
    }
  })

  output$selectNumTopPaths <- renderUI({
    if (!is.null(input$pathwayGeneLists) && input$pathwayGeneLists == "ALL" && input$genelistSource == "MSigDB c2 (Human, Entrez ID only)"){
      sliderInput("pickNtopPaths", "Number of top pathways:", min = 5,
                  max = length(c2BroadSets), value = 25, step = 5)
    }
  })

  observeEvent(input$pathwayRun, {
    if (is.null(vals$counts)){
      alert("Warning: Upload data first!")
    } else {
      withBusyIndicatorServer("pathwayRun", {
        vals$gsva_res <- GSVA_sce(SCEdata = vals$counts,
                                  use_assay = input$pathwayAssay,
                                  pathway_source = input$genelistSource,
                                  pathway_names = input$pathwayGeneLists)
      })
    }
  })

  observe({
    if(length(input$pathwayPlotVar) == 1 & !(is.null(vals$gsva_res))){
      fit <- limma::lmFit(vals$gsva_res, stats::model.matrix(~factor(colData(vals$counts)[, input$pathwayPlotVar])))
      fit <- limma::eBayes(fit)
      toptableres <- limma::topTable(fit, number = nrow(vals$gsva_res))
      temptable <- cbind(rownames(toptableres), toptableres)
      rownames(temptable) <- NULL
      colnames(temptable)[1] <- "Pathway"
      vals$gsva_limma <- temptable
    } else {
      vals$gsva_limma <- NULL
    }
  })

  output$pathwaytable <- DT::renderDataTable({
    if (!is.null(vals$gsva_limma)){
      if (!is.null(input$pathwayGeneLists) && input$pathwayGeneLists == "ALL" && input$genelistSource == "MSigDB c2 (Human, Entrez ID only)"){
        vals$gsva_limma[1:min(input$pickNtopPaths, nrow(vals$gsva_limma)), , drop = F]
      } else {
        vals$gsva_limma
      }
    }
  }, options = list(scrollX = TRUE, pageLength = 30))

  output$pathwayPlot <- renderPlot({
    if (!(is.null(vals$gsva_res))){
      if(input$genelistSource == "MSigDB c2 (Human, Entrez ID only)" & "ALL" %in% input$pathwayGeneLists & !(is.null(vals$gsva_limma))){
        tempgsvares <- vals$gsva_res[vals$gsva_limma$Pathway[1:min(input$pickNtopPaths, nrow(vals$gsva_limma))], , drop = F]
      } else {
        tempgsvares <- vals$gsva_res
      }
      if (input$pathwayOutPlot == "Violin" && length(input$pathwayPlotVar) > 0){
        tempgsvares <- tempgsvares[1:min(49, input$pickNtopPaths, nrow(tempgsvares)), , drop = F]
      }
      GSVA_plot(SCEdata = vals$counts,
                gsva_data = tempgsvares,
                plot_type = input$pathwayOutPlot,
                condition = input$pathwayPlotVar)
    }
  })

  #save pathawy activity results in the colData
  observeEvent(input$savePathway, {
    if (!(is.null(vals$gsva_res))){
      if (all(colnames(vals$counts) == colnames(vals$gsva_res))){
        colData(vals$counts) <- cbind(colData(vals$counts), data.frame(t(vals$gsva_res)))
        updateColDataNames()
      }
    } else {
      alert("Run pathway first!")
    }
  })

  #download mast results
  output$downloadPathway <- downloadHandler(
    filename = function() {
      paste("pathway_results-", Sys.Date(), ".csv", sep = "")
    },
    content = function(file) {
      write.csv(vals$gsva_res, file)
    }
  )

  #-----------------------------------------------------------------------------
  # Page 7: Subsampling
  #-----------------------------------------------------------------------------

  #Run subsampling analysis
  observeEvent(input$runSubsampleDepth, {
    if (is.null(vals$counts)){
      alert("Warning: Upload data first!")
    }
    else{
      withBusyIndicatorServer("runSubsampleDepth", {
        vals$subDepth <- DownsampleDepth(originalData = vals$counts,
                                         minCount = input$minCount,
                                         minCells = input$minCells,
                                         maxDepth = 10 ^ input$maxDepth,
                                         realLabels = input$select_ReadDepth_Condition,
                                         depthResolution = input$depthResolution,
                                         iterations = input$iterations)

        output$DepthDone <- renderPlot({
          plot(apply(vals$subDepth[, , 1], 2, median)~
                 seq(from = 0, to = input$maxDepth, length.out = input$depthResolution),
               lwd = 4, xlab = "log10(Total read counts)", ylab = "Number of detected genes",
               main = "Number of dected genes by sequencing depth")
          lines(apply(vals$subDepth[, , 1], 2, function(x){quantile(x, 0.25)})~
                  seq(from = 0, to = input$maxDepth, length.out = input$depthResolution), lty = 2, lwd = 3)
          lines(apply(vals$subDepth[, , 1], 2, function(x){quantile(x, 0.25)})~
                  seq(from = 0, to = input$maxDepth, length.out = input$depthResolution), lty = 2, lwd = 3)
        })
        output$MinEffectDone <- renderPlot({
          plot(apply(vals$subDepth[, , 2], 2, median)~
                 seq(from = 0, to = input$maxDepth, length.out = input$depthResolution),
               lwd = 4, xlab = "log10(Total read counts)", ylab = "Average significant effect size",
               ylim = c(0, 2))
          lines(apply(vals$subDepth[, , 2], 2, function(x){quantile(x, 0.25)})~
                  seq(from = 0, to = input$maxDepth, length.out = input$depthResolution), lty = 2, lwd = 3)
          lines(apply(vals$subDepth[, , 2], 2, function(x){quantile(x, 0.75)})~
                  seq(from = 0, to = input$maxDepth, length.out = input$depthResolution), lty = 2, lwd = 3)
        })
        output$sigNumDone <- renderPlot({
          plot(apply(vals$subDepth[, , 3], 2, median)~
                 seq(from = 0, to = input$maxDepth, length.out = input$depthResolution),
               lwd = 4, xlab = "log10(Total read counts)", ylab = "Number of significantly DiffEx genes")
          lines(apply(vals$subDepth[, , 3], 2, function(x){quantile(x, 0.25)})~
                  seq(from = 0, to = input$maxDepth, length.out = input$depthResolution), lty = 2, lwd = 3)
          lines(apply(vals$subDepth[, , 3], 2, function(x){quantile(x, 0.75)})~
                  seq(from = 0, to = input$maxDepth, length.out = input$depthResolution), lty = 2, lwd = 3)
        })
      })
    }
  })

  observeEvent(input$runSubsampleCells, {
    if (is.null(vals$counts)){
      alert("Warning: Upload data first!")
    }
    else{
      withBusyIndicatorServer("runSubsampleCells", {
        if(input$useReadCount){
          vals$subCells <- DownsampleCells(originalData = vals$counts,
                                           realLabels = input$select_CellNum_Condition,
                                           totalReads = sum(counts(vals$counts)),
                                           minCellnum = input$minCellNum,
                                           maxCellnum = input$maxCellNum,
                                           minCountDetec = input$minCount,
                                           minCellsDetec = input$minCells,
                                           depthResolution = input$depthResolution,
                                           iterations = input$iterations)
        }
        else{
          vals$subCells <- DownsampleCells(originalData = vals$counts,
                                           realLabels = input$select_CellNum_Condition,
                                           totalReads = input$totalReads,
                                           minCellnum = input$minCellNum,
                                           maxCellnum = input$maxCellNum,
                                           minCountDetec = input$minCount,
                                           minCellsDetec = input$minCells,
                                           depthResolution = input$depthResolution,
                                           iterations = input$iterations)
        }
        output$CellsDone <- renderPlot({
          plot(apply(vals$subCells[, , 1], 2, median)~
                 seq(from = input$minCellNum, to = input$maxCellNum, length.out = input$depthResolution),
               lwd = 4, xlab = "Number of virtual cells", ylab = "Number of detected genes",
               main = "Number of dected genes by cell number")
          lines(apply(vals$subCells[, , 1], 2, function(x){quantile(x, 0.25)})~
                  seq(from = input$minCellNum, to = input$maxCellNum, length.out = input$depthResolution), lty = 2, lwd = 3)
          lines(apply(vals$subCells[, , 1], 2, function(x){quantile(x, 0.25)})~
                  seq(from = input$minCellNum, to = input$maxCellNum, length.out = input$depthResolution), lty = 2, lwd = 3)
        })
        output$MinEffectCells <- renderPlot({
          plot(apply(vals$subCells[, , 2], 2, median)~
                 seq(from = input$minCellNum, to = input$maxCellNum, length.out = input$depthResolution),
               lwd = 4, xlab = "Number of virtual cells", ylab = "Average significant effect size",
               ylim = c(0, 2))
          lines(apply(vals$subCells[, , 2], 2, function(x){quantile(x, 0.25)})~
                  seq(from = input$minCellNum, to = input$maxCellNum, length.out = input$depthResolution), lty = 2, lwd = 3)
          lines(apply(vals$subCells[, , 2], 2, function(x){quantile(x, 0.75)})~
                  seq(from = input$minCellNum, to = input$maxCellNum, length.out = input$depthResolution), lty = 2, lwd = 3)
        })
        output$sigNumCells <- renderPlot({
          plot(apply(vals$subCells[, , 3], 2, median)~
                 seq(from = input$minCellNum, to = input$maxCellNum, length.out = input$depthResolution),
               lwd = 4, xlab = "Number of vitual cells", ylab = "Number of significantly DiffEx genes")
          lines(apply(vals$subCells[, , 3], 2, function(x){quantile(x, 0.25)})~
                  seq(from = input$minCellNum, to = input$maxCellNum, length.out = input$depthResolution), lty = 2, lwd = 3)
          lines(apply(vals$subCells[, , 3], 2, function(x){quantile(x, 0.75)})~
                  seq(from = input$minCellNum, to = input$maxCellNum, length.out = input$depthResolution), lty = 2, lwd = 3)
        })
      })
    }
  })

  #Run differential power analysis
  observeEvent(input$runSnapshot, {
    if (is.null(vals$counts)){
      alert("Warning: Upload data first!")
    }
    else{
      withBusyIndicatorServer("runSnapshot", {
        vals$snapshot <- iterateSimulations(originalData = vals$counts,
                                            realLabels = input$select_Snapshot_Condition,
                                            totalReads = input$numReadsSnap,
                                            cells = input$numCellsSnap,
                                            iterations = input$iterationsSnap)
        vals$effectSizes <- calcEffectSizes(countMatrix = counts(vals$counts), condition = colData(vals$counts)[, input$select_Snapshot_Condition])
        output$Snaplot <- renderPlot({
          plot(apply(vals$snapshot, 1, function(x){sum(x <= 0.05) / length(x)}) ~ vals$effectSizes,
               xlab = "Cohen's d effect size", ylab = "Detection power", lwd = 4, main = "Power to detect diffex by effect size")
        })
      })
    }
  })
})
