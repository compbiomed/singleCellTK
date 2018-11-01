#1GB max upload size
options(shiny.maxRequestSize = 1000 * 1024 ^ 2)

internetConnection <- suppressWarnings(Biobase::testBioCConnection())

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
    gsvaRes = NULL,
    gsvaLimma = NULL,
    visplotobject = NULL,
    enrichRes = NULL,
    diffexheatmapplot = NULL,
    diffexBmName = NULL
  )

  #Update all of the columns that depend on pvals columns
  updateColDataNames <- function(){
    pdataOptions <- colnames(colData(vals$counts))
    updateSelectInput(session, "filteredSample",
                      choices = c("none", pdataOptions))
    updateSelectInput(session, "deleterowdatacolumn",
                      choices = pdataOptions)
    updateSelectInput(session, "colorBy",
                      choices = c("No Color", "Gene Expression", pdataOptions))
    updateSelectInput(session, "shapeBy",
                      choices = c("No Shape", pdataOptions))
    updateSelectInput(session, "selectDiffexCondition",
                      choices = pdataOptions)
    updateSelectInput(session, "subCovariate",
                      choices = pdataOptions)
    updateSelectInput(session, "batchVarPlot",
                      choices = c("none", pdataOptions))
    updateSelectInput(session, "conditionVarPlot",
                      choices = c("none", pdataOptions))
    updateSelectInput(session, "combatBatchVar",
                      choices = pdataOptions)
    updateSelectInput(session, "combatConditionVar",
                      choices = pdataOptions)
    updateSelectInput(session, "hurdlecondition",
                      choices = pdataOptions)
    updateSelectInput(session, "pathwayPlotVar",
                      choices = pdataOptions)
    updateSelectInput(session, "selectReadDepthCondition",
                      choices = pdataOptions)
    updateSelectInput(session, "selectCellNumCondition",
                      choices = pdataOptions)
    updateSelectInput(session, "selectSnapshotCondition",
                      choices = pdataOptions)
    updateSelectInput(session, "annotModifyChoice",
                      choices = c("none", pdataOptions))
    updateSelectInput(session, "visCondn",
                      choices = c("none", pdataOptions))
  }

  updateGeneNames <- function(){
    selectthegenes <- rownames(vals$counts)
    updateSelectizeInput(session, "colorGenes",
                         choices = selectthegenes, server = TRUE)
    updateSelectizeInput(session, "selectvisGenes",
                         choices = selectthegenes, server = TRUE)
    updateSelectizeInput(session, "enrichGenes",
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
    updateSelectInput(session, "modifyAssaySelect", choices = currassays)
    updateSelectInput(session, "filterAssaySelect", choices = currassays)
    updateSelectInput(session, "visAssaySelect", choices = currassays)
    updateSelectInput(session, "enrichAssay", choices = currassays)
    updateSelectInput(session, "depthAssay", choices = currassays)
    updateSelectInput(session, "cellsAssay", choices = currassays)
    updateSelectInput(session, "snapshotAssay", choices = currassays)
  }

  updateReddimInputs <- function(){
    currreddim <- names(reducedDims(vals$counts))
    updateSelectInput(session, "delRedDimType", choices = currreddim)
  }

  updateEnrichDB <- function(){
    if (internetConnection){
      enrDB <- enrichR::listEnrichrDbs()$libraryName
    } else {
      enrDB <- ""
    }
    updateSelectInput(session, "enrichDb", choices = c("ALL", enrDB))
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
        vals$original <- createSCE(assayFile = input$countsfile$datapath,
                                   annotFile = input$annotFile$datapath,
                                   featureFile = input$featureFile$datapath,
                                   assayName = input$inputAssayType,
                                   createLogCounts = input$createLogcounts)
      } else if (input$uploadChoice == "example"){
        if (input$selectExampleData == "mouseBrainSubset"){
          data(list = paste0(input$selectExampleData, "SCE"))
          vals$original <- base::eval(parse(text = paste0(input$selectExampleData, "SCE")))
        } else if (input$selectExampleData == "maits"){
          data(maits, package = "MAST")
          vals$original <- createSCE(assayFile = t(maits$expressionmat),
                                     annotFile = maits$cdat,
                                     featureFile = maits$fdat,
                                     assayName = "logtpm",
                                     inputDataFrames = TRUE,
                                     createLogCounts = FALSE)
          rm(maits)
        } else if (input$selectExampleData == "fluidigm_pollen_et_al") {
          data(fluidigm, package = "scRNAseq")
          tempsce <- as(fluidigm, "SingleCellExperiment")
          vals$original <- as(tempsce, "SCtkExperiment")
          rm(fluidigm, tempsce)
        } else if (input$selectExampleData == "th2_mahata_et_al") {
          data(th2, package = "scRNAseq")
          tempsce <- as(th2, "SingleCellExperiment")
          vals$original <- as(tempsce, "SCtkExperiment")
          rm(th2, tempsce)
        } else if (input$selectExampleData == "allen_tasic_et_al") {
          data(allen, package = "scRNAseq")
          tempsce <- as(allen, "SingleCellExperiment")
          vals$original <- as(tempsce, "SCtkExperiment")
          rm(allen, tempsce)
        }
      } else if (input$uploadChoice == "rds") {
        importedrds <- readRDS(input$rdsFile$datapath)
        if (methods::is(importedrds, "SummarizedExperiment")) {
          vals$original <- importedrds
        } else {
          vals$original <- NULL
        }
      }
      if (!is.null(vals$original)) {
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
          ui = tags$div(
            class = "alert alert-success alert-dismissible",
            HTML("<span class='glyphicon glyphicon-ok' aria-hidden='true'> \
                  </span> Successfully Uploaded! <button type='button' \
                  class='close' data-dismiss='alert'>&times;</button>"))
        )
      } else {
        shinyalert::shinyalert("Error!", "The data upload failed!",
                               type = "error")
      }
      vals$diffexheatmapplot <- NULL
      vals$combatstatus <- ""
      vals$diffexgenelist <- NULL
      vals$gsvaRes <- NULL
      vals$gsvaLimma <- NULL
      vals$visplotobject <- NULL
      vals$enrichRes <- NULL
      vals$diffexheatmapplot <- NULL
      vals$diffexBmName <- NULL
      myValues$dList <- NULL
    })
  })

  #-----------------------------------------------------------------------------
  # Page 2: Data Summary and Filtering
  #-----------------------------------------------------------------------------

  #Sidebar buttons functionality - not an accordion
  shinyjs::onclick("f_hideAllSections", allSections(
    "hide", c(paste("f_collapse", 1:7, sep = ""))), add = TRUE)
  shinyjs::onclick("f_showAllSections", allSections(
    "show", c(paste("f_collapse", 1:7, sep = ""))), add = TRUE)
  shinyjs::onclick("f_button1", shinyjs::toggle(id = "f_collapse1",
                                                anim = TRUE), add = TRUE)
  shinyjs::onclick("f_button2", shinyjs::toggle(id = "f_collapse2",
                                                anim = TRUE), add = TRUE)
  shinyjs::onclick("f_button3", shinyjs::toggle(id = "f_collapse3",
                                                anim = TRUE), add = TRUE)
  shinyjs::onclick("f_button4", shinyjs::toggle(id = "f_collapse4",
                                                anim = TRUE), add = TRUE)
  shinyjs::onclick("f_button5", shinyjs::toggle(id = "f_collapse5",
                                                anim = TRUE), add = TRUE)
  shinyjs::onclick("f_button6", shinyjs::toggle(id = "f_collapse6",
                                                anim = TRUE), add = TRUE)
  shinyjs::onclick("f_button7", shinyjs::toggle(id = "f_collapse7",
                                                anim = TRUE), add = TRUE)
  #Button styling
  shinyjs::addClass(id = "f_button1", class = "btn-block")
  shinyjs::addClass(id = "f_button2", class = "btn-block")
  shinyjs::addClass(id = "f_button3", class = "btn-block")
  shinyjs::addClass(id = "f_button4", class = "btn-block")
  shinyjs::addClass(id = "f_button5", class = "btn-block")
  shinyjs::addClass(id = "f_button6", class = "btn-block")
  shinyjs::addClass(id = "f_button7", class = "btn-block")
  shinyjs::addClass(id = "filterData", class = "btn-block")
  shinyjs::addClass(id = "resetData", class = "btn-block")
  shinyjs::addClass(id = "convertGenes", class = "btn-block")
  shinyjs::addClass(id = "deleterowDatabutton", class = "btn-block")
  shinyjs::addClass(id = "downsampleGo", class = "btn-block")

  #Render data table if there are fewer than 50 samples
  output$contents <- DT::renderDataTable({
    req(vals$counts)
    if (!is.null(getShinyOption("inputSCEset"))){
      updateGeneNames()
    }
    if (!(is.null(vals$counts)) & ncol(vals$counts) < 50){
      temptable <- cbind(rownames(vals$counts), assay(vals$counts, input$filterAssaySelect))
      colnames(temptable)[1] <- "Gene"
      temptable
    }
  }, options = list(scrollX = TRUE), rownames = FALSE)

  #Render histogram of read counts per cell
  output$countshist <- renderPlotly({
    if (!(is.null(vals$counts))){
      f <- list(family = "Arial", size = 14, color = "#7f7f7f")
      x <- list(title = "Reads per cell", titlefont = f)
      y <- list(title = "Number of cells", titlefont = f)
      plotly::plot_ly(x = apply(assay(vals$counts, input$filterAssaySelect), 2, function(x) sum(x)),
                      type = "histogram") %>%
        plotly::layout(xaxis = x, yaxis = y)
    } else {
      plotly::plotly_empty(type = "scatter") %>% plotly::add_trace(mode = "lines")
    }
  })

  #Render histogram of genes detected per cell
  output$geneshist <- renderPlotly({
    if (!(is.null(vals$counts))){
      f <- list(family = "Arial", size = 14, color = "#7f7f7f")
      x <- list(title = "Genes detected per cell", titlefont = f)
      y <- list(title = "Number of cells", titlefont = f)
      plotly::plot_ly(x = apply(assay(vals$counts, input$filterAssaySelect), 2,
                                function(x) sum(x > 0)), type = "histogram") %>%
        plotly::layout(xaxis = x, yaxis = y)
    } else {
      plotly::plotly_empty(type = "scatter") %>% plotly::add_trace(mode = "lines")
    }
  })

  #random downsample of samples
  observeEvent(input$downsampleGo, {
    req(vals$counts)
    withBusyIndicatorServer("downsampleGo", {
      vals$counts <- vals$counts[, sample(ncol(vals$counts), input$downsampleNum)]
      updateNumSamples()
      vals$diffexheatmapplot <- NULL
      vals$combatstatus <- ""
      vals$diffexgenelist <- NULL
      vals$gsvaRes <- NULL
      vals$gsvaLimma <- NULL
      vals$visplotobject <- NULL
      vals$enrichRes <- NULL
      vals$diffexBmName <- NULL
      myValues$dList <- NULL
    })
  })

  #Render summary table
  output$summarycontents <- renderTable({
    req(vals$counts)
    singleCellTK::summarizeTable(inSCE = vals$counts,
                                 useAssay = input$filterAssaySelect,
                                 expressionCutoff = input$minDetectGene)
  })

  #Filter the data based on the options
  observeEvent(input$filterData, {
    if (is.null(vals$original)){
      shinyalert::shinyalert("Error!", "Upload data first.", type = "error")
    }
    else{
      withBusyIndicatorServer("filterData", {
        deletesamples <- input$deletesamplelist
        vals$counts <- filterSCData(inSCE = vals$counts,
                                    useAssay = input$filterAssaySelect,
                                    deletesamples = deletesamples,
                                    removeNoExpress = input$removeNoexpress,
                                    removeBottom = 0.01 * input$LowExpression,
                                    minimumDetectGenes = input$minDetectGene) #TODO: user decides to filter spikeins
        vals$diffexheatmapplot <- NULL
        vals$combatstatus <- ""
        vals$diffexgenelist <- NULL
        vals$gsvaRes <- NULL
        vals$gsvaLimma <- NULL
        vals$visplotobject <- NULL
        vals$enrichRes <- NULL
        vals$diffexBmName <- NULL
        myValues$dList <- NULL
        #Refresh things for the clustering tab
        updateGeneNames()
        updateEnrichDB()
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
      shinyalert::shinyalert("Error!", "Upload data first.", type = "error")
    }
    else{
      vals$counts <- vals$original
      updateSelectInput(session, "deletesamplelist",
                        choices = colnames(vals$counts))
      vals$diffexheatmapplot <- NULL
      vals$combatstatus <- ""
      vals$diffexgenelist <- NULL
      vals$gsvaRes <- NULL
      vals$gsvaLimma <- NULL
      vals$visplotobject <- NULL
      vals$enrichRes <- NULL
      vals$diffexBmName <- NULL
      myValues$dList <- NULL
      #Refresh things for the clustering tab
      updateColDataNames()
      updateNumSamples()
      updateAssayInputs()
      updateGeneNames()
      updateEnrichDB()
    }
  })

  #Delete a column from the colData annotations
  observeEvent(input$deleterowDatabutton, {
    if (is.null(vals$original)){
      shinyalert::shinyalert("Error!", "Upload data first.", type = "error")
    }
    else{
      colData(vals$counts) <- colData(vals$counts)[, !(colnames(colData(vals$counts)) %in% input$deleterowdatacolumn), drop = FALSE]
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
                                                      choices = unique(colData(vals$counts)[, input$filteredSample]))),
                              tags$h5(tags$i("Note: the Reset button is in 'Delete Outliers' tab above."))
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
      vals$gsvaRes <- NULL
      vals$enrichRes <- NULL
      vals$visplotobject <- NULL
      vals$diffexheatmapplot <- NULL
      vals$combatstatus <- ""
      vals$gsvaLimma <- NULL
      vals$diffexBmName <- NULL
      myValues$dList <- NULL
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
    req(vals$counts)
    withBusyIndicatorServer("convertGenes", {
      vals$counts <- convertGeneIDs(inSCE = vals$counts,
                                    inSymbol = input$orgFromCol,
                                    outSymbol = input$orgToCol,
                                    database = input$orgOrganism)
      updateGeneNames()
      vals$diffexgenelist <- NULL
      vals$gsvaRes <- NULL
      vals$enrichRes <- NULL
      vals$visplotobject <- NULL
      vals$diffexheatmapplot <- NULL
      vals$diffexBmName <- NULL
      myValues$dList <- NULL
    })
  })

  #Filter the selected features
  observeEvent(input$runFilterFeature, {
    filter <- rowData(vals$counts)[, input$filteredFeature] %in% input$filterFeatureChoices
    vals$counts <- vals$counts[filter, ]
    updateGeneNames()
    vals$diffexgenelist <- NULL
    vals$gsvaRes <- NULL
    vals$enrichRes <- NULL
    vals$visplotobject <- NULL
    vals$diffexheatmapplot <- NULL
    vals$diffexBmName <- NULL
    myValues$dList <- NULL
  })

  #disable the downloadSCE button if no object is loaded
  isAssayResult <- reactive(is.null(vals$counts))
  observe({
    if (isAssayResult()) {
      shinyjs::disable("downloadSCE")
    } else {
      shinyjs::enable("downloadSCE")
    }
  })

  output$downloadSCE <- downloadHandler(
    filename <- function() {
      paste("SCE-", Sys.Date(), ".rds", sep = "")
    },
    content <- function(file) {
      saveRDS(vals$counts, file)
    })

  output$assayList <- renderTable({
    req(vals$counts)
    if (!is.null(vals$counts) & length(names(assays(vals$counts))) > 0){
      data.table(assays = names(assays(vals$counts)))
    }
  })

  output$reducedDimsList <- renderTable({
    req(vals$counts)
    if (!is.null(vals$counts) & length(names(reducedDims(vals$counts))) > 0){
      data.table("Reduced Dimension" = names(reducedDims(vals$counts)))
    }
  })

  observeEvent(input$modifyAssay, {
    req(vals$counts)
    withBusyIndicatorServer("modifyAssay", {
      if (input$assayModifyAction == "log") {
        if (!(input$modifyAssaySelect %in% names(assays(vals$counts)))){
          shinyalert::shinyalert("Error!", "Assay does not exist!",
                                 type = "error")
        } else if (input$modifyAssayOutname == "") {
          shinyalert::shinyalert("Error!", "Invalid output name!",
                                 type = "error")
        } else if (input$modifyAssayOutname %in% names(assays(vals$counts))) {
          shinyalert::shinyalert("Error!", "Output name already exists! Delete to Rename.",
                                 type = "error")
        } else {
          assay(vals$counts, input$modifyAssayOutname) <- log2(assay(vals$counts, input$modifyAssaySelect) + 1)
          updateAssayInputs()
        }
      } else if (input$assayModifyAction == "cpm") {
        if (!(input$modifyAssaySelect %in% names(assays(vals$counts)))){
          shinyalert::shinyalert("Error!", "Assay does not exist!",
                                 type = "error")
        } else if (input$modifyAssayOutname == "") {
          shinyalert::shinyalert("Error!", "Invalid output name!",
                                 type = "error")
        } else if (input$modifyAssayOutname %in% names(assays(vals$counts))) {
          shinyalert::shinyalert("Error!", "Output name already exists! Delete to Rename.",
                                 type = "error")
        } else {
          assay(vals$counts, input$modifyAssayOutname) <- apply(assay(vals$counts, input$modifyAssaySelect), 2, function(x) { x / (sum(x) / 1000000) })
          updateAssayInputs()
        }
      } else if (input$assayModifyAction == "rename") {
        if (!(input$modifyAssaySelect %in% names(assays(vals$counts)))){
          shinyalert::shinyalert("Error!", "Assay does not exist!",
                                 type = "error")
        } else if (input$modifyAssayOutname == "") {
          shinyalert::shinyalert("Error!", "Invalid output name!",
                                 type = "error")
        } else if (input$modifyAssayOutname %in% names(assays(vals$counts))) {
          shinyalert::shinyalert("Error!", "Output name already exists! Delete to Rename.",
                                 type = "error")
        } else {
          assay(vals$counts, input$modifyAssayOutname) <- assay(vals$counts, input$modifyAssaySelect)
          assay(vals$counts, input$modifyAssaySelect) <- NULL
          updateAssayInputs()
        }
      } else if (input$assayModifyAction == "delete") {
        if (!(input$modifyAssaySelect %in% names(assays(vals$counts)))){
          shinyalert::shinyalert("Error!", "Assay does not exist!",
                                 type = "error")
        } else {
          assay(vals$counts, input$modifyAssaySelect) <- NULL
          updateAssayInputs()
        }
      } else {
        shinyalert::shinyalert("Error!", "Assay Modify Action Does Not Exist!",
                               type = "error")
      }
    })
  })

  observeEvent(input$delRedDim, {
    req(vals$counts)
    if (!(input$delRedDimType %in% names(reducedDims(vals$counts)))){
      shinyalert::shinyalert("Error!", "reducedDim does not exist!",
                             type = "error")
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

  #disable downloadcolData button if the data is not present
  isColDataResult <- reactive(is.null(vals$counts))
  observe({
    if (isColDataResult()) {
      shinyjs::disable("downloadcolData")
    } else {
      shinyjs::enable("downloadcolData")
    }
  })

  #download colData
  output$downloadcolData <- downloadHandler(
    filename = function() {
      paste("colData-", Sys.Date(), ".csv", sep = "")
    },
    content = function(file) {
      write.csv(data.frame(colData(vals$counts)), file)
    }
  )

  #upload and replace colData
  observeEvent(input$newAnnotFile, {
    req(input$newAnnotFile)
    req(vals$counts)
    indata <- read.csv(input$newAnnotFile$datapath, header = TRUE, row.names = 1)
    colData(vals$counts) <- DataFrame(indata)
    updateColDataNames()
  })

  #render UI for factor vs numeric
  output$annotModifyUI <- renderUI({
    if (!is.null(vals$counts)){
      if (input$annotModifyChoice != "none"){
        if (is.factor(colData(vals$counts)[, input$annotModifyChoice])){
          radioButtons("annotTypeSelect", "Field Type:", choices = c("factor", "numeric"), selected = "factor")
        } else {
          radioButtons("annotTypeSelect", "Field Type:", choices = c("factor", "numeric"), selected = "numeric")
        }
      }
    }
  })

  #update factor vs numeric for colData
  observeEvent(input$annotTypeSelect, {
    if (input$annotTypeSelect == "factor" & !is.factor(colData(vals$counts)[, input$annotModifyChoice])){
      colData(vals$counts)[, input$annotModifyChoice] <- as.factor(colData(vals$counts)[, input$annotModifyChoice])
    } else if (input$annotTypeSelect == "numeric" & is.factor(colData(vals$counts)[, input$annotModifyChoice])){
      f <- colData(vals$counts)[, input$annotModifyChoice]
      if (any(is.na(as.numeric(levels(f))[f]))){
        shinyalert::shinyalert("Error!", "Cannot convert to numeric.", type = "error")
      } else {
        colData(vals$counts)[, input$annotModifyChoice] <- as.numeric(levels(f))[f]
      }
    }
  })
  output$visOptions <- renderUI({
    if (!is.null(vals$counts)){
      if (input$visPlotMethod != "heatmap") {
        tagList(
          checkboxInput("visFWrap", "Plot genes individually?", value = TRUE)
        )
      } else {
        tagList(
          checkboxInput("visScaleHMap", "Scale expression values?", value = TRUE)
        )
      }
    }
  })
  observeEvent(input$plotvis, {
    if (is.null(vals$counts)){
      shinyalert::shinyalert("Error!", "Upload data first.", type = "error")
    } else{
      withBusyIndicatorServer("plotvis", {
        tryCatch({
          if (input$visCondn == "none"){
            incondition <- NULL
          } else {
            incondition <- input$visCondn
          }
          vals$visplotobject <- visPlot(inSCE = vals$counts,
                                        useAssay = input$visAssaySelect,
                                        method =  input$visPlotMethod,
                                        condition = incondition,
                                        glist = input$selectvisGenes,
                                        facetWrap = input$visFWrap,
                                        scaleHMap = input$visScaleHMap)
        },
        error = function(e){
          shinyalert::shinyalert("Error!", e$message, type = "error")
        })
      })
    }
  })

  output$visPlot <- renderPlot({
    req(vals$visplotobject)
    vals$visplotobject
  }, height = 600)

  #-----------------------------------------------------------------------------
  # Page 3: DR & Clustering
  #-----------------------------------------------------------------------------

  output$clusterPlot <- renderPlotly({
    if (is.null(vals$counts)){
      plotly::ggplotly(ggplot2::ggplot())
    } else{
      if (input$dimRedPlotMethod == "PCA"){
        pcadimname <- paste0("PCA", "_", input$dimRedAssaySelect)
        if (is.null(reducedDim(vals$counts, pcadimname))) {
          vals$counts <- getPCA(inSCE = vals$counts,
                                useAssay = input$dimRedAssaySelect,
                                reducedDimName = pcadimname)
          updateReddimInputs()
        }
        if (!is.null(reducedDim(vals$counts, pcadimname))){
          if (input$colorBy != "Gene Expression") {
            g <- singleCellTK::plotPCA(inSCE = vals$counts,
                                       colorBy = input$colorBy,
                                       shape = input$shapeBy, pcX = input$pcX,
                                       pcY = input$pcY,
                                       useAssay = input$dimRedAssaySelect,
                                       reducedDimName = pcadimname)
          } else if (input$colorGenes == ""){
            g <- singleCellTK::plotPCA(vals$counts, "No Color", "No Shape",
                                       input$pcX, input$pcY,
                                       useAssay = input$dimRedAssaySelect,
                                       reducedDimName = pcadimname)
          }
          plotly::ggplotly(g)
        } else {
          plotly::ggplotly(ggplot2::ggplot() + ggplot2::geom_point())
        }
      } else if (input$dimRedPlotMethod == "tSNE"){
        tsnedimname <- paste0("TSNE", "_", input$dimRedAssaySelect)
        if (is.null(reducedDim(vals$counts, tsnedimname))) {
          vals$counts <- getTSNE(inSCE = vals$counts,
                                 useAssay = input$dimRedAssaySelect,
                                 reducedDimName = tsnedimname)
          updateReddimInputs()
        }
        if (!is.null(reducedDim(vals$counts, tsnedimname))){
          if (input$colorBy != "Gene Expression") {
            g <- singleCellTK::plotTSNE(vals$counts, input$colorBy, input$shapeBy,
                                        useAssay = input$dimRedAssaySelect,
                                        reducedDimName = tsnedimname)
          } else if (input$colorGenes == ""){
            g <- singleCellTK::plotTSNE(vals$counts, "No Color", "No Shape",
                                        useAssay = input$dimRedAssaySelect,
                                        reducedDimName = tsnedimname)
          }
          plotly::ggplotly(g)
        } else {
          plotly::ggplotly(ggplot2::ggplot() + ggplot2::geom_point())
        }
      } else{
        plotly::ggplotly(ggplot2::ggplot() + ggplot2::theme_bw() +
                           ggplot2::theme(plot.background = ggplot2::element_rect(fill = "white")) +
                           ggplot2::theme(panel.border = ggplot2::element_rect(colour = "white")))
      }
    }
  })

  #TODO: this doesn't work with multiple pca dims
  output$pctable <- renderTable({
    if (is.null(vals$counts) | !(class(vals$counts) == "SCtkExperiment")){
    } else{
      if (input$dimRedPlotMethod == "PCA") {
        if (nrow(pcaVariances(vals$counts)) == ncol(vals$counts)){
          data.frame(PC = paste("PC", seq_len(ncol(vals$counts)), sep = ""),
                     Variances = pcaVariances(vals$counts)$percentVar * 100)[1:10, ]
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
          vals$counts <- getPCA(inSCE = vals$counts,
                                useAssay = input$dimRedAssaySelect,
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
                ggplot2::ggplot() + ggplot2::theme_bw() +
                  ggplot2::theme(plot.background = ggplot2::element_rect(fill = "white")) +
                  ggplot2::theme(panel.border = ggplot2::element_rect(colour = "white"))
              } else {
                biomarkers <- data.frame(eval(parse(text = paste("rowData(vals$counts)[,'", input$colorGenesBiomarker, "']", sep = ""))))
                rownames(biomarkers) <- rowData(vals$counts)[, "Gene"]
                biomarkers <- rownames(subset(biomarkers, biomarkers[, 1] == 1))
                g <- plotBiomarker(vals$counts, biomarkers, input$colorBinary,
                                   "PCA", input$shapeBy, input$pcX, input$pcY,
                                   useAssay = input$dimRedAssaySelect,
                                   reducedDimName = pcadimname)
                g
              }
            } else if (input$colorGeneBy == "Manual Input") {
              if (is.null(input$colorGenes)){
                ggplot2::ggplot() + ggplot2::theme_bw() +
                  ggplot2::theme(plot.background = ggplot2::element_rect(fill = "white")) +
                  ggplot2::theme(panel.border = ggplot2::element_rect(colour = "white"))
              } else {
                g <- plotBiomarker(vals$counts, input$colorGenes,
                                   input$colorBinary, "PCA", input$shapeBy,
                                   input$pcX, input$pcY,
                                   useAssay = input$dimRedAssaySelect,
                                   reducedDimName = pcadimname)
                g
              }
            }
          }
        }
      } else if (input$dimRedPlotMethod == "tSNE"){
        tsnedimname <- paste0("TSNE", "_", input$dimRedAssaySelect)
        if (is.null(reducedDim(vals$counts, tsnedimname))){
          vals$counts <- getTSNE(inSCE = vals$counts,
                                 useAssay = input$dimRedAssaySelect,
                                 reducedDimName = tsnedimname)
          updateReddimInputs()
        }
        if (!is.null(reducedDim(vals$counts, tsnedimname))){
          if (input$colorBy == "Gene Expression") {
            if (input$colorGeneBy == "Biomarker (from DE tab)"){
              if (input$colorGenesBiomarker == ""){
                ggplot2::ggplot() +
                  ggplot2::theme_bw() +
                  ggplot2::theme(plot.background = ggplot2::element_rect(fill = "white")) +
                  ggplot2::theme(panel.border = ggplot2::element_rect(colour = "white"))
              } else {
                biomarkers <- data.frame(eval(parse(text = paste("rowData(vals$counts)[,'", input$colorGenesBiomarker, "']", sep = ""))))
                rownames(biomarkers) <- rowData(vals$counts)[, "Gene"]
                biomarkers <- rownames(subset(biomarkers, biomarkers[, 1] == 1))
                g <- plotBiomarker(vals$counts, biomarkers, input$colorBinary,
                                   "tSNE", input$shapeBy,
                                   useAssay = input$dimRedAssaySelect,
                                   reducedDimName = tsnedimname)
                g
              }
            } else if (input$colorGeneBy == "Manual Input") {
              if (is.null(input$colorGenes)){
                ggplot2::ggplot() + ggplot2::theme_bw() +
                  ggplot2::theme(plot.background = ggplot2::element_rect(fill = "white")) +
                  ggplot2::theme(panel.border = ggplot2::element_rect(colour = "white"))
              } else {
                g <- plotBiomarker(vals$counts, input$colorGenes,
                                   input$colorBinary, "tSNE", input$shapeBy,
                                   useAssay = input$dimRedAssaySelect,
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
      shinyalert::shinyalert("Error!", "Upload data first.", type = "error")
    } else {
      if (input$dimRedPlotMethod == "Dendrogram" & paste0("PCA", "_", input$dimRedAssaySelect) %in% names(reducedDims(vals$counts))){
        data <- getClusterInputData(inSCE = vals$counts,
                                    inputData = "PCA Components",
                                    useAssay = input$dimRedAssaySelect,
                                    reducedDimName = paste0("PCA", "_", input$dimRedAssaySelect))
        d <- stats::dist(data)
        h <- stats::hclust(d, input$dendroDistanceMetric)
        if (input$clusteringAlgorithmD == "Phylogenetic Tree") {
          g <- ggtree::ggtree(as.phylo(h), layout = "circular", open.angle = 360) + ggtree::geom_tiplab2(size = 2)
        } else if (input$clusteringAlgorithmD == "Hierarchical") {
          g <- ggtree::ggtree(as.phylo(h)) + ggtree::theme_tree2() + ggtree::geom_tiplab(size = 2)
        } else {
          stop("Input clustering algorithm not found ", input$clusteringAlgorithmD)
        }
        g
      }
    }
  }, height = 600)

  observeEvent(input$clusterData, {
    if (is.null(vals$counts)){
      shinyalert::shinyalert("Error!", "Upload data first.", type = "error")
    } else if (input$clusterName == "") {
      shinyalert::shinyalert("Error!", "Cluster name required.", type = "error")
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
                                      useAssay = input$dimRedAssaySelect,
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
          data <- getClusterInputData(inSCE = vals$counts,
                                      inputData = input$selectClusterInputData,
                                      useAssay = input$dimRedAssaySelect,
                                      reducedDimName = currdimname)
          coutput <- cluster::clara(data, input$Cnumber)
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
      if (input$dimRedPlotMethod == "PCA"){
        pcadimname <- paste0("PCA", "_", input$dimRedAssaySelect)
        if (!is.null(reducedDim(vals$counts, pcadimname))) {
          currPcs <- colnames(reducedDim(vals$counts, "PCA_counts"))
          updateSelectInput(session, "pcX", choices = currPcs,
                            selected = currPcs[1])
          updateSelectInput(session, "pcY", choices = currPcs,
                            selected = currPcs[2])
        }
      }
    }
  })

  observeEvent(input$reRunTSNE, {
    if (is.null(vals$original)){
      shinyalert::shinyalert("Error!", "Upload data first.", type = "error")
    }
    else{
      withBusyIndicatorServer("reRunTSNE", {
        vals$counts <- getTSNE(inSCE = vals$counts,
                               useAssay = input$dimRedAssaySelect,
                               reducedDimName = paste0("TSNE", "_",
                                                       input$dimRedAssaySelect))
        updateReddimInputs()
      })
    }
  })

  observeEvent(input$reRunPCA, {
    if (is.null(vals$original)){
      shinyalert::shinyalert("Error!", "Upload data first.", type = "error")
    }
    else{
      withBusyIndicatorServer("reRunPCA", {
        vals$counts <- getPCA(inSCE = vals$counts,
                              useAssay = input$dimRedAssaySelect,
                              reducedDimName = paste0("PCA", "_",
                                                      input$dimRedAssaySelect))
        updateReddimInputs()
      })
    }
  })

  #-----------------------------------------------------------------------------
  # Page 4: Batch Correction
  #-----------------------------------------------------------------------------

  output$selectCombatRefBatchUI <- renderUI({
    if (!is.null(vals$counts)){
      if (input$combatRef){
        selectInput("combatRefBatch", "Choose Reference Batch:",
                    unique(sort(colData(vals$counts)[, input$combatBatchVar])))
      }
    }
  })

  observeEvent(input$combatRun, {
    if (is.null(vals$counts)){
      shinyalert::shinyalert("Error!", "Upload data first.", type = "error")
    }
    else{
      withBusyIndicatorServer("combatRun", {
        if (input$batchMethod == "ComBat"){
          #check for zeros
          if (any(rowSums(assay(vals$counts, input$combatAssay)) == 0)){
            shinyalert::shinyalert("Error!", "Rows with a sum of zero found. Filter data to continue.", type = "error")
          } else {
            saveassayname <- gsub(" ", "_", input$combatSaveAssay)
            if (input$combatRef){
              assay(vals$counts, saveassayname) <-
                ComBatSCE(inSCE = vals$counts, batch = input$combatBatchVar,
                          useAssay = input$combatAssay,
                          par.prior = input$combatParametric,
                          covariates = input$combatConditionVar,
                          mean.only = input$combatMeanOnly,
                          ref.batch = input$combatRefBatch)
            } else {
              assay(vals$counts, saveassayname) <-
                ComBatSCE(inSCE = vals$counts, batch = input$combatBatchVar,
                          useAssay = input$combatAssay,
                          par.prior = input$combatParametric,
                          covariates = input$combatConditionVar,
                          mean.only = input$combatMeanOnly)
            }
            updateAssayInputs()
            vals$combatstatus <- "ComBat Complete"
          }
        } else {
          shinyalert::shinyalert("Error!", "Unsupported Batch Correction Method", type = "error")
        }
      })
    }
  })

  output$combatStatus <- renderUI({
    h2(vals$combatstatus)
  })

  output$combatBoxplot <- renderPlot({
    if (!is.null(vals$counts) &
        !is.null(input$batchVarPlot) &
        !is.null(input$conditionVarPlot) &
        input$batchVarPlot != "none" &
        input$conditionVarPlot != "none" &
        input$batchVarPlot != input$conditionVarPlot){
      plotBatchVariance(inSCE = vals$counts,
                        useAssay = input$combatAssay,
                        batch = input$batchVarPlot,
                        condition = input$conditionVarPlot)
    }
  }, height = 600)

  #-----------------------------------------------------------------------------
  # Page 5.1: Differential Expression
  #-----------------------------------------------------------------------------
  shinyjs::onclick("Diffex_hideAllSections", allSections(
    "hide", c(paste("de", 1:7, sep = ""))), add = TRUE)
  shinyjs::onclick("Diffex_showAllSections", allSections(
    "show", c(paste("de", 1:7, sep = ""))), add = TRUE)
  shinyjs::onclick("diffex1",
                   shinyjs::toggle(id = "de1",
                                   anim = TRUE), add = TRUE)
  shinyjs::onclick("diffex2",
                   shinyjs::toggle(id = "de2",
                                   anim = TRUE), add = TRUE)
  shinyjs::onclick("diffex3",
                   shinyjs::toggle(id = "de3",
                                   anim = TRUE), add = TRUE)
  shinyjs::onclick("diffex4",
                   shinyjs::toggle(id = "de4",
                                   anim = TRUE), add = TRUE)
  shinyjs::onclick("diffex5",
                   shinyjs::toggle(id = "de5",
                                   anim = TRUE), add = TRUE)
  shinyjs::onclick("diffex6",
                   shinyjs::toggle(id = "de6",
                                   anim = TRUE), add = TRUE)
  shinyjs::onclick("diffex7",
                   shinyjs::toggle(id = "de7",
                                   anim = TRUE), add = TRUE)
  shinyjs::addClass(id = "diffex1", class = "btn-block")
  shinyjs::addClass(id = "diffex2", class = "btn-block")
  shinyjs::addClass(id = "diffex3", class = "btn-block")
  shinyjs::addClass(id = "diffex4", class = "btn-block")
  shinyjs::addClass(id = "diffex5", class = "btn-block")
  shinyjs::addClass(id = "diffex6", class = "btn-block")
  shinyjs::addClass(id = "diffex7", class = "btn-block")

  output$selectDiffexConditionUI <- renderUI({
    if (!is.null(vals$counts)){
      if (input$selectDiffex == "ANOVA") {
        tagList(
          selectInput("selectDiffexCondition", "Select Condition(s):",
                      colnames(colData(vals$counts)), multiple = TRUE)
        )
      } else {
        tagList(
          selectInput("selectDiffexCondition",
                      "Select Condition:",
                      colnames(colData(vals$counts))),
          selectInput("selectDiffexCovariates",
                      "Select Additional Covariates:",
                      colnames(colData(vals$counts)), multiple = TRUE)
        )
      }
    }
  })

  #For conditions with more than two factors, select the factor of interest
  output$selectDiffexConditionLevelUI <- renderUI({
    req(vals$counts)
    if (length(colnames(colData(vals$counts))) > 0){
      if (length(unique(colData(vals$counts)[, input$selectDiffexCondition])) > 2 & input$selectDiffex == "DESeq2"){
        tagList(
          radioButtons("selectDiffexConditionMethod", "Select Analysis Method:",
                       choiceNames = c("Biomarker (1 vs all)", "Factor of Interest vs. Control Factor",
                                       "Entire Factor (Full/Reduced)"),
                       choiceValues = c("biomarker", "contrast", "fullreduced")
          ),
          conditionalPanel(
            condition = "input.selectDiffexConditionMethod != 'fullreduced'",
            selectInput("selectDiffexConditionOfInterest",
                        "Select Factor of Interest",
                        unique(sort(colData(vals$counts)[, input$selectDiffexCondition])))
          ),
          conditionalPanel(
            condition = "input.selectDiffexConditionMethod == 'contrast' && input.selectDiffex == 'DESeq2'",
            selectInput("selectDiffexControlCondition",
                        "Select Control Factor",
                        unique(sort(colData(vals$counts)[, input$selectDiffexCondition])))
          )
        )
      } else if (length(unique(colData(vals$counts)[, input$selectDiffexCondition])) > 2 & input$selectDiffex == "limma") {
        tagList(
          radioButtons("selectDiffexConditionMethod", "Select Analysis Method:",
                       choiceNames = c("Biomarker (1 vs all)", "Factor of Interest",
                                       "Entire Factor"),
                       choiceValues = c("biomarker", "coef", "allcoef")
          ),
          conditionalPanel(
            condition = "input.selectDiffexConditionMethod != 'allcoef'",
            selectInput("selectDiffexConditionOfInterest",
                        "Select Factor of Interest",
                        unique(sort(colData(vals$counts)[, input$selectDiffexCondition])))
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
      shinyalert::shinyalert("Error!", "Upload data first.", type = "error")
    }
    else{
      withBusyIndicatorServer("runDiffex", {
        vals$diffexheatmapplot <- NULL
        #run diffex to get gene list and pvalues
        if (input$selectDiffex == "ANOVA"){
          useCovariates <- input$anovaCovariates
        } else {
          useCovariates <- input$selectDiffexCovariates
        }
        vals$diffexgenelist <- scDiffEx(inSCE = vals$counts,
                                        useAssay = input$diffexAssay,
                                        condition = input$selectDiffexCondition,
                                        covariates = useCovariates,
                                        significance = input$selectPval,
                                        ntop = nrow(vals$counts),
                                        usesig = FALSE,
                                        diffexmethod = input$selectDiffex,
                                        levelofinterest = input$selectDiffexConditionOfInterest,
                                        analysisType = input$selectDiffexConditionMethod,
                                        controlLevel = input$selectDiffexControlCondition,
                                        adjust = input$selectCorrection)
      })
    }
  })

  output$colorBarConditionUI <- renderUI({
    if (is.null(vals$counts)){
      selectInput("colorBarCondition", "Select Condition", NULL)
    } else {
      selectInput("colorBarCondition", "Select Condition",
                  colnames(colData(vals$counts)), multiple = TRUE)
    }
  })

  annotationColors <- reactiveValues(cols = list())

  output$heatmapSampleAnnotations <- renderUI({
    if (!is.null(input$colorBarCondition)) {
      if (!is.null(vals$counts) & length(input$colorBarCondition) > 0){
        if (all(input$colorBarCondition %in% colnames(colData(vals$counts)))) {
          h <- input$colorBarCondition
          L <- lapply(seq_along(h), function(i) colourGroupInput(paste0("colorGroup", i)))
          annotationColors$cols <- lapply(
            seq_along(h),
            function(i) {
              callModule(colourGroup, paste0("colorGroup", i), heading = h[i],
                         options = unique(unlist(colData(vals$counts)[, h[i]])))
            }
          )
          return(L)
        }
      }
    }
  })

  output$diffexNgenes <- renderUI({
    req(vals$diffexgenelist)
    HTML(paste(em("Max genes: "), nrow(vals$diffexgenelist), sep = ""))
  })

  output$logFCDiffexRange <- renderUI({
    req(vals$diffexgenelist)
    if (input$selectDiffex != 'ANOVA') {
      logFCIndex <- which(grepl("*log*", colnames(vals$diffexgenelist)))
      if (length(logFCIndex) == 0) {
        #for DESeq2 with more than 1 covariate, choose the first column
        logFCIndex <- 1
      }
      minlogFC <- paste(em("Min logFC : "), round(min(na.omit(vals$diffexgenelist[, logFCIndex])), digits = 6))
      maxlogFC <- paste(em("Max logFC : "), round(max(na.omit(vals$diffexgenelist[, logFCIndex])), digits = 6))
      HTML(paste(minlogFC, maxlogFC, sep = '<br/>'))
    }
  })

  #Plot the differential expression results
  observeEvent(input$runPlotDiffex, {
    req(vals$diffexgenelist)
    withBusyIndicatorServer("runPlotDiffex", {
      tryCatch ({
        #logFC or abs(logFC)
        if (input$applyAbslogFCDiffex == TRUE) {
          absLogFCDiffex <- abs(input$selectlogFCDiffex)
        } else {
          absLogFCDiffex <- input$selectlogFCDiffex
        }
        #for convenience, index logFC and p-val columns for all the methods
        pvalIndex <- which(grepl("*padj*", colnames(vals$diffexgenelist)))
        logFCIndex <- which(grepl("*log*", colnames(vals$diffexgenelist)))
        if (input$selectNGenes > nrow(vals$diffexgenelist)) {
          stop("Max value exceeded for Input.")
        }
        #p-Val and logFC cutoff
        if (input$applyCutoff == TRUE & input$applylogFCCutoff == TRUE) {
          if (input$selectDiffex == 'ANOVA') {
            stop("logFC is not applicable for ANOVA")
          } else {
            if (min(na.omit(vals$diffexgenelist[, pvalIndex])) > input$selectPval) {
              diffexFilterRes <- vals$diffexgenelist
              stop("the min/least p-value in the results is greater than the selected p-val range")
            } else if (min(na.omit(vals$diffexgenelist[, logFCIndex])) > absLogFCDiffex) {
              diffexFilterRes <- vals$diffexgenelist
              stop("the min/least logFC in the results is greater than the selected logFC range")
            } else {
              diffexFilterRes <-  vals$diffexgenelist[(vals$diffexgenelist[, pvalIndex] <= input$selectPval &
                                                         vals$diffexgenelist[, logFCIndex] <= absLogFCDiffex), ]
            }
          }
        }
        #p-Val cutoff
        else if (input$applyCutoff == TRUE) {
          if (min(na.omit(vals$diffexgenelist[, pvalIndex])) > input$selectPval) {
            diffexFilterRes <- vals$diffexgenelist
            stop("the min/least p-value in the results is greater than the selected p-val range")
          } else {
            diffexFilterRes <-  vals$diffexgenelist[(vals$diffexgenelist[, pvalIndex] <= input$selectPval), ]
          }
        }
        #logFC cutoff
        else if (input$applylogFCCutoff == TRUE) {
          if (input$selectDiffex == 'ANOVA') {
            stop("logFC is not applicable for ANOVA")
          } else  {
            if (min(na.omit(vals$diffexgenelist[, logFCIndex])) > absLogFCDiffex) {
              diffexFilterRes <- vals$diffexgenelist
              stop("the min/least logFC in the results is greater than the selected logFC range")
            } else {
              diffexFilterRes <-  vals$diffexgenelist[(vals$diffexgenelist[, logFCIndex] <= absLogFCDiffex), ]
            }
          }
        } else {
          diffexFilterRes <- vals$diffexgenelist
        }
        if (is.null(diffexFilterRes)){
          diffexFilterRes <- vals$diffexgenelist
        }
        rowLengthFiltered <- nrow(diffexFilterRes)
        if (rowLengthFiltered == 0) {
          stop("You've got 0 genes after filtering.. adjust your filters accordingly")
        }
        if (rowLengthFiltered < input$selectNGenes) {
          diffexFilterRes <- diffexFilterRes[seq_len(rowLengthFiltered), ]
        } else {
          diffexFilterRes <- diffexFilterRes[seq_len(input$selectNGenes), ]
        }
        #run plotDiffex
        if (!is.null(vals$diffexgenelist)){
          if (input$displayHeatmapColorBar){
            if (is.null(input$colorBarCondition)){
              colors <- NULL
            } else {
              colors <- lapply(annotationColors$cols, function(col) col())
              names(colors) <- input$colorBarCondition
              if (is.null(colors[[length(colors)]][[1]])){
                colors <- NULL
              }
            }
          } else {
            colors <- NULL
          }
          vals$diffexheatmapplot <- plotDiffEx(inSCE = vals$counts,
                                               useAssay = input$diffexAssay,
                                               condition = input$colorBarCondition,
                                               geneList = rownames(diffexFilterRes),
                                               clusterRow = input$clusterRows,
                                               clusterCol = input$clusterColumns,
                                               displayRowLabels = input$displayHeatmapRowLabels,
                                               displayColumnLabels = input$displayHeatmapColumnLabels,
                                               displayRowDendrograms = input$displayHeatmapRowDendrograms,
                                               displayColumnDendrograms = input$displayHeatmapColumnDendrograms,
                                               annotationColors = colors,
                                               scaleExpression = input$applyScaleDiffex,
                                               columnTitle = input$heatmapColumnsTitle)
        }
      }, error = function(e){
        shinyalert::shinyalert("Error!", e$message, type = "error")
      })
    })
  })

  output$diffPlot <- renderPlot({
    req(vals$diffexheatmapplot)
    ComplexHeatmap::draw(vals$diffexheatmapplot)
  }, height = 600)

  #Create the differential expression results table
  output$diffextable <- DT::renderDataTable({
    if (!is.null(vals$diffexgenelist)){
      temptable <- cbind(rownames(vals$diffexgenelist), data.frame(vals$diffexgenelist))
      colnames(temptable)[1] <- "Gene"
      temptable
    }
  }, rownames = FALSE)

  #disable downloadGeneList button if the result is not null
  isDiffExResult <- reactive(is.null(vals$diffexgenelist))
  observe({
    if (isDiffExResult()) {
      shinyjs::disable("downloadGeneList")
    } else {
      shinyjs::enable("downloadGeneList")
    }
  })

  #create custom name for the results
  customName <- reactive(paste(input$selectDiffexCondition, input$selectDiffex, sep = "_"))
  # Download the differential expression results table
  output$downloadGeneList <- downloadHandler(
    filename = function() {
      paste(customName(), Sys.Date(), ".csv", sep = "")
    },
    content = function(file) {
      results <- vals$diffexgenelist
      colnames(results) <- paste(customName(), colnames(results), sep = "_")
      utils::write.csv(results, file)
    }
  )

  #reactive list to store names of results given by the user.
  myValues <- reactiveValues(
    index = 0
  )

  #save results wrt to custom name
  observeEvent(input$saveResults, {
    if (input$ResultsName == ""){
      shinyalert::shinyalert("Error!", "Specify name of the results.", type = "error")
    } else {
      withBusyIndicatorServer("saveResults", {
        ResultsName <- gsub(" ", "_", input$ResultsName)
        if (length(myValues$dList) >= 1) {
          if (anyDuplicated(myValues$dList)) {
            shinyalert::shinyalert("Error", "name already exists. Please use a unique result name",
                                   type = "error")
          } else {
            myValues$index <- myValues$index + 1
            #myValues$dList[myValues$index] <- isolate(input$ResultsName)
            myValues$dList[myValues$index] <- ResultsName
          }
        } else {
          myValues$index <- myValues$index + 1
          #myValues$dList[myValues$index] <- isolate(input$ResultsName)
          myValues$dList[myValues$index] <- ResultsName
        }
        vals$counts <- saveDiffExResults(inSCE = vals$counts,
                                         diffex = vals$diffexgenelist,
                                         name = input$ResultsName,
                                         method = input$selectDiffex)
      })
    }
  })

  #dynamically create a list of names of the results
  output$savedRes <- renderUI({
    if (!is.null(vals$counts)) {
      if (is.null(myValues$dList)) {
        savedObjResults <- gsub("_padj$", "", colnames(rowData(vals$counts))[grepl("_padj$", colnames(rowData(vals$counts)))])
        #myValues$index <- myValues$index + 1
        myValues$index <- length(savedObjResults)
        myValues$dList[seq_len(myValues$index)] <- savedObjResults[seq_len(myValues$index)]
      }
      selectizeInput("savedDiffExResults", "Select available results",
                     choices = myValues$dList)
    }
  })

  output$saveDiffResultsNote <- renderUI({
    req(vals$diffexgenelist)
    HTML(paste(em("Note: Use a unique name to save results each time.")))
  })

  #load specific result according to users' input
  observeEvent(input$loadResults, {
    if (!is.null(input$savedDiffExResults)) {
      df <- data.frame(rowData(vals$counts))
      #extract all columns except biomarker columns by checking for unique values == 2
      listColNames <- names(which(apply(df, 2, function(a) length(unique(a)) == 2) == FALSE))
      #arrange your df according to these extracted names
      df <- df[, listColNames]
      #select columns based on the saved result selected
      selected_cols <- sub(listColNames[1], "_", input$savedDiffExResults, fixed = TRUE)
      #now choose the unique result
      str_match <- unique(sub("^_[^_]+_$", "", selected_cols, fixed = TRUE))
      df <- df[, which(grepl(paste0(str_match, "_"), listColNames) == TRUE)]
      diffexRow <- rownames(df)[seq_len(nrow(df))]
      rownames(df) <- rownames(vals$counts)[as.integer(diffexRow)]
      colnames(df) <- gsub(paste0(str_match, "_"), "", colnames(df))
      vals$diffexgenelist <- df
      vals$diffexheatmapplot <- NULL
    }
  })

  output$BioNgenes <- renderUI({
    req(vals$diffexgenelist)
    HTML(paste(em("Max genes: "), nrow(vals$diffexgenelist), sep = ""))
  })

  output$logFCBioRange <- renderUI({
    req(vals$diffexgenelist)
    if (input$selectDiffex != 'ANOVA') {
      logFCIndex <- which(grepl("*log*", colnames(vals$diffexgenelist)))
      if (length(logFCIndex)) {
        #for DESeq2 with more than 1 covariate, choose the first column
        logFCIndex <- 1
      }
      minlogFC <- paste(em("Min logFC : "), round(min(na.omit(vals$diffexgenelist[, logFCIndex])), digits = 6))
      maxlogFC <- paste(em("Max logFC : "), round(max(na.omit(vals$diffexgenelist[, logFCIndex])), digits = 6))
      HTML(paste(minlogFC, maxlogFC, sep = '<br/>'))
    }
  })

  #save biomarker in rowData() wrt name and conditions.
  observeEvent(input$saveBiomarker, {
    if (input$biomarkerName == ""){
      shinyalert::shinyalert("Error!", "Specify biomarker name.", type = "error")
    } else {
      withBusyIndicatorServer("saveBiomarker", {
        req(vals$diffexgenelist)
        biomarkerName <- gsub(" ", "_", input$biomarkerName)
        if (anyDuplicated(biomarkerName)) {
          shinyalert::shinyalert("Error", "name already exists. Please use a unique result name",
                                 type = "error")
        }
        if (input$applyAbslogFC == TRUE) {
          absLogFC <- abs(input$selectlogFC)
        } else {
          absLogFC <- input$selectlogFC
        }
        if (input$selectBioNGenes > nrow(vals$diffexgenelist)) {
          stop("Max value exceeded for Input.")
        }
        if (input$applyBioCutoff1 == TRUE & input$applyBioCutoff2 == TRUE) {
          vals$counts <- saveBiomarkerRes(inSCE = vals$counts,
                                          diffex = vals$diffexgenelist,
                                          biomarkerName = biomarkerName,
                                          method = input$selectDiffex,
                                          ntop = input$selectBioNGenes,
                                          logFC = absLogFC,
                                          pVal = input$selectAdjPVal)
        } else if (input$applyBioCutoff1 == TRUE) {
          vals$counts <- saveBiomarkerRes(inSCE = vals$counts,
                                          diffex = vals$diffexgenelist,
                                          biomarkerName = biomarkerName,
                                          method = input$selectDiffex,
                                          ntop = input$selectBioNGenes,
                                          logFC = NULL,
                                          pVal = input$selectAdjPVal)
        } else if (input$applyBioCutoff2 == TRUE) {
          vals$counts <- saveBiomarkerRes(inSCE = vals$counts,
                                          diffex = vals$diffexgenelist,
                                          biomarkerName = biomarkerName,
                                          method = input$selectDiffex,
                                          ntop = input$selectBioNGenes,
                                          logFC = absLogFC,
                                          pVal = NULL)
        } else {
          vals$counts <- saveBiomarkerRes(inSCE = vals$counts,
                                          diffex = vals$diffexgenelist,
                                          biomarkerName = input$biomarkerName,
                                          method = input$selectDiffex,
                                          ntop = input$selectBioNGenes,
                                          logFC = NULL,
                                          pVal = NULL)
        }
        vals$diffexBmName <- TRUE
      })
    }
  })

  observe({
    output$bioMarkerNote <- renderUI({
      req(vals$counts)
      req(isolate(input$biomarkerName))
      if (vals$diffexBmName) {
        isolate({
          biomarkerName <- gsub(" ", "_", input$biomarkerName)
          countBioGenes <- count(rowData(vals$counts)[, biomarkerName] == 1)
          HTML(paste("Saved ", countBioGenes, " genes after applying the selected filter(s)", sep = ""))
        })
      }
    })
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
      shinyalert::shinyalert("Error!", "Upload data first.", type = "error")
    } else {
      withBusyIndicatorServer("runDEhurdle", {
        #run diffex to get gene list and pvalues
        vals$mastgenelist <- MAST(inSCE = vals$counts,
                                  useAssay = input$mastAssay,
                                  condition = input$hurdlecondition,
                                  interest.level = input$hurdleconditionofinterest,
                                  freqExpressed = input$hurdlethresh,
                                  fcThreshold = input$FCthreshold,
                                  p.value = input$hurdlepvalue,
                                  useThresh = input$useAdaptThresh)
      })
    }
  })

  observeEvent(input$runThreshPlot, {
    if (is.null(vals$counts)){
      shinyalert::shinyalert("Error!", "Upload data first.", type = "error")
    }
    else{
      withBusyIndicatorServer("runThreshPlot", {
        output$threshplot <- renderPlot({
          vals$thres <- thresholdGenes(inSCE = vals$counts,
                                       useAssay = input$mastAssay)
          par(mfrow = c(5, 4))
          plot(vals$thres)
          par(mfrow = c(1, 1))
        }, height = 600)
      })
    }
  })

  output$hurdleviolin <- renderPlot({
    if (!(is.null(vals$mastgenelist))){
      MASTviolin(inSCE = vals$counts, useAssay = input$mastAssay,
                 fcHurdleSig = vals$mastgenelist,
                 condition = input$hurdlecondition,
                 threshP = input$useAdaptThresh)
    }
  }, height = 600)

  output$hurdlelm <- renderPlot({
    if (!(is.null(vals$mastgenelist))){
      MASTregression(inSCE = vals$counts, useAssay = input$mastAssay,
                     fcHurdleSig = vals$mastgenelist,
                     condition = input$hurdlecondition,
                     threshP = input$useAdaptThresh)
    }
  }, height = 600)

  output$hurdleHeatmap <- renderPlot({
    if (!(is.null(vals$mastgenelist))){
      draw(plotDiffEx(vals$counts, useAssay = input$mastAssay,
                      condition = input$hurdlecondition,
                      geneList = vals$mastgenelist$Gene,
                      annotationColors = "auto", columnTitle = "MAST"))
    }
  }, height = 600)

  #Create the MAST results table
  output$mastresults <- DT::renderDataTable({
    if (!is.null(vals$mastgenelist)){
      vals$mastgenelist
    }
  })

  #disable mast dowload button if the mastgenelist data is null
  isMastGeneListResult <- reactive(is.null(vals$mastgenelist))
  observe({
    if (isMastGeneListResult()) {
      shinyjs::disable("downloadHurdleResult")
    } else {
      shinyjs::enable("downloadHurdleResult")
    }
  })

  #download mast results
  output$downloadHurdleResult <- downloadHandler(
    filename = function() {
      paste("mast_results-", Sys.Date(), ".csv", sep = "")
    },
    content = function(file) {
      utils::write.csv(vals$mastgenelist, file)
    }
  )

  #-----------------------------------------------------------------------------
  # Page 6: Pathway Activity Analysis
  #-----------------------------------------------------------------------------

  output$selectPathwayGeneLists <- renderUI({
    if (input$genelistSource == "Manual Input"){
      if (!is.null(vals$counts)){
        #fn to check if each column is 1 and 0 only
        biomarkercols <- names(which(apply(rowData(vals$counts), 2, function(a) length(unique(a)) == 2) == TRUE))
        selectizeInput("pathwayGeneLists", "Select Gene List(s):",
                       biomarkercols, multiple = TRUE)
      } else {
        h4("Note: upload data.")
      }
    } else {
      selectInput("pathwayGeneLists", "Select Gene List(s):",
                  c("ALL", names(c2BroadSets)), multiple = TRUE)
    }
  })

  output$selectNumTopPaths <- renderUI({
    if (!is.null(input$pathwayGeneLists)) {
      if ("ALL" %in% input$pathwayGeneLists & input$genelistSource == "MSigDB c2 (Human, Entrez ID only)"){
        sliderInput("pickNtopPaths", "Number of top pathways:", min = 5,
                    max = length(c2BroadSets), value = 25, step = 5)
      }
    }
  })

  observeEvent(input$pathwayRun, {
    if (is.null(vals$counts)){
      shinyalert::shinyalert("Error!", "Upload data first.", type = "error")
    } else {
      withBusyIndicatorServer("pathwayRun", {
        vals$gsvaRes <- gsvaSCE(inSCE = vals$counts,
                                useAssay = input$pathwayAssay,
                                pathwaySource = input$genelistSource,
                                pathwayNames = input$pathwayGeneLists)
      })
    }
  })

  observe({
    if (length(input$pathwayPlotVar) == 1 & !(is.null(vals$gsvaRes))){
      fit <- limma::lmFit(vals$gsvaRes, stats::model.matrix(~factor(colData(vals$counts)[, input$pathwayPlotVar])))
      fit <- limma::eBayes(fit)
      toptableres <- limma::topTable(fit, number = nrow(vals$gsvaRes))
      temptable <- cbind(rownames(toptableres), toptableres)
      rownames(temptable) <- NULL
      colnames(temptable)[1] <- "Pathway"
      vals$gsvaLimma <- temptable
    } else {
      vals$gsvaLimma <- NULL
    }
  })

  output$pathwaytable <- DT::renderDataTable({
    if (!is.null(vals$gsvaLimma)){
      if (!is.null(input$pathwayGeneLists) & "ALL" %in% input$pathwayGeneLists & input$genelistSource == "MSigDB c2 (Human, Entrez ID only)"){
        vals$gsvaLimma[1:min(input$pickNtopPaths, nrow(vals$gsvaLimma)), , drop = FALSE]
      } else {
        vals$gsvaLimma
      }
    }
  }, options = list(scrollX = TRUE, pageLength = 30))

  output$pathwayPlot <- renderPlot({
    if (!(is.null(vals$gsvaRes))){
      if (input$genelistSource == "MSigDB c2 (Human, Entrez ID only)" & "ALL" %in% input$pathwayGeneLists & !(is.null(vals$gsvaLimma))){
        tempgsvares <- vals$gsvaRes[as.character(vals$gsvaLimma$Pathway[1:min(input$pickNtopPaths, nrow(vals$gsvaLimma))]), , drop = FALSE]
      } else if (input$genelistSource == "MSigDB c2 (Human, Entrez ID only)" & !("ALL" %in% input$pathwayGeneLists)) {
        tempgsvares <- vals$gsvaRes
      } else {
        tempgsvares <- vals$gsvaRes[1:input$pickNtopPaths, , drop = FALSE]
      }
      if (input$pathwayOutPlot == "Violin" & length(input$pathwayPlotVar) > 0){
        tempgsvares <- tempgsvares[1:min(49, input$pickNtopPaths, nrow(tempgsvares)), , drop = FALSE]
        gsvaPlot(inSCE = vals$counts,
                 gsvaData = tempgsvares,
                 plotType = input$pathwayOutPlot,
                 condition = input$pathwayPlotVar)
      } else if (input$pathwayOutPlot == "Heatmap"){
        gsvaPlot(inSCE = vals$counts,
                 gsvaData = tempgsvares,
                 plotType = input$pathwayOutPlot,
                 condition = input$pathwayPlotVar)
      }
    }
  })

  #save pathawy activity results in the colData
  observeEvent(input$savePathway, {
    if (!(is.null(vals$gsvaRes))){
      if (all(colnames(vals$counts) == colnames(vals$gsvaRes))){
        #if we have limma results
        if (!(is.null(vals$gsvaLimma))){
          tempdf <- DataFrame(t(vals$gsvaRes[vals$gsvaLimma$Pathway[1:input$pickNtopPaths], , drop = FALSE]))
        } else {
          tempdf <- DataFrame(t(vals$gsvaRes[1:input$pickNtopPaths, , drop = FALSE]))
        }
        tempdf <- tempdf[, !(colnames(tempdf) %in% colnames(colData(vals$counts))), drop = FALSE]
        colData(vals$counts) <- cbind(colData(vals$counts), tempdf)
        updateColDataNames()
      }
    } else {
      shinyalert::shinyalert("Error!", "Run pathway first.", type = "error")
    }
  })

  #disable downloadPathway button if the pathway data doesn't exist
  isPathwayResult <- reactive(is.null(vals$gsvaRes))
  observe({
    if (isPathwayResult()) {
      shinyjs::disable("downloadPathway")
    } else {
      shinyjs::enable("downloadPathway")
    }
  })

  #download mast results
  output$downloadPathway <- downloadHandler(
    filename = function() {
      paste("pathway_results-", Sys.Date(), ".csv", sep = "")
    },
    content = function(file) {
      utils::write.csv(vals$gsvaRes, file)
    }
  )

  #-----------------------------------------------------------------------------
  # Page 6.2 : Enrichment Analysis - EnrichR
  #-----------------------------------------------------------------------------

  enrichRfile <- reactive(read.csv(input$enrFile$datapath,
                                   header = input$header,
                                   sep = input$sep,
                                   quote = input$quote,
                                   row.names = 1))
  output$enrBioGenes <- renderUI({
    if (!is.null(vals$counts)) {
      selectInput("selEnrBioGenes", "Select Gene List(s):",
                  names(which(apply(rowData(vals$counts), 2, function(a) length(unique(a)) == 2) == TRUE)))
    }
  })
  dbs <- reactive({
    if (internetConnection){
      enrDatabases <- enrichR::listEnrichrDbs()$libraryName
    } else {
      enrDatabases <- ""
    }
    if (is.null(input$enrichDb)){
      dbs <- enrDatabases
    } else {
      if (any(input$enrichDb %in% "ALL")){
        dbs <- enrDatabases
      } else {
        dbs <- input$enrichDb
      }
    }
  })
  #count_db <- reactive(length(dbs()))
  observeEvent (input$enrichRun, {
    if (is.null(vals$counts)){
      shinyalert::shinyalert("Error!", "Upload data first.", type = "error")
    } else {
      withBusyIndicatorServer ("enrichRun", {
        tryCatch ({
          if (input$geneListChoice == "selectGenes"){
            genes <- input$enrichGenes
          } else if (input$geneListChoice == "geneFile"){
            req(input$enrFile)
            genes <- rownames(enrichRfile())
          } else  {
            genes <- rownames(vals$counts)[SingleCellExperiment::rowData(vals$counts)[, input$selEnrBioGenes] == 1]
          }
          vals$enrichRes <- enrichRSCE(inSCE = vals$counts,
                                       glist = genes,
                                       db = dbs())
        }, error = function(e){
          shinyalert::shinyalert("Error!", e$message, type = "error")
        })
      })
    }
  })

  output$enrTabs <- renderUI({
    req(vals$enrichRes)
    nTabs <- length(dbs())
    #create tabPanel with datatable in it
    myTabs <- lapply(seq_len((nTabs)), function(i) {
      tabPanel(paste0(isolate(dbs()[i])),
               DT::dataTableOutput(paste0(isolate(dbs()[i])))
      )
    })
    do.call(tabsetPanel, myTabs)
  })

  #create datatables
  observe({
    req(vals$enrichRes)
    enrResults <- vals$enrichRes[, c(1:10)] %>%
      mutate(Database_selected =
               paste0("<a href='", vals$enrichRes[, 11],
                      "' target='_blank'>",
                      vals$enrichRes[, 1], "</a>"))
    lapply(seq_len(length(dbs())), function(i){
      output[[paste0(dbs()[i])]] <- DT::renderDataTable({
        DT::datatable({
          enr <- enrResults[which(vals$enrichRes[, 1] %in% dbs()[i]), ]
        }, escape = FALSE, options = list(scrollX = TRUE, pageLength = 30), rownames = FALSE)
      })
    })
  })

  #disable the downloadEnrichR button if the result doesn't exist
  isResult <- reactive(is.null(vals$enrichRes))
  observe({
    if (isResult()) {
      shinyjs::disable("downloadEnrichR")
    } else {
      shinyjs::enable("downloadEnrichR")
    }
  })

  output$downloadEnrichR <- downloadHandler(
    filename = function() {
      paste("enrichR-results-", Sys.Date(), ".csv", sep = "")
    },
    content = function(file) {
      utils::write.csv(vals$enrichRes, file)
    },
    contentType = "text/csv"
  )

  #-----------------------------------------------------------------------------
  # Page 7: Subsampling
  #-----------------------------------------------------------------------------

  #Run subsampling analysis
  observeEvent(input$runSubsampleDepth, {
    if (is.null(vals$counts)){
      shinyalert::shinyalert("Error!", "Upload data first.", type = "error")
    } else{
      withBusyIndicatorServer("runSubsampleDepth", {
        vals$subDepth <- DownsampleDepth(originalData = vals$counts,
                                         useAssay = input$depthAssay,
                                         minCount = input$minCount,
                                         minCells = input$minCells,
                                         maxDepth = 10 ^ input$maxDepth,
                                         realLabels = input$selectReadDepthCondition,
                                         depthResolution = input$depthResolution,
                                         iterations = input$iterations)

        output$depthDone <- renderPlot({
          plot(apply(vals$subDepth[, , 1], 2, median)~
                 seq(from = 0, to = input$maxDepth, length.out = input$depthResolution),
               lwd = 4, xlab = "log10(Total read counts)", ylab = "Number of detected genes",
               main = "Number of dected genes by sequencing depth")
          lines(apply(vals$subDepth[, , 1], 2, function(x){quantile(x, 0.25)})~
                  seq(from = 0, to = input$maxDepth, length.out = input$depthResolution), lty = 2, lwd = 3)
          lines(apply(vals$subDepth[, , 1], 2, function(x){quantile(x, 0.75)})~
                  seq(from = 0, to = input$maxDepth, length.out = input$depthResolution), lty = 2, lwd = 3)
        })
        output$minEffectDone <- renderPlot({
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
      shinyalert::shinyalert("Error!", "Upload data first.", type = "error")
    } else{
      withBusyIndicatorServer("runSubsampleCells", {
        if (input$useReadCount){
          vals$subCells <- DownsampleCells(originalData = vals$counts,
                                           useAssay = input$cellsAssay,
                                           realLabels = input$selectCellNumCondition,
                                           totalReads = sum(SummarizedExperiment::assay(vals$counts, input$cellsAssay)),
                                           minCellnum = input$minCellNum,
                                           maxCellnum = input$maxCellNum,
                                           minCountDetec = input$minCount,
                                           minCellsDetec = input$minCells,
                                           depthResolution = input$depthResolution,
                                           iterations = input$iterations)
        }
        else{
          vals$subCells <- DownsampleCells(originalData = vals$counts,
                                           useAssay = input$cellsAssay,
                                           realLabels = input$selectCellNumCondition,
                                           totalReads = input$totalReads,
                                           minCellnum = input$minCellNum,
                                           maxCellnum = input$maxCellNum,
                                           minCountDetec = input$minCount,
                                           minCellsDetec = input$minCells,
                                           depthResolution = input$depthResolution,
                                           iterations = input$iterations)
        }
        output$cellsDone <- renderPlot({
          plot(apply(vals$subCells[, , 1], 2, median)~
                 seq(from = input$minCellNum, to = input$maxCellNum, length.out = input$depthResolution),
               lwd = 4, xlab = "Number of virtual cells", ylab = "Number of detected genes",
               main = "Number of dected genes by cell number")
          lines(apply(vals$subCells[, , 1], 2, function(x){quantile(x, 0.25)})~
                  seq(from = input$minCellNum, to = input$maxCellNum, length.out = input$depthResolution), lty = 2, lwd = 3)
          lines(apply(vals$subCells[, , 1], 2, function(x){quantile(x, 0.75)})~
                  seq(from = input$minCellNum, to = input$maxCellNum, length.out = input$depthResolution), lty = 2, lwd = 3)
        })
        output$minEffectCells <- renderPlot({
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
      shinyalert::shinyalert("Error!", "Upload data first.", type = "error")
    } else{
      withBusyIndicatorServer("runSnapshot", {
        vals$snapshot <- iterateSimulations(originalData = vals$counts,
                                            useAssay = input$snapshotAssay,
                                            realLabels = input$selectSnapshotCondition,
                                            totalReads = input$numReadsSnap,
                                            cells = input$numCellsSnap,
                                            iterations = input$iterationsSnap)
        vals$effectSizes <- calcEffectSizes(countMatrix = SummarizedExperiment::assay(vals$counts, input$snapshotAssay), condition = colData(vals$counts)[, input$selectSnapshotCondition])
        output$Snaplot <- renderPlot({
          plot(apply(vals$snapshot, 1, function(x){sum(x <= 0.05) / length(x)}) ~ vals$effectSizes,
               xlab = "Cohen's d effect size", ylab = "Detection power", lwd = 4, main = "Power to detect diffex by effect size")
        })
      })
    }
  })
})
