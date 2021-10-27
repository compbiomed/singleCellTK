#1GB max upload size
options(shiny.maxRequestSize = 1000 * 1024 ^ 2)
options(useFancyQuotes = FALSE)
options(shiny.autoreload = TRUE)

internetConnection <- suppressWarnings(Biobase::testBioCConnection())
source("partials.R", local = TRUE) # creates several smaller UI components
# R.utils::sourceDirectory("qc_help_pages")
source("qc_help_pages/ui_decontX_help.R", local = TRUE) # creates several smaller UI components
source("qc_help_pages/ui_cxds_help.R", local = TRUE) # creates several smaller UI components
source("qc_help_pages/ui_bcds_help.R", local = TRUE) # creates several smaller UI components
source("qc_help_pages/ui_cxds_bcds_hybrid_help.R", local = TRUE) # creates several smaller UI components
source("qc_help_pages/ui_doubletFinder_help.R", local = TRUE) # creates several smaller UI components
source("qc_help_pages/ui_scrublet_help.R", local = TRUE) # creates several smaller UI components
source("qc_help_pages/ui_dc_and_qcm_help.R", local = TRUE) # creates several smaller UI components
source("qc_help_pages/ui_scDblFinder_help.R", local = TRUE) # creates several smaller UI components
# source("server_partials/server_01_data.R", local = TRUE) # functions for Data section

# Define server logic required to draw a histogram
shinyServer(function(input, output, session) {

  # PushBar setup
  # setup_pushbar(blur = FALSE, overlay = FALSE)

  # library(fs)
  # library(shinyFiles)

  #-----------------------------------------------------------------------------
  # MISC - Used throughout app
  #-----------------------------------------------------------------------------

  #reactive values object
  vals <- reactiveValues(
    counts = getShinyOption("inputSCEset"),
    original = getShinyOption("inputSCEset"),
    batchRes = NULL,
    gsvaRes = NULL,
    gsvaLimma = NULL,
    vamRes = NULL,
    vamCdf = NULL,
    vamResults = NULL,
    vamScore = NULL,
    gsvaScore = NULL,
    gsvaResults = NULL,
    visplotobject = NULL,
    enrichRes = NULL,
    celdaMod = NULL,
    celdaList = NULL,
    celdaListAll = NULL,
    celdaListAllNames = NULL,
    celdatSNE = NULL,
    celdaModuleFeature = NULL,
    dimRedPlot = NULL,
    dimRedPlot_geneExp = NULL,
    dendrogram = NULL,
    pcX = NULL,
    pcY = NULL,
    showAssayDetails = FALSE,
    hmCSPresets = list("RWB" = c("blue", "white", "red"),
                       "RdBu_r" = c("#2971b1", "#f7f6f6", "#b92732"),
                       "BrBG" = c("#995d12", "#f4f4f4", "#0c7068"),
                       "Blues" = c("#dae8f5", "#6daed4", "#0b559f"),
                       "Greens" = c("#E1F3F6", "#69C2A1", "#04702F")),
    hmCSURL = NULL,
    hmTmpColData = NULL,
    hmTmpRowData = NULL,
    hvgCalculated = list(status = FALSE, method = NULL),
    fmHMshowHide = FALSE
  )

  #Update all of the columns that depend on pvals columns
  updateColDataNames <- function(){
    pdataOptions <- colnames(colData(vals$counts))

    updateSelectInput(session, "qcSampleSelect", choices = pdataOptions)
    updateSelectInput(session, "filteredSample",
                      choices = c("none", pdataOptions))
    updateSelectInput(session, "deleterowdatacolumn",
                      choices = pdataOptions)
    updateSelectInput(session, "colorBy",
                      choices = c("No Color", "Gene Expression", pdataOptions))
    updateSelectInput(session, "shapeBy",
                      choices = c("No Shape", pdataOptions))
    updateSelectInput(session, "scMergeCT",
                      choices = c(pdataOptions))
    updateSelectInput(session, "combatCond",
                      choices = pdataOptions)
    updateSelectInput(session, "combatBioCond",
                      choices = c("None", pdataOptions))
    updateSelectInput(session, "batchCorrVar",
                      choices = pdataOptions)
    updateSelectInput(session, "batchCheckVar",
                      choices = pdataOptions)
    updateSelectInput(session, "batchCheckCond",
                      choices = c("None", pdataOptions))
    updateSelectInput(session, "clustVisCol", choices = pdataOptions)
    updateSelectInput(session, "deC1Class",
                      choices = pdataOptions)
    updateSelectInput(session, "deC2G1Col",
                      choices = pdataOptions)
    updateSelectInput(session, "deC2G2Col",
                      choices = pdataOptions)
    updateSelectInput(session, 'deCovar', choices = pdataOptions)
    updateSelectInput(session, "deHMcolData",
                      choices = pdataOptions)
    updateSelectInput(session, "deHMSplitCol",
                      choices = c('condition', pdataOptions),
                      selected = 'condition')
    updateSelectInput(session, "fmCluster", choices = pdataOptions)
    updateSelectInput(session, "fmHMcolData",
                      choices = pdataOptions)
    updateSelectInput(session, "hmCellAnn", choices = pdataOptions)
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
    updateSelectInput(session, "hmCellCol",
                      choices = pdataOptions)
    updateSelectInput(session, "hmCellTextBy",
                      choices = c("Row Names", pdataOptions))
    updateSelectInput(session, 'hmAddCellLabel',
                      choices = c("Default cell IDs", pdataOptions))
    updateSelectInput(session, "ctLabelByCluster",
                      choices = pdataOptions)
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
    selectRowData <- colnames(rowData(vals$counts))
    updateSelectInput(session, "gsByParam",
                      choices = c("rownames", selectRowData))
    updateSelectInput(session, "importFeatureDispOpt",
                      choices = c("Rownames (Default)", selectRowData))
    updateSelectInput(session, "filteredFeature",
                      choices = c("none", selectRowData))
    updateSelectInput(session, "deHMrowData",
                      choices = selectRowData)
    updateSelectInput(session, "deHMSplitRow",
                      choices = c('regulation', selectRowData),
                      selected = 'regulation')
    updateSelectInput(session, 'deVioLabel',
                      choices = c('Default ID', selectRowData))
    updateSelectInput(session, 'deRegLabel',
                      choices = c('Default ID', selectRowData))
    updateSelectInput(session, "fmHMrowData",
                      choices = selectRowData)
    updateSelectInput(session, "hmGeneCol",
                      choices = selectRowData)
    updateSelectInput(session, "hmGeneTextBy",
                      choices = c("Row Names", selectRowData))
    updateSelectInput(session, 'hmGeneAnn', choices = selectRowData)
    updateSelectInput(session, 'hmAddGeneLabel',
                      choices = c("Default feature IDs", selectRowData))
  }

  updateNumSamples <- function(){
    numsamples <- ncol(vals$counts)
    updateNumericInput(session, "downsampleNum", value = numsamples,
                       max = numsamples)
  }

  
  updateSelectInputTag <- function(session, inputId, choices = NULL, selected = NULL,
                                   label = "Select assay:", tags = NULL, recommended = NULL, showTags = TRUE,
                                   redDims = FALSE){
    choices <- expTaggedData(vals$counts, tags, redDims = redDims, showTags = showTags, recommended = recommended)
    updateSelectizeInput(session = session, inputId = inputId, label = label, choices = choices, selected = selected)
  }
  
  

  observeEvent(input$hvgMethodFS,{
    req(vals$counts)
    updateAssayInputs()
  })

  updateAssayInputs <- function(){
    currassays <- names(assays(vals$counts))
    updateSelectInputTag(session, "dimRedAssaySelect",
                         label = "Select Input Matrix:",
                         recommended = c("hvg"),
                         choices = expDataNames(vals$counts))
    updateSelectInputTag(session, "dimRedAssaySelect_tsneUmap",
                         label = "Select Input Matrix:",
                         recommended = c("redDims"),
                         redDims = TRUE)
    updateSelectInputTag(session, "batchCheckAssay", choices = currassays)
    updateSelectInputTag(session, "batchCheckOrigAssay", choices = currassays)
    updateSelectInputTag(session, "clustScranSNNMat", label = "Select Input Matrix:",
                         choices = expDataNames(vals$counts),
                         recommended = "redDims", redDims = TRUE)
    if (is.null(input$deMethod)) {
      updateSelectInputTag(session, "deAssay", tags = c("raw", "transformed", "uncategorized", "normalized", "scaled"), recommended = c("transformed"))
    } else if (input$deMethod == "DESeq2") {
      updateSelectInputTag(session, "deAssay", tags = c("raw", "transformed", "uncategorized", "normalized", "scaled"), recommended = c("raw"))
    } else {
      updateSelectInputTag(session, "deAssay", tags = c("raw", "transformed", "uncategorized", "normalized", "scaled"), recommended = c("transformed"))
    }
    if (is.null(input$fmMethod)) {
      updateSelectInputTag(session, "fmAssay", recommended = c("transformed"))
    } else if (input$fmMethod == "DESeq2") {
      updateSelectInputTag(session, "fmAssay", recommended = c("raw"))
    } else {
      updateSelectInputTag(session, "fmAssay", recommended = c("transformed"))
    }
    updateSelectInputTag(session, "fmHMAssay", choices = currassays, selected = input$fmAssay)
    updateSelectInputTag(session, "pathwayAssay", recommended = c("transformed", "normalized", "scaled"))
    updateSelectInputTag(session, "vamAssay", recommended = c("transformed", "normalized", "scaled"))

    updateSelectInputTag(session, "modifyAssaySelect")
    updateSelectInputTag(session, "normalizeAssaySelect", label = "Select assay to normalize:", recommended = "raw")

    updateSelectInputTag(session, "seuratSelectNormalizationAssay", choices = currassays, showTags = FALSE)
    if(input$hvgMethodFS == "vst"){
      updateSelectInputTag(session, "assaySelectFS_Norm", recommended = c("raw"))
    }
    else{
      updateSelectInputTag(session, "assaySelectFS_Norm", recommended = c("transformed", "normalized"))
    }
    updateSelectInputTag(session, "filterAssaySelect", choices = currassays)
    updateSelectInputTag(session, "qcAssaySelect", recommended = "raw")
    updateSelectInputTag(session, "celdaAssay", choices = currassays)
    updateSelectInputTag(session, "celdaAssayGS", choices = currassays)
    updateSelectInputTag(session, "celdaAssaytSNE", choices = currassays)
    updateSelectInputTag(session, "celdaAssayProbabilityMap",
                         choices = currassays)
    updateSelectInputTag(session, "celdaAssayModuleHeatmap",
                         choices = currassays)
    updateSelectInputTag(session, "depthAssay", choices = currassays)
    updateSelectInputTag(session, "cellsAssay", choices = currassays)
    updateSelectInputTag(session, "snapshotAssay", choices = currassays)
    updateSelectInputTag(session, "exportAssay", choices = currassays)
    updateSelectInputTag(session, "hmAssay", recommended = "transformed")
    updateSelectInputTag(session, "ctLabelAssay", choices = currassays, recommended = c("transformed"))
    # batch correction assay conditions
    bc.recommended <- NULL
    method.log <- c("FastMNN", "Limma", "MNN")
    method.scale <- c("BBKNN")
    method.raw <- c("ZINBWaVE", "ComBatSeq")
    if (is.null(input$batchCorrMethods)) {
      bc.recommended <- "raw"
    } else if (input$batchCorrMethods %in% method.log) {
      bc.recommended <- c("transformed")
    } else if (input$batchCorrMethods %in% method.raw) {
      bc.recommended <- "raw"
    } else if (input$batchCorrMethods %in% method.scale) {
      bc.recommended <- "scaled"
    }
    updateSelectInputTag(session, "batchCorrAssay",
                         label = "Select Assay to Correct:",
                         choices = currassays,
                         recommended = bc.recommended)
    updateSelectInputTag(session, "AdvancedMethodSelect_Colorby",
                         label = h5("Advanced Method"),
                         choices = currassays)
    updateSelectInputTag(session, "AdvancedMethodSelect_Xaxis",
                         label = h5("Advanced Method"),
                         choices = currassays)
    updateSelectInputTag(session, "AdvancedMethodSelect_Yaxis",
                         label = h5("Advanced Method"),
                         choices = currassays)
  }


  observeEvent(vals$counts, {
    # vals$counts
    if (!is.null(vals$counts)) {
      updateAssayInputs()
    }
  })

  observeEvent(vals$original, {
    if (!is.null(vals$original)) {
      #if (!is.null(metadata(vals$original)$sctk$genesets)) {
        #newGSchoices <- sctkListGeneSetCollections(vals$original)
        #updateSelectInput(session, "gsExisting", choices = c("None", newGSchoices))
        #updateSelectInput(session, "QCMgeneSets", choices =c("None", newGSchoices))
        #shinyjs::show(id = "gsAddToExisting", anim = FALSE)
      #} else {
        #shinyjs::hide(id = "gsAddToExisting", anim = FALSE)
        #updateSelectInput(session, "gsExisting", choices = c("None"), selected = "None")
        #updateSelectInput(session, "QCMgeneSets", choices =c("None"), selected = "None")
      #}
      shinyjs::show(id="combineOptions")
      #gsByChoices <- c("None", "rownames", names(rowData(vals$original)))
      #updateSelectInput(session, "gsByParam", choices = gsByChoices, selected = "rownames")
    } else {
      shinyjs::hide(id="combineOptions")
    }
  })

  updateReddimInputs <- function(){
    currreddim <- names(reducedDims(vals$counts))
    updateSelectInput(session, "FastMNNReddim", choices = currreddim)
    updateSelectInput(session, "HarmonyReddim", choices = currreddim)
    updateSelectInput(session, "clustVisReddim", choices = currreddim)
    updateSelectInput(session, "clustKMeansReddim", choices = currreddim)
    updateSelectInput(session, "clustSeuratReddim", choices = currreddim)
    updateSelectInput(session, "QuickAccess",
                      choices = c(currreddim, "Custom"))
    updateSelectInput(session, "ApproachSelect_Xaxis", choices = currreddim)
    updateSelectInput(session, "ApproachSelect_Yaxis", choices = currreddim)
    updateSelectInput(session, "ApproachSelect_Colorby", choices = currreddim)
  }

  updateEnrichDB <- function(){
    if (internetConnection){
      enrDB <- enrichR::listEnrichrDbs()$libraryName
    } else {
      enrDB <- ""
    }
    updateSelectInput(session, "enrichDb", choices = c("ALL", enrDB))
  }

  observeEvent(input$consoleToggle, {
    toggle(id = "console")
  })


  # js$disableTabs()
  # Close app on quit
  # session$onSessionEnded(stopApp)

  #-----------------------------------------------------------------------------
  # Page 1: Upload ####
  #-----------------------------------------------------------------------------
  sysname <- Sys.info()[['sysname']]
  if (sysname == "Windows") {
    roots <- getVolumes()()
  } else {
    roots <- c(home = "~/")
  }
  dirPaths <- reactiveValues(
    bDirectory = ".",
    sDirectory = ".",
    directory = ".",
    outputDirectory = "."
  )

  # Upload data through shiny app

  allImportEntries <- reactiveValues(samples=list(), id_count=0)

  shinyDirChoose(input, "bDirectory", roots = roots)
  shinyDirChoose(input, "sDirectory", roots = roots)
  shinyDirChoose(input, 'directory', roots = roots)

  output$bDirectoryPath <- renderText({
    dirPaths$bDirectory
  })
  output$sDirectoryPath <- renderText({
    dirPaths$sDirectory
  })
  output$directoryPath <- renderText({
    dirPaths$directory
  })

  # event listener for the base directory modal (need to populate table for sample names)
  # see https://github.com/wleepang/shiny-directory-input
  observeEvent(
    ignoreNULL = TRUE,
    eventExpr = {
      input$bDirectory
    },
    handlerExpr = {
      if ("path" %in% names(input$bDirectory)) {
        # condition prevents handler execution on initial app launch
        #path = choose.dir(default = readDirectoryInput(session, 'bDirectory'),
        #                  caption="Choose a directory")
        #updateDirectoryInput(session, 'bDirectory', value = path)

        vol <- roots[[input$bDirectory$root]]
        dirPaths$bDirectory <- paste0(vol, paste(unlist(input$bDirectory$path[-1]),
                                                 collapse = .Platform$file.sep))
        path <- dirPaths$bDirectory
        # clear the previous table of sample names
        prevPath <- path
        count <- 0
        for (prev in list.dirs(prevPath, recursive = FALSE)) {
          count <- count+1
          removeUI(
            selector = paste0("#sampleRow", count)
          )
        }
        # create a new table for the selected directory
        count <- 0
        if (!is.na(path)) {
          counts <- vector()
          for (sample in list.dirs(path, recursive = FALSE)) {
            count <- count+1
            counts <- c(counts, count)
            insertUI(
              selector = "#bDirTable",
              ui = fluidRow(
                id = paste0("sampleRow", count),
                column(6, basename(sample)),
                column(6, textAreaInput(paste0("sampleName", count), "Sample Name", resize = "none", value = basename(sample)))
              )
            )
          }
        }
      }
    }
  )

  # for sample directory modal
  observeEvent(
    ignoreNULL = TRUE,
    eventExpr = {
      input$sDirectory
    },
    handlerExpr = {
      #if (input$sDirectory > 0) {
      #  # condition prevents handler execution on initial app launch
      #  path = choose.dir(default = readDirectoryInput(session, 'sDirectory'),
      #                    caption="Choose a directory")
      #  updateDirectoryInput(session, 'sDirectory', value = path)
      #  if (!is.na(path)) {
      #    updateTextInput(session, "sSampleID", value = basename(path))
      #  }
      #}
      if ("path" %in% names(input$sDirectory)) {
        vol <- roots[[input$sDirectory$root]]
        dirPaths$sDirectory <- paste0(vol, paste(unlist(input$sDirectory$path[-1]),
                                                 collapse = .Platform$file.sep))
      }
    }
  )

  observeEvent(
    ignoreNULL = TRUE,
    eventExpr = {
      input$directory
    },
    handlerExpr = {
      #if (input$directory > 0) {
      #  # condition prevents handler execution on initial app launch
      #  path = choose.dir(default = readDirectoryInput(session, 'directory'),
      #                    caption="Choose a directory")
      #  updateDirectoryInput(session, 'directory', value = path)
      #}
      if ("path" %in% names(input$directory)) {
        vol <- roots[[input$directory$root]]
        dirPaths$directory <- paste0(vol, paste(unlist(input$directory$path[-1]),
                                                collapse = .Platform$file.sep))
      }
    }
  )

  # event listeners for "Add Sample" buttons
  observeEvent(input$addCR2Sample, {
    showModal(importCRModal())
  })
  observeEvent(input$crOpt1, {
    removeModal()
    showModal(importCRBDir())
  })
  observeEvent(input$crOpt2, {
    removeModal()
    showModal(importCRSDir())
  })
  observeEvent(input$crOpt3, {
    removeModal()
    showModal(importCRDDir())
  })
  observeEvent(input$addCR3Sample, {
    showModal(importCRModal())
  })
  observeEvent(input$addSSSample, {
    showModal(importModal())
  })
  observeEvent(input$addBUSSample, {
    showModal(importModal())
  })
  observeEvent(input$addSEQSample, {
    showModal(importModal(needsDir = TRUE))
  })
  observeEvent(input$addOptSample, {
    showModal(importModal())
  })

  # event listener for "Remove Sample" buttons
  observeEvent(input$clearAllImport, {
    for (entry in allImportEntries$samples) {
      removeUI(selector = paste0("#", entry$id))
    }
    allImportEntries$samples <- list()
  })

  # base directory
  observeEvent(input$BDirOK, {
    basePath <- dirPaths$bDirectory
    # if the user doesn't specify a base directory, show the modal again with the warning message
    if (identical(basePath, character(0))) {
      showModal(importCRBDir(failed = TRUE))
    } else {
      allDirs <- list.dirs(basePath, recursive = FALSE)
      # if we are adding a new CellRangerV2 sample
      if (input$algoChoice == "cellRanger2") {
        allUI <- vector()
        allIDs <- vector()
        count <- 0
        for (sample in allDirs) {
          count <- count + 1
          name <- input[[paste0("sampleName", count)]]
          if (!nzchar(name)) {
            name <- basename(sample)
          }
          id <- paste0("bnewSampleCR2", allImportEntries$id_count)
          entry <- list(type="cellRanger2", id=id, params=list(cellRangerDirs = basePath, sampleDirs = basename(sample), sampleNames = name))
          allImportEntries$samples <- c(allImportEntries$samples, list(entry))
          fluidRowStyle <- paste0(paste0("#", id), "{border-bottom: 1px solid #bababa; padding-top: .9%; padding-bottom: .5%}")
          removeBtnStyle <- paste0(paste0("#remove", id), "{padding-top: 0; padding-bottom: 0;}")
          ui_i <- fluidRow(
            id = id,
            tags$style(HTML(paste0(fluidRowStyle, removeBtnStyle))),
            column(3, basePath),
            column(3, basename(sample)),
            column(3, name),
            column(3, actionButton(paste0("remove", id), "X"))
          )
          allImportEntries$id_count <- allImportEntries$id_count + 1
          allUI <- c(allUI, list(ui_i))
          allIDs <- c(allIDs, id)
        }
      } else { # if we are adding a new CellRangerV3 sample
        allUI <- vector()
        allIDs <- vector()
        count <- 0
        for (sample in allDirs) {
          count <- count + 1
          name <- input[[paste0("sampleName", count)]]
          if (!nzchar(name)) {
            name <- basename(sample)
          }
          id <- paste0("bnewSampleCR3", allImportEntries$id_count)
          entry <- list(type="cellRanger3", id=id, params=list(cellRangerDirs = basePath, sampleDirs = basename(sample), sampleNames = name))
          allImportEntries$samples <- c(allImportEntries$samples, list(entry))
          fluidRowStyle <- paste0(paste0("#", id), "{border-bottom: 1px solid #bababa; padding-top: .9%; padding-bottom: .5%}")
          removeBtnStyle <- paste0(paste0("#remove", id), "{padding-top: 0; padding-bottom: 0;}")
          ui_i <- fluidRow(
            id = id,
            tags$style(HTML(paste0(fluidRowStyle, removeBtnStyle))),
            column(3, basePath),
            column(3, basename(sample)),
            column(3, name),
            column(3, actionButton(paste0("remove", id), "X"))
          )
          allImportEntries$id_count <- allImportEntries$id_count + 1
          allUI <- c(allUI, list(ui_i))
          allIDs <- c(allIDs, id)
        }
      }
      # insert all the new sample rows
      for (i in seq_along(allUI)) {
        insertUI(
          selector = "#newSampleImport",
          ui = allUI[i]
        )
      }
      # create event handlers for all the remove buttons
      # from: https://stackoverflow.com/questions/40038749/r-shiny-how-to-write-loop-for-observeevent
      lapply(
        X = allIDs,
        FUN = function(id_i){
          observeEvent(input[[paste0("remove", id_i)]], {
            removeUI(
              selector = paste0("#", id_i)
            )
            toRemove <- vector()
            for (entry in allImportEntries$samples) {
              if (entry$id == id_i) {
                toRemove <- c(toRemove, FALSE)
              } else {
                toRemove <- c(toRemove, TRUE)
              }
            }
            allImportEntries$samples <- allImportEntries$samples[toRemove]
          })
        }
      )
      removeModal()
    }
  })

  # event listeners for Cell Ranger import modals' OK buttons
  # sample directory
  observeEvent(input$SDirOK, {
    samplePath <- dirPaths$sDirectory
    # make sure a directory is selected
    if (identical(samplePath, character(0))) {
      showModal(importCRSDir(failed = TRUE))
    } else {
      # add the files to the appropriate reactiveValues
      if (input$algoChoice == "cellRanger2") {
        id <- paste0("snewSampleCR2", allImportEntries$id_count)
        entry <- list(type="cellRanger2", id=id, params=list(cellRangerDirs = dirname(samplePath), sampleDirs = basename(samplePath), sampleNames = input$sSampleID))
        allImportEntries$samples <- c(allImportEntries$samples, list(entry))
        allImportEntries$id_count <- allImportEntries$id_count + 1
      } else {
        id <- paste0("snewSampleCR3", allImportEntries$id_count)
        entry <- list(type="cellRanger3", id=id, params=list(cellRangerDirs = paste0(dirname(samplePath), "/"), sampleDirs = basename(samplePath), sampleNames = input$sSampleID))
        allImportEntries$samples <- c(allImportEntries$samples, list(entry))
        allImportEntries$id_count <- allImportEntries$id_count + 1
      }
      # add new row to table
      addToGeneralSampleTable(input$algoChoice, id, samplePath, input$sSampleID)
      # handler to remove the sample that was just added
      observeEvent(input[[paste0("remove", id)]],{
        removeUI(
          selector = paste0("#", id)
        )
        toRemove <- vector()
        for (entry in allImportEntries$samples) {
          if (entry$id == id) {
            toRemove <- c(toRemove, FALSE)
          } else {
            toRemove <- c(toRemove, TRUE)
          }
        }
        allImportEntries$samples <- allImportEntries$samples[toRemove]
      })
      removeModal()
    }
  })

  # data directory
  observeEvent(input$DDirOK, {
    dataPath <- dirPaths$directory
    if ((!nzchar(input$dSampleID)) || (identical(dataPath, character(0)))) {
      showModal(importCRDDir(failed = TRUE))
    } else {
      if (input$algoChoice == "cellRanger2") {
        id <- paste0("dnewSampleCR2", allImportEntries$id_count)
        entry <- list(type="cellRanger2", id=id, params=list(dataDir = dataPath, sampleName = input$dSampleID))
        allImportEntries$samples <- c(allImportEntries$samples, list(entry))
        allImportEntries$id_count <- allImportEntries$id_count + 1
      } else {
        id <- paste0("dnewSampleCR3", allImportEntries$id_count)
        entry <- list(type="cellRanger3", id=id, params=list(dataDir = dataPath, sampleName = input$dSampleID))
        allImportEntries$samples <- c(allImportEntries$samples, list(entry))
        allImportEntries$id_count <- allImportEntries$id_count + 1
      }
      # add new row to table
      addToGeneralSampleTable(input$algoChoice, id, dataPath, input$dSampleID)
      observeEvent(input[[paste0("remove", id)]],{
        removeUI(
          selector = paste0("#", id)
        )
        toRemove <- vector()
        for (entry in allImportEntries$samples) {
          if (entry$id == id) {
            toRemove <- c(toRemove, FALSE)
          } else {
            toRemove <- c(toRemove, TRUE)
          }
        }
        allImportEntries$samples <- allImportEntries$samples[toRemove]
      })
      removeModal()
    }
  })

  # event handler for pressing OK on the import modal
  observeEvent(input$modalOk, {
    basePath <- dirPaths$directory
    curFiles <- list()
    if ((!nzchar(input$sampleName)) || (identical(basePath, character(0)))) {
      showModal(importModal(failed = TRUE))
    } else {
      entry <- list()
      if (input$algoChoice == "starSolo") {
        id <- paste0("newSampleSS", allImportEntries$id_count)
        entry <- list(type="starSolo", id = id, params=list(STARsoloDirs = basePath, samples = input$sampleName))
        allImportEntries$samples <- c(allImportEntries$samples, list(entry))
        allImportEntries$id_count <- allImportEntries$id_count+1
      } else if (input$algoChoice == "busTools") {
        id <- paste0("newSampleBUS", allImportEntries$id_count)
        entry <- list(type="busTools", id = id, params=list(BUStoolsDirs = basePath, samples = input$sampleName))
        allImportEntries$samples <- c(allImportEntries$samples, list(entry))
        allImportEntries$id_count <- allImportEntries$id_count+1
      } else if (input$algoChoice == "seqc") {
        id <- paste0("newSampleSEQ", allImportEntries$id_count)
        entry <- list(type="seqc", id = id, params=list(seqcDirs = basePath, prefix = input$sampleID, samples = input$sampleName))
        updateTextInput(session, "sampleID", value = "")
        allImportEntries$samples <- c(allImportEntries$samples, list(entry))
        allImportEntries$id_count <- allImportEntries$id_count+1
      } else if (input$algoChoice == "optimus") {
        id <- paste0("newSampleOpt", allImportEntries$id_count)
        entry <- list(type="optimus", id = id, params=list(OptimusDirs = basePath, samples = input$sampleName))
        allImportEntries$samples <- c(allImportEntries$samples, list(entry))
        allImportEntries$id_count <- allImportEntries$id_count+1
      }
      addToGeneralSampleTable(input$algoChoice, id, basePath, input$sampleName)
      observeEvent(input[[paste0("remove", id)]],{
        removeUI(
          selector = paste0("#", id)
        )
        toRemove <- vector()
        for (entry in allImportEntries$samples) {
          if (entry$id == id) {
            toRemove <- c(toRemove, FALSE)
          } else {
            toRemove <- c(toRemove, TRUE)
          }
        }
        allImportEntries$samples <- allImportEntries$samples[toRemove]
      })
      removeModal()
    }
  })

  # Event handler to import a file input
  observeEvent(input$addFilesImport, {
    id <- paste0("newSampleFiles", allImportEntries$id_count)
    entry <- list(type="files", id = id, params=list(assayFile = input$countsfile$datapath, annotFile = input$annotFile$datapath,
                                                     featureFile = input$featureFile$datapath, assayName = input$inputAssayType))
    allImportEntries$samples <- c(allImportEntries$samples, list(entry))
    allImportEntries$id_count <- allImportEntries$id_count+1
    assayFileCol <- ""
    annotFileCol <- ""
    featureFileCol <- ""
    if (!is.null(input$countsfile$datapath)) {
      assayFileCol <- paste0("Assay: ", input$countsfile$datapath)
    }
    if (!is.null(input$annotFile$datapath)) {
      annotFileCol <- paste0("Annotation: ", input$annotFile$datapath)
    }
    if (!is.null(input$featureFile$datapath)) {
      featureFileCol <- paste0("Features: ", input$featureFile$datapath)
    }

    locCol <- paste(c(assayFileCol, annotFileCol, featureFileCol), collapse = "\n")

    addToGeneralSampleTable("files", id, locCol, input$inputAssayType)

    observeEvent(input[[paste0("remove", id)]],{
      removeUI(
        selector = paste0("#", id)
      )
      toRemove <- vector()
      for (entry in allImportEntries$samples) {
        if (entry$id == id) {
          toRemove <- c(toRemove, FALSE)
        } else {
          toRemove <- c(toRemove, TRUE)
        }
      }
      allImportEntries$samples <- allImportEntries$samples[toRemove]
    })
  })

  # Event handler to import an example input
  observeEvent(input$addExampleImport, {
    id <- paste0("newSampleExample", allImportEntries$id_count)
    entry <- list(type="example", id = id, params=list(dataset = input$selectExampleData))
    allImportEntries$samples <- c(allImportEntries$samples, list(entry))
    allImportEntries$id_count <- allImportEntries$id_count+1

    scRNAseqDatasets <- c("fluidigm_pollen", "allen_tasic")
    tenxPbmcDatasets <- c("pbmc3k", "pbmc4k", "pbmc6k", "pbmc8k", "pbmc33k", "pbmc68k")
    locCol <- ""
    if (input$selectExampleData %in% scRNAseqDatasets) {
      locCol <- "scRNA"
    } else {
      locCol <- "TENx"
    }

    addToGeneralSampleTable("example", id, locCol, input$selectExampleData)


    observeEvent(input[[paste0("remove", id)]],{
      removeUI(
        selector = paste0("#", id)
      )
      toRemove <- vector()
      for (entry in allImportEntries$samples) {
        if (entry$id == id) {
          toRemove <- c(toRemove, FALSE)
        } else {
          toRemove <- c(toRemove, TRUE)
        }
      }
      allImportEntries$samples <- allImportEntries$samples[toRemove]
    })
  })

  # Event handler to import an RDS input
  observeEvent(input$addRDSImport, {
    id <- paste0("newSampleRDS", allImportEntries$id_count)
    entry <- list(type="rds", id = id, params=list(rdsFile=input$rdsFile$datapath))
    allImportEntries$samples <- c(allImportEntries$samples, list(entry))
    allImportEntries$id_count <- allImportEntries$id_count+1

    addToGeneralSampleTable("rds", id, input$rdsFile$datapath, "")

    observeEvent(input[[paste0("remove", id)]],{
      removeUI(
        selector = paste0("#", id)
      )
      toRemove <- vector()
      for (entry in allImportEntries$samples) {
        if (entry$id == id) {
          toRemove <- c(toRemove, FALSE)
        } else {
          toRemove <- c(toRemove, TRUE)
        }
      }
      allImportEntries$samples <- allImportEntries$samples[toRemove]
    })
  })

  # Event handler for "Upload" button on import page
  observeEvent(input$uploadData, {
    withBusyIndicatorServer("uploadData", {
      if (length(allImportEntries$samples) == 0) {
        stop("You have not selected any samples to import.")
      }
      sceObj <- importMultipleSources(allImportEntries)
      if (input$combineSCEChoice == "addToExistingSCE") {
        if(!is.null(vals$original)) {
          sceList <- list(vals$original, sceObj)
          vals$original <- combineSCE(sceList = sceList,
                                      by.r = NULL,
                                      by.c = Reduce(intersect, lapply(sceList, function(x) { colnames(colData(x))})),
                                      combined = TRUE)
        } else {
          vals$original <- sceObj
        }
      } else {
        vals$original <- sceObj
      }
      # clear table and empty reactive
      for (entry in allImportEntries$samples) {
        removeUI(selector = paste0("#", entry$id))
      }
      allImportEntries$samples <- list()

      # Add sample variable if it was not included
      if (!"sample" %in% names(colData(vals$original)) &&
          !"Sample" %in% names(colData(vals$original))) {
        sampleVar <- "sample"
        colData(vals$original)$sample = sampleVar
      } else if ("sample" %in% names(colData(vals$original))) {
        sampleVar <- "sample"
      } else {
        sampleVar <- "Sample"
      }

      if (!is.null(vals$original)) {
        vals$counts <- vals$original
        #store assayType information in the metadata
        # if (!"assayType" %in% names(metadata(vals$counts))) {
        #   vals$counts <- expSetDataTag(
        #     inSCE = vals$counts,
        #     assayType = "raw",
        #     assays = assayNames(vals$counts))
        # }
        if (any(duplicated(rownames(vals$counts)))) {
          warning("Duplicated rownames detected, making them unique...")
          vals$counts <- dedupRowNames(vals$counts)
        }
        # ToDo: Remove these automatic updates and replace with
        # observeEvents functions that activate upon the tab selection
        updateColDataNames()
        updateFeatureAnnots()
        updateNumSamples()
        # updateAssayInputs()
        updateGeneNames()
        updateReddimInputs()
        shinyjs::show(id="annotationData")
        js$enableTabs();
      } else {
        shinyalert::shinyalert("Error!", "The data upload failed!",
                               type = "error")
      }
      vals$gsvaRes <- NULL
      vals$vamRes <- NULL
      vals$vamResults <- NULL
      vals$gsvaResults <- NULL
      vals$gsvaLimma <- NULL
      vals$vamScore <- NULL
      vals$gsvaScore <- NULL
      vals$visplotobject <- NULL
      vals$enrichRes <- NULL
      vals$dimRedPlot <- NULL
      vals$dimRedPlot_geneExp <- NULL
      vals$dendrogram <- NULL
      vals$pcX <- NULL
      vals$pcY <- NULL
      vals$batchRes <- NULL
      vals$hvgCalculated <- list(status = FALSE, method = NULL)
      dbList <- getMSigDBTable()
      geneSetDBChoices <- formatGeneSetDBChoices(dbIDs = dbList$ID, dbCats = dbList$Category_Description)
      updateCheckboxGroupInput(session, 'geneSetDB', choices = geneSetDBChoices)

      updateSeuratUIFromRDS(vals$counts)
      cleanGSTable()
      # TODO: There are more things that need to be cleaned when uploading new
      # dataset, including any plots, tables that are origined from the old
      # datasets. Otherwise, errors may pop out when Shiny listens to the new
      # object but cannot find the old result.
    })
    callModule(module = nonLinearWorkflow, id = "nlw-import", parent = session, qcf = TRUE)
  })

  updateSeuratUIFromRDS <- function(inSCE){
    if(!is.null(metadata(inSCE)$seurat$plots)){
      showNotification(HTML("Computation from Seurat Report detected in the input object, therefore the toolkit will now populate the Seurat tab with computated data & plots for further inspection. Click on the button below to directly go the the Seurat tab of the toolkit now! <br><br>"),
                       type = "message", duration = 0, action = actionBttn(
        inputId = "goToSeurat",
        label = "Go to Seurat Curated Workflow",
        style = "bordered",
        color = "royal",
        size = "s",
        icon = icon("arrow-right")
      ), id = "goSeuratNotification")

      #Normalize Data
      shinyjs::enable(selector = "div[value='Normalize Data']")
      updateCollapse(session = session, "SeuratUI", style = list("Normalize Data" = "success"))
      normalizeParams <- metadata(vals$counts)$seurat$sctk$report$normalizeParams
      updateSelectInput(session, "normalization_method", selected = normalizeParams$normalizationMethod)
      updateTextInput(session, "scale_factor", value = normalizeParams$scaleFactor)

      #Scale Data
      shinyjs::enable(selector = "div[value='Scale Data']")
      updateCollapse(session = session, "SeuratUI", style = list("Scale Data" = "success"))
      scaleParams <- metadata(vals$counts)$seurat$sctk$report$scaleParams
      updateSelectInput(session, "model.use", selected = scaleParams$model)

      #HVG
      hvgParams <- metadata(vals$counts)$seurat$sctk$report$hvgParams
      output$plot_hvg <- renderPlotly({
        isolate({
          plotly::ggplotly(seuratPlotHVG(vals$counts, labelPoints = hvgParams$labelPoints))
        })
      })
      shinyjs::enable(selector = "div[value='Highly Variable Genes']")
      updateCollapse(session = session, "SeuratUI", style = list("Highly Variable Genes" = "success"))
      updateSelectInput(session, "hvg_method", selected = hvgParams$hvgMethod)
      updateTextInput(session, "hvg_no_features", value = hvgParams$hvgNumber)
      updateTextInput(session, "hvg_no_features_view", value = hvgParams$labelPoints)

      #DR
      pcaParams <- metadata(vals$counts)$seurat$sctk$report$pcaParams
      shinyjs::enable(selector = "div[value='Dimensionality Reduction']")
      updateCollapse(session = session, "SeuratUI", style = list("Dimensionality Reduction" = "success"))

      removeTab(inputId = "seuratPCAPlotTabset", target = "PCA Plot")
      removeTab(inputId = "seuratPCAPlotTabset", target = "Elbow Plot")
      removeTab(inputId = "seuratPCAPlotTabset", target = "JackStraw Plot")
      removeTab(inputId = "seuratPCAPlotTabset", target = "Heatmap Plot")

      shinyjs::show(selector = ".seurat_pca_plots")

      appendTab(inputId = "seuratPCAPlotTabset", tabPanel(title = "PCA Plot",
                                                          panel(heading = "PCA Plot",
                                                                plotlyOutput(outputId = "plot_pca")
                                                          )
      ), select = TRUE)
      appendTab(inputId = "seuratPCAPlotTabset", tabPanel(title = "Elbow Plot",
                                                          panel(heading = "Elbow Plot",
                                                                plotlyOutput(outputId = "plot_elbow_pca")
                                                          )
      ))
      appendTab(inputId = "seuratPCAPlotTabset", tabPanel(title = "JackStraw Plot",
                                                          panel(heading = "JackStraw Plot",
                                                                plotlyOutput(outputId = "plot_jackstraw_pca")
                                                          )
      ))
      appendTab(inputId = "seuratPCAPlotTabset", tabPanel(title = "Heatmap Plot",
                                                          panel(heading = "Heatmap Plot",
                                                                panel(heading = "Plot Options",
                                                                      fluidRow(
                                                                        column(6,
                                                                               pickerInput(inputId = "picker_dimheatmap_components_pca", label = "Select principal components to plot:", choices = c(), options = list(`actions-box` = TRUE, size = 10, `selected-text-format` = "count > 3"), multiple = TRUE)
                                                                        ),
                                                                        column(6,
                                                                               sliderInput(inputId = "slider_dimheatmap_pca", label = "Number of columns for the plot: ", min = 1, max = 4, value = 2)
                                                                        )
                                                                      ),
                                                                      actionButton(inputId = "plot_heatmap_pca_button", "Plot")
                                                                ),
                                                                panel(heading = "Plot",
                                                                      shinyjqui::jqui_resizable(plotOutput(outputId = "plot_heatmap_pca"), options = list(maxWidth = 700))
                                                                )
                                                          )
      ))


        # output$plot_pca <- renderPlotly({
        #   plotly::ggplotly(metadata(inSCE)$seurat$plots$pca)
        # })

      output$plot_pca <- renderPlotly({
        isolate({
          plotly::ggplotly(seuratReductionPlot(inSCE = vals$counts, useReduction = "pca"))
        })
      })


          #updateNumericInput(session = session, inputId = "pca_significant_pc_counter", value = singleCellTK:::.computeSignificantPC(vals$counts))
          # output$plot_elbow_pca <- renderPlotly({
          #   metadata(inSCE)$seurat$plots$elbow
          # })

          #update parameters from seurat report
          output$plot_elbow_pca <- renderPlotly({
            isolate({
              plotly::ggplotly(seuratElbowPlot(inSCE = vals$counts))
            })
          })

          output$pca_significant_pc_output <- renderText({
            isolate({
              paste("<p>Number of significant components suggested by ElbowPlot: <span style='color:red'>", pcaParams$significant_PC," </span> </p> <hr>")
            })
          })

          # output$plot_jackstraw_pca <- renderPlotly({
          #   plotly::ggplotly(metadata(inSCE)$seurat$plots$jackstraw)
          # })

          output$plot_jackstraw_pca <- renderPlotly({
            isolate({
              plotly::ggplotly(seuratJackStrawPlot(vals$counts))
            })
          })


          # output$plot_heatmap_pca <- renderPlot({
          #   metadata(inSCE)$seurat$plots$heatmap
          # })

          updateTextInput(session, "pca_no_components", value = pcaParams$nPCs)
          updateMaterialSwitch(session, "pca_compute_jackstraw", value = TRUE)
          updateNumericInput(session, "pca_significant_pc_counter", value = pcaParams$significant_PC)

          pcHeatmapParams <- metadata(inSCE)$seurat$plots$heatmap
          pcHeatmapParams$inSCE <- vals$counts
          output$plot_heatmap_pca <- renderPlot({
            isolate({
              do.call("seuratComputeHeatmap", pcHeatmapParams)
            })
          })

          updatePickerInput(session = session, inputId = "picker_dimheatmap_components_pca", choices = singleCellTK:::.getComponentNames(vals$counts@metadata$seurat$count_pc, "PC"))



      #tSNE/UMAP
          shinyjs::enable(selector = "div[value='tSNE/UMAP']")
          updateCollapse(session = session, "SeuratUI", style = list("tSNE/UMAP" = "success"))

          # output$plot_tsne <- renderPlotly({
          #   metadata(inSCE)$seurat$plots$tsne
          # })
          #
          # output$plot_umap <- renderPlotly({
          #   metadata(inSCE)$seurat$plots$umap
          # })

          output$plot_tsne <- renderPlotly({
            isolate({
              plotly::ggplotly(seuratReductionPlot(vals$counts, useReduction = "tsne"))
            })
          })

          output$plot_umap <- renderPlotly({
            isolate({
              plotly::ggplotly(seuratReductionPlot(vals$counts, useReduction = "umap"))
            })
          })


      #Clustering
          clusterParams <- metadata(vals$counts)$seurat$sctk$report$clusterParams
          shinyjs::enable(selector = "div[value='Clustering']")
          updateCollapse(session = session, "SeuratUI", style = list("Clustering" = "success"))

          removeTab(inputId = "seuratClusteringPlotTabset", target = "PCA Plot")
          removeTab(inputId = "seuratClusteringPlotTabset", target = "ICA Plot")
          removeTab(inputId = "seuratClusteringPlotTabset", target = "tSNE Plot")
          removeTab(inputId = "seuratClusteringPlotTabset", target = "UMAP Plot")

          appendTab(inputId = "seuratClusteringPlotTabset", tabPanel(title = "PCA Plot",
                                                                     panel(heading = "PCA Plot",
                                                                           plotlyOutput(outputId = "plot_pca_clustering")
                                                                     )
          ), select = TRUE

          )

          output$plot_pca_clustering <- renderPlotly({
            plotly::ggplotly(seuratReductionPlot(vals$counts, useReduction = "pca", showLegend = TRUE))
          })

          appendTab(inputId = "seuratClusteringPlotTabset", tabPanel(title = "tSNE Plot",
                                                                     panel(heading = "tSNE Plot",
                                                                           plotlyOutput(outputId = "plot_tsne_clustering")
                                                                     )
          )
          )

          output$plot_tsne_clustering <- renderPlotly({
            plotly::ggplotly(seuratReductionPlot(vals$counts, useReduction = "tsne", showLegend = TRUE))
          })

          appendTab(inputId = "seuratClusteringPlotTabset", tabPanel(title = "UMAP Plot",
                                                                     panel(heading = "UMAP Plot",
                                                                           plotlyOutput(outputId = "plot_umap_clustering")
                                                                     )
          )
          )

          output$plot_umap_clustering <- renderPlotly({
            plotly::ggplotly(seuratReductionPlot(vals$counts, useReduction = "umap", showLegend = TRUE))
          })

          shinyjs::show(selector = ".seurat_clustering_plots")

          updateNumericInput(session, "resolution_clustering", value = clusterParams$resolution)


      #Find Markers
          shinyjs::enable(selector = "div[value='Find Markers']")
          updateCollapse(session = session, "SeuratUI", style = list("Find Markers" = "success"))

          shinyjs::show(selector = ".seurat_findmarker_table")
          shinyjs::show(selector = ".seurat_findmarker_jointHeatmap")
          shinyjs::show(selector = ".seurat_findmarker_plots")

          removeTab(inputId = "seuratFindMarkerPlotTabset", target = "Ridge Plot")
          removeTab(inputId = "seuratFindMarkerPlotTabset", target = "Violin Plot")
          removeTab(inputId = "seuratFindMarkerPlotTabset", target = "Feature Plot")
          removeTab(inputId = "seuratFindMarkerPlotTabset", target = "Dot Plot")
          removeTab(inputId = "seuratFindMarkerPlotTabset", target = "Heatmap Plot")

          appendTab(inputId = "seuratFindMarkerPlotTabset", tabPanel(title = "Ridge Plot",
                                                                     panel(heading = "Ridge Plot",
                                                                           shinyjqui::jqui_resizable(
                                                                             plotOutput(outputId = "findMarkerRidgePlot")
                                                                           )
                                                                     )
          )
          )
          appendTab(inputId = "seuratFindMarkerPlotTabset", tabPanel(title = "Violin Plot",
                                                                     panel(heading = "Violin Plot",
                                                                           shinyjqui::jqui_resizable(
                                                                             plotOutput(outputId = "findMarkerViolinPlot")
                                                                           )
                                                                     )
          )
          )
          appendTab(inputId = "seuratFindMarkerPlotTabset", tabPanel(title = "Feature Plot",
                                                                     panel(heading = "Feature Plot",
                                                                           shinyjqui::jqui_resizable(
                                                                             plotOutput(outputId = "findMarkerFeaturePlot")
                                                                           )
                                                                     )
          )
          )
          appendTab(inputId = "seuratFindMarkerPlotTabset", tabPanel(title = "Dot Plot",
                                                                     panel(heading = "Dot Plot",
                                                                           shinyjqui::jqui_resizable(
                                                                             plotOutput(outputId = "findMarkerDotPlot")
                                                                           )
                                                                     )
          )
          )
          appendTab(inputId = "seuratFindMarkerPlotTabset", tabPanel(title = "Heatmap Plot",
                                                                     panel(heading = "Heatmap Plot",
                                                                           fluidRow(
                                                                             column(12, align = "center",
                                                                                    panel(
                                                                                      plotOutput(outputId = "findMarkerHeatmapPlot")
                                                                                    )
                                                                             )
                                                                           )
                                                                     )
          )

          )

          showTab(inputId = "seuratFindMarkerPlotTabset", target = "Joint Heatmap Plot")
          updateTabsetPanel(session = session, inputId = "seuratFindMarkerPlotTabset", selected = "Ridge Plot")
          shinyjs::show(selector = ".seurat_findmarker_plots")

          groupHeatmapParams <- metadata(vals$counts)$seurat$plots$groupHeatmapParams
          groupHeatmapParams$inSCE <- vals$counts
          output$findMarkerHeatmapPlotFull <- renderPlot({
            isolate({
              do.call("seuratGenePlot", groupHeatmapParams)
            })
          })

          output$findMarkerHeatmapPlotFullTopText <- renderUI({
            h6(paste("Heatmap plotted across all groups against genes with adjusted p-values <", input$seuratFindMarkerPValAdjInput))
          })

          updateSelectInput(session, "seuratFindMarkerSelectPhenotype", choices = colnames(colData(vals$counts)), selected = metadata(vals$counts)$seurat$plots$group)


          vals$fts <- callModule(
            module = filterTableServer,
            id = "filterSeuratFindMarker",
            dataframe = metadata(vals$counts)$seurat$plots$top9
          )

      #Downstream Analysis
          shinyjs::show(selector = "div[value='Downstream Analysis']")
          updateCollapse(session = session, "SeuratUI", style = list("Downstream Analysis" = "info"))
    }
  }

  observeEvent(input$goToSeurat,{
    updateTabsetPanel(session, "navbar",
                      selected = "Seurat")
    removeNotification(id = "goSeuratNotification", session = session)
  })

  observeEvent(input$importFeatureDipSet, {
    if (!is.null(vals$counts)) {
      withBusyIndicatorServer("importFeatureDipSet", {
        if (!input$importFeatureDispOpt == "Rownames (Default)") {
          vals$counts <- setSCTKDisplayRow(vals$counts,
                                           input$importFeatureDispOpt)
        }
      })
    }
  })

  #-----------#
  # Gene Sets ####
  #-----------#
  numGS <- reactiveValues(id_count = 0)

  addToGSTable <- function(nameCol, locCol) {
    numGS$id_count <- numGS$id_count + 1
    id <- paste0("geneSet", numGS$id_count)
    fluidRowStyle <- paste0(paste0("#", id), "{border-bottom: 1px solid #bababa; padding-top: .9%; padding-bottom: .5%}")
    insertUI(
      selector = "#newGSImport",
      ui = fluidRow(
        id = id,
        tags$style(HTML(fluidRowStyle)),
        column(3, nameCol),
        column(9, locCol),
      )
    )
  }

  vals$defaultQCGS <- c("None" = "none",
                        "Human Mitochondrial Genes (Ensembl)" = "he",
                        "Human Mitochondrial Genes (Symbol)" = "hs",
                        "Mouse Mitochondrial Genes (Ensembl)" = "me",
                        "Mouse Mitochondrial Genes (Symbol)" = "ms")
  cleanGSTable <- function() {
    for (i in seq(numGS$id_count)) {
      removeUI(
        selector = paste0("#geneSet", i)
      )
    }
    numGS$id_count <- 0
    if (!is.null(vals$counts)) {
      existGS <- sctkListGeneSetCollections(vals$counts)
      if (length(existGS) > 0) {
        for (i in existGS) {
          addToGSTable(i, "SCE Object")
        }
        updateSelectInput(session, "gsExisting", choices = c("None", existGS))
        names(existGS) <- existGS
        updateSelectInput(session, "QCMgeneSets", choices =c(vals$defaultQCGS, existGS),
                          selected = "none")
        shinyjs::show(id = "gsAddToExisting", anim = FALSE)
      } else {
        shinyjs::hide(id = "gsAddToExisting", anim = FALSE)
      }
    } else {
      updateSelectInput(session, "gsExisting", choices = "None")
      updateSelectInput(session, "QCMgeneSets", choices = vals$defaultQCGS,
                        selected = "none")
      shinyjs::hide(id = "gsAddToExisting", anim = FALSE)
    }
  }

  handleGSPasteOption <- function(byParam) {
    if (!nzchar(input$geneSetText)) {
      shinyjs::show(id = "gsUploadError", anim = FALSE)
    } else if ((!nzchar(input$gsCollectionNameText)) && (input$gsExisting == "None")) {
      shinyjs::show(id = "gsUploadError", anim = FALSE)
    } else {
      shinyjs::hide(id = "gsUploadError", anim = FALSE)
      setList <- formatGeneSetList(input$geneSetText)
      if (nzchar(input$gsCollectionNameText)) {
        vals$counts <- importGeneSetsFromList(vals$counts,
                                                setList,
                                                by = byParam,
                                                collectionName = input$gsCollectionNameText)
        addToGSTable(input$gsCollectionNameText, "Paste-In")
      } else if (input$gsExisting != "None") {
        vals$counts <- importGeneSetsFromList(vals$counts,
                                                setList,
                                                by = byParam,
                                                collectionName = input$gsExisting)
        addToGSTable(input$gsExisting, "Paste-In")
      }
    }
  }

  observeEvent(input$uploadGS, {
    withBusyIndicatorServer("uploadGS", {
      byParam = NULL
      if (input$gsByParam != "None") {
        byParam <- input$gsByParam
      }
      if (input$geneSetSourceChoice == "gsGMTUpload") {
        if (is.null(input$geneSetGMT)) {
          shinyjs::show(id = "gsUploadError", anim = FALSE)
        } else if (!nzchar(input$gsCollectionNameGMT)){
          shinyjs::show(id = "gsUploadError", anim = FALSE)
        } else {
          shinyjs::hide(id = "gsUploadError", anim = FALSE)
          vals$counts <- importGeneSetsFromGMT(vals$counts,
                                                 input$geneSetGMT$datapath,
                                                 by = byParam,
                                                 collectionName = input$gsCollectionNameGMT)
          addToGSTable(input$gsCollectionNameGMT, input$geneSetGMT$datapath)
        }

      } else if (input$geneSetSourceChoice == "gsDBUpload") {
        if (is.null(input$geneSetDB)) {
          shinyjs::show(id = "gsUploadError", anim = FALSE)
        } else {
          shinyjs::hide(id = "gsUploadError", anim = FALSE)
          vals$counts <- importGeneSetsFromMSigDB(vals$counts,
                                                    input$geneSetDB,
                                                    by = byParam)
          for(i in input$geneSetDB){
            # Handling multiple selections from the checkboxInput
            addToGSTable(i, "Database")
          }
        }

      } else if (input$geneSetSourceChoice == "gsMito") {
        vals$counts <- importMitoGeneSet(vals$counts,
                                         reference = input$geneSetMitoSpecies,
                                         id = input$geneSetMitoID,
                                         by = byParam,
                                         collectionName = input$geneSetMitoName)
        addToGSTable(input$geneSetMitoName, "SCTK Curated Geneset")
      } else if (input$geneSetSourceChoice == "gsPasteUpload") {
        handleGSPasteOption(byParam)
      }

      allGS <- sctkListGeneSetCollections(vals$counts)
      updateSelectInput(session, "gsExisting", choices = c("None", allGS))
      names(allGS) <- allGS
      updateSelectInput(session, "QCMgeneSets", choices =c(vals$defaultQCGS, allGS),
                        selected = "none")
      shinyjs::show(id = "gsAddToExisting", anim = FALSE)
    })
  })

  #-----------------------------------------------------------------------------
  # Page 2: Data Summary and Filtering ####
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

  #----#
  # QC #####
  #----#
  # Hide and show parameters for QC functions
  shinyjs::onclick("QCMetrics", shinyjs::toggle(id = "QCMetricsParams",
                                                anim = FALSE), add = TRUE)
  shinyjs::onclick("decontX", shinyjs::toggle(id = "decontXParams",
                                              anim = FALSE), add = TRUE)
  shinyjs::onclick("scDblFinder", shinyjs::toggle(id = "scDblFinderParams",
                                                   anim = FALSE), add = TRUE)
  shinyjs::onclick("cxds", shinyjs::toggle(id = "cxdsParams",
                                           anim = FALSE), add = TRUE)
  shinyjs::onclick("bcds", shinyjs::toggle(id = "bcdsParams",
                                           anim = FALSE), add = TRUE)
  shinyjs::onclick("cxds_bcds_hybrid", shinyjs::toggle(id = "cxds_bcds_hybridParams",
                                                       anim = FALSE), add = TRUE)
  shinyjs::onclick("scrublet", shinyjs::toggle(id = "scrubletParams",
                                               anim = FALSE), add = TRUE)
  shinyjs::onclick("doubletFinder", shinyjs::toggle(id = "doubletFinderParams",
                                                    anim = FALSE), add = TRUE)

  qc_choice_list <- list("scDblFinder", "cxds", "bcds",
                         "cxds_bcds_hybrid", "decontX", "QCMetrics", "scrublet", "doubletFinder")
  # holds all the input ids for the QC algorithm parameters by algorithm name
  qc_input_ids <- list(scDblFinder = list(nNeighbors="DCnNeighbors", simDoublets="DCsimDoublets"),

                       cxds = list(ntop="CXntop", binThresh="CXbinThresh", verb="CXverb", retRes="CXretRes", estNdbl="CXestNdbl"),

                       bcds = list(ntop="BCntop", srat="BCsrat", verb="BCverb", retRes="BCretRes", nmax="BCnmax", varImp="BCvarImp", estNdbl="BCestNdbl"),

                       cxds_bcds_hybrid = list(cxdsArgs=list(ntop="CX2ntop", binThresh="CX2binThresh", retRes="CX2retRes"),
                                               bcdsArgs=list(ntop="BC2ntop", srat="BC2srat", retRes="BC2retRes", nmax="BC2nmax", varImp="BC2varImp"),
                                               verb="CXBCverb", estNdbl="CXBCestNdbl"),

                       decontX = list(maxIter="DXmaxIter", estimateDelta="DXestimateDelta", convergence="DXconvergence",
                                      iterLogLik="DXiterLogLik", varGenes="DXvarGenes", dbscanEps="DXdbscanEps", verbose="DXverbose"),

                       doubletFinder = list(seuratNfeatures="DFseuratNfeatures", seuratRes="DFseuratRes", formationRate="DFformationRate", verbose="DFverbose"),
                       scrublet = list(simDoubletRatio="SsimDoubletRatio", nNeighbors="SnNeighbors", minDist="SminDist", expectedDoubletRate="SexpectedDoubletRate",
                                       stdevDoubletRate='SstdevDoubletRate', syntheticDoubletUmiSubsampling="SsyntheticDoubletUmiSubsampling",
                                       useApproxNeighbors="SuseApproxNeighbors", distanceMetric="SdistanceMetric", getDoubletNeighborParents="SgetDoubletNeighborParents", minCounts="SminCounts",
                                       minCells="SminCells", minGeneVariabilityPctl="SminGeneVariabilityPctl", logTransform="SlogTransform", meanCenter="SmeanCenter",
                                       normalizeVariance="SnormalizeVariance", nPrinComps="SnPrinComps", tsneAngle="StsneAngle", tsnePerplexity="StsnePerplexity", verbose="Sverbose")
  )
  # to keep track of whether an algo has already been run
  qc_algo_status = reactiveValues(scDblFinder=NULL, cxds=NULL, bcds=NULL, cxds_bcds_hybrid=NULL, decontX=NULL,
                                  QCMetrics=NULL, scrublet=NULL, doubletFinder=NULL)

  qc_plot_ids = reactiveValues(scDblFinder="DCplots", cxds="CXplots", bcds="BCplots", cxds_bcds_hybrid="CXBCplots", decontX="DXplots",
                               QCMetrics="QCMplots", scrublet="Splots", doubletFinder="DFplots")


  # event handlers to open help pages for each qc algorithm
  observeEvent(input$DXhelp, {
    showModal(decontXHelpModal())
  })
  observeEvent(input$CXhelp, {
    showModal(cxdsHelpModal())
  })
  observeEvent(input$BChelp, {
    showModal(bcdsHelpModal())
  })
  observeEvent(input$CXBChelp, {
    showModal(cxdsBcdsHybridHelpModal())
  })
  observeEvent(input$DFhelp, {
    showModal(doubletFinderHelpModal())
  })
  observeEvent(input$Shelp, {
    showModal(scrubletHelpModal())
  })
  observeEvent(input$DChelp, {
    showModal(scDblFinderHelpModal())
  })
  observeEvent(input$QCMhelp, {
    showModal(QCMHelpModal())
  })
  observeEvent(input$QCImportGS, {
    showTab(inputId = "navbar",
            target = "Import Gene Sets",
            select = TRUE,
            session = session)
  })


  # format the parameters for decontX
  prepDecontXParams <- function(paramsList) {
    inputIds <- qc_input_ids[["decontX"]]
    dxParams <- list()
    # put in all the params from the list (the straightforward ones)
    for (key in names(inputIds)) {
      dxParams[[key]] = input[[inputIds[[key]]]]
    }

    # put in the delta params (c-bind the two priors)
    dxParams[["delta"]] <- c(input$DXnativePrior, input$DXcontPrior)

    # add to master params list
    paramsList[["decontX"]] = dxParams
    return(paramsList)
  }

  # format the parameters for doubletFinder
  prepDoubletFinderParams <- function(paramsList) {
    inputIds <- qc_input_ids[["doubletFinder"]]
    dfParams <- list()
    # put in all the params from the list (the straightforward ones)
    for (key in names(inputIds)) {
      dfParams[[key]] = input[[inputIds[[key]]]]
    }

    # put in the seuratPcs param (range from 1 to given value)
    dfParams[["seuratPcs"]] <- 1:input$DFseuratPcs

    # add to master params list
    paramsList[["doubletFinder"]] = dfParams
    return(paramsList)
  }

  qcInputExists <- function() {
    for (algo in qc_choice_list) {
      if (isTRUE(input[[algo]])) {
        return(TRUE)
      }
    }
    return(FALSE)
  }

  updateQCPlots <- function() {
    # get selected sample from run QC section
    if (!is.null(vals$counts)) {
      qcSample <- input$qcSampleSelect
      if (qcSample == "None") {
        qcSample <- NULL
      } else {
        qcSample <- colData(vals$counts)[,input$qcSampleSelect]
      }
      # build list of selected algos
      algoList = list()
      for (algo in qc_choice_list) {
        if (isTRUE(input[[algo]])) {
          algoList <- c(algoList, algo)
        }
      }
      # only run getUMAP if there are no reducedDimNames
      # redDimName <- input$qcPlotRedDim
      # show the tabs for the result plots  output[[qc_plot_ids[[a]]]]

      showQCResTabs(vals, algoList, qc_algo_status, qc_plot_ids)
      arrangeQCPlots(vals$counts, input, output, algoList,
                     colData(vals$counts)[[input$qcSampleSelect]], qc_plot_ids,
                     qc_algo_status, input$QCUMAPName)

      uniqueSampleNames = unique(colData(vals$counts)[[input$qcSampleSelect]])
      for (algo in algoList) {
        qc_algo_status[[algo]] <- list(self="done")
        if (length(uniqueSampleNames) > 1) {
          for (s in uniqueSampleNames) {
            qc_algo_status[[algo]][[s]] = TRUE
          }
        }
      }
    }
  }

  observeEvent(input$runQC, withConsoleMsgRedirect({
    #withBusyIndicatorServer("runQC", {
      if (!qcInputExists()) {
        insertUI(
          selector = "#qcPageErrors",
          ui = wellPanel(id = "noSelected", tags$b("Please select at least one algorithm.", style = "color: red;"))
        )
      } else if (is.null(vals$counts)) {
        insertUI(
          selector = "#qcPageErrors",
          ui = wellPanel(id = "noSCE", tags$b("Please upload a sample first.", style = "color: red;"))
        )
      } else if (is.null(input$qcAssaySelect)) {
        insertUI(
          selector = "#qcPageErrors",
          ui = wellPanel(id = "noQCAssay", tags$b("Please select an assay.", style = "color: red;"))
        )
      } else {
        removeUI(
          selector = "#noSelected"
        )
        removeUI(
          selector = "#noSCE"
        )
        removeUI(
          selector = "#noQCAssay"
        )
        useAssay <- input$qcAssaySelect
        qcSample <- colData(vals$original)[,input$qcSampleSelect]
        if (length(qcSample)==1 && qcSample == "None") {
          qcSample <- NULL
        }
        qcCollName <- NULL

        if (input$QCMgeneSets != "none") {
          if (input$QCMgeneSets %in% c("he", "hs", "me", "ms")) {
            if (input$QCMgeneSets == "he") {
              # Import Human Mito Ensembl
              mgsRef <- "human"
              mgsId <- "ensembl"
            } else if (input$QCMgeneSets == "hs") {
              # Import Human Mito Symbol
              mgsRef <- "human"
              mgsId <- "symbol"
            } else if (input$QCMgeneSets == "me") {
              # Import Mouse Mito Ensembl
              mgsRef <- "mouse"
              mgsId <- "ensembl"
            } else if (input$QCMgeneSets == "ms") {
              # Import Mouse Mito Symbol
              mgsRef <- "mouse"
              mgsId <- "symbol"
            }
            vals$counts <- importMitoGeneSet(inSCE = vals$counts,
                                             reference = mgsRef,
                                             id = mgsId,
                                             by = "rownames",
                                             collectionName = "Mito")
            qcCollName <- "Mito"
          } else {
            qcCollName <- input$QCMgeneSets
          }
        }
        algoList = list()
        paramsList <- list()
        for (algo in qc_choice_list) {
          if (isTRUE(input[[algo]])) {
            algoList <- c(algoList, algo)
            # use the specific prep functions for decontX and doubletFinder
            if (algo == "decontX") {
              paramsList <- prepDecontXParams(paramsList)
              next
            }
            if (algo == "doubletFinder") {
              paramsList <- prepDoubletFinderParams(paramsList)
              next
            }
            # everything else can go through the rest of the loop
            inputIds <- qc_input_ids[[algo]]
            algoParams <- list()
            for (key in names(inputIds)) {
              if(typeof(inputIds[[key]]) == "list") {
                paramSubList <- list()
                for (key2 in names(inputIds[[key]])) {
                  paramSubList[[key2]] <- input[[inputIds[[key]][[key2]]]]
                }
                algoParams[[key]] = paramSubList
              } else {
                if (nzchar(input[[inputIds[[key]]]])) {
                  algoParams[[key]] = NULL
                } else {
                  algoParams[[key]] = input[[inputIds[[key]]]]
                }
              }
            }
            paramsList[[algo]] = algoParams
          }
        }
        # run selected cell QC algorithms
        vals$counts <- runCellQC(inSCE = vals$counts,
                                 algorithms = algoList,
                                 sample = qcSample,
                                 collectionName = qcCollName,
                                 useAssay = input$qcAssaySelect,
                                 paramsList = paramsList)
        updateColDataNames()
        # redDimList <- strsplit(reducedDimNames(vals$counts), " ")
        # run getUMAP if doublet/ambient RNA detection conducted
        if(length(intersect(c("scDblFinder", "cxds", "bcds",
             "cxds_bcds_hybrid", "decontX",
             "scrublet", "doubletFinder"), algoList))){
          message(paste0(date(), " ... Running 'UMAP'"))
          vals$counts <- getUMAP(inSCE = vals$counts,
                                 sample = qcSample,
                                 useAssay = input$qcAssaySelect,
                                 nNeighbors = input$UnNeighbors,
                                 nIterations = input$UnIterations,
                                 alpha = input$Ualpha,
                                 minDist = input$UminDist,
                                 spread = input$Uspread,
                                 initialDims = input$UinitialDims,
                                 reducedDimName = input$QCUMAPName
                                 )
        }
        updateQCPlots()

        # Show downstream analysis options
        callModule(module = nonLinearWorkflow, id = "nlw-qcf", parent = session, nbc = TRUE, cw = TRUE, cv = TRUE)
      }
    #})

  }))

  #-----------#
  # FILTERING #####
  #-----------#
  shinyjs::onclick("colGT", shinyjs::toggle(id = "filterThreshGT",
                                            anim = FALSE), add = TRUE)
  shinyjs::onclick("colLT", shinyjs::toggle(id = "filterThreshLT",
                                            anim = FALSE), add = TRUE)
  filteringParams <- reactiveValues(params = list(), id_count = 0)
  rowFilteringParams <- reactiveValues(params = list(), id_count = 0)

  observeEvent(input$addFilteringParam, {
    if (!is.null(vals$counts)) {
      showModal(filteringModal(colNames = names(colData(vals$counts))))
    }
  })

  observeEvent(input$addRowFilteringParam, {
    if (!is.null(vals$counts) &&
        !is.null(names(assays(vals$counts)))) {
      showModal(rowFilteringModal(assayInput = names(assays(vals$counts))))
    }
  })

  observeEvent(input$filterColSelect, {
    # prep the modal - remove the threshold div and hide the categorical option
    shinyjs::hide("convertFilterType")
    removeUI(selector = "#newThresh")
    removeUI(selector = "div:has(>> #convertToCat)")
    # check if column contains numerical values
    isNum <- is.numeric(vals$counts[[input$filterColSelect]][0])
    if (length(vals$counts[[input$filterColSelect]]) > 0) {
      if (isTRUE(isNum)) {
        # (from partials) insertUI for choosing greater than and less than params
        addFilteringThresholdOptions(vals$counts[[input$filterColSelect]])
        # if less than 25 unique categories, give categorical option
        if (length(unique(vals$counts[[input$filterColSelect]])) < 25) {
          insertUI(
            selector = "#convertFilterType",
            ui = checkboxInput("convertToCat", "Convert to categorical filter?")
          )
          shinyjs::show("convertFilterType")
        }

      } else { # if non-numerical values, create checkbox input
        insertUI(
          selector = "#filterCriteria",
          ui = tags$div(id="newThresh",
                        checkboxGroupInput("filterThresh", "Please select which columns to keep:",
                                           choices = as.vector(unique(vals$counts[[input$filterColSelect]])),
                        ),
          )
        )
      }
    } else { # if no values in column, show error
      insertUI(
        selector = "#filterCriteria",
        ui = tags$div(id="newThresh", tags$b("This column does not have any filtering criteria", style = "color: red;"))
      )
    }
  })

  observeEvent(input$convertToCat, {
    if (!is.null(input$filterColSelect)) {
      removeUI(selector = "#newThresh")
      if (input$convertToCat) {
        insertUI(
          selector = "#filterCriteria",
          ui = tags$div(id="newThresh",
                        checkboxGroupInput("filterThresh", "Please select which columns to keep:",
                                           choices = as.vector(unique(vals$counts[[input$filterColSelect]])),
                        )
          )
        )
      } else {
        addFilteringThresholdOptions(vals$counts[[input$filterColSelect]])
        if (length(unique(vals$counts[[input$filterColSelect]])) < 25) {
          shinyjs::show("convertFilterType")
        }
      }
    }
  })

  observeEvent(input$filterAssaySelect, {
    removeUI(selector = "#newThresh")
    insertUI(
      selector = "#rowFilterCriteria",
      ui = tags$div(id="newThresh",
                    numericInput("filterThreshX", "Keep features with this many counts:", 0),
                    numericInput("filterThreshY", "In at least this many cells:", 0),
      )
    )

  })

  observeEvent(input$filtModalOK, {
    if (is.null(input$filterThresh) && is.null(input$filterThreshGT) && is.null(input$filterThreshLT)) {
      showModal(filteringModal(failed=TRUE, colNames = names(colData(vals$counts))))
    } else {
      id <- paste0("filteringParam", filteringParams$id_count)
      # figure out which options the user selected
      criteriaGT <- NULL
      criteriaLT <- NULL
      categoricalCol = FALSE
      if (isTRUE(input$colGT)) {
        criteriaGT = input$filterThreshGT
      }
      if (isTRUE(input$colLT)) {
        criteriaLT = input$filterThreshLT
      }
      if (!is.null(input$filterThresh)) {
          categoricalCol = TRUE
      }
      if (isTRUE(input$colLT) && isTRUE(input$colGT)) {
        if (criteriaGT > criteriaLT) {
          insertUI(
            selector = "#filterCrErrors",
            ui = wellPanel(id = "voidRange",
                           tags$b("Please set a valid range.",
                                  style = "color: red;"))
          )
          return()
        }
      }
      # new row in parameters table
      addToColFilterParams(name = input$filterColSelect,
                           categorial = categoricalCol,
                           criteria = input$filterThresh,
                           criteriaGT = criteriaGT,
                           criteriaLT = criteriaLT,
                           id = id,
                           paramsReactive = filteringParams)
      threshStr <- ""
      if (isTRUE(categoricalCol)) {
        threshStr <- paste(input$filterThresh, collapse = ', ')
      } else {
        if (is.null(criteriaGT)) {
          threshStr <- sprintf("< %.5f", input$filterThreshLT)
        } else if (is.null(criteriaLT)) {
          threshStr <- sprintf("> %.5f", input$filterThreshGT)
        } else {
          threshStr <- sprintf("> %.5f & < %.5f", input$filterThreshGT, input$filterThreshLT)
        }
      }
      make3ColTableRow("#newFilteringParams", id, input$filterColSelect, threshStr)
      observeEvent(input[[paste0("remove", id)]],{
        removeUI(
          selector = paste0("#", id)
        )
        toRemove <- vector()
        for (entry in filteringParams$params) {
          if (entry$id == id) {
            toRemove <- c(toRemove, FALSE)
          } else {
            toRemove <- c(toRemove, TRUE)
          }
        }
        filteringParams$params <- filteringParams$params[toRemove]
      })
      removeModal()
    }
  })

  observeEvent(input$rowFiltModalOK, {
    if ((is.null(input$filterThreshX)) || (is.null(input$filterThreshY)) || (is.null(input$filterAssaySelect))) {
      showModal(rowFilteringModal(failed=TRUE, assayInput = names(assays(vals$counts))))
    } else {
      id <- paste0("rowFilteringParam", rowFilteringParams$id_count)
      # new row in parameters table
      threshStr <- sprintf("> %i counts in > %i cells", input$filterThreshX, input$filterThreshY)

      addToRowFilterParams(input$filterAssaySelect, input$filterThreshX, input$filterThreshY, id, rowFilteringParams)
      make3ColTableRow("#newRowFilteringParams", id, input$filterAssaySelect, threshStr)
      observeEvent(input[[paste0("remove", id)]],{
        removeUI(
          selector = paste0("#", id)
        )
        toRemove <- vector()
        for (entry in rowFilteringParams$params) {
          if (entry$id == id) {
            toRemove <- c(toRemove, FALSE)
          } else {
            toRemove <- c(toRemove, TRUE)
          }
        }
        rowFilteringParams$params <- rowFilteringParams$params[toRemove]
      })
      removeModal()
    }
  })

  observeEvent(input$clearAllFilters, {
    for (entry in filteringParams$params) {
      removeUI(selector = paste0("#", entry$id))
    }
    filteringParams$params <- list()
  })

  observeEvent(input$clearAllRowParams, {
    for (entry in rowFilteringParams$params) {
      removeUI(selector = paste0("#", entry$id))
    }
    rowFilteringParams$params <- list()
  })

  observeEvent(input$filterSCE, {
    withBusyIndicatorServer("filterSCE", {
      # handle column filtering (pull out the criteria strings first)
      colInput <- formatFilteringCriteria(filteringParams$params)
      if (length(colInput) > 0) {
        vals$counts <- subsetSCECols(vals$counts, colData = colInput)
      }

      # handle row filtering (enter information as rows first, then pull out criteria strings)
      vals$counts <- addRowFiltersToSCE(vals$counts, rowFilteringParams)
      rowInput <- formatFilteringCriteria(rowFilteringParams$params)
      if (length(rowInput) > 0) {
        temp <- subsetSCERows(vals$counts, rowData = rowInput, returnAsAltExp = FALSE)
        if (nrow(temp) == 0) {
          stop("This filter will clear all rows. Filter has not been applied.")
        } else {
          vals$counts <- temp
        }
      }
      shinyjs::show(id="filteringSummary")

      # Show downstream analysis options
      shinyjs::show(selector = ".nlw-qcf")
    })
  })

  #Render summary table
  output$beforeFiltering <- renderTable({
      req(vals$original)
      if ("Sample" %in% names(colData(vals$counts))) {
        sampleVar <- "Sample"
      } else if ("sample" %in% names(colData(vals$counts))) {
        sampleVar <- "sample"
      } else {
        sampleVar <- NULL
      }
      # Setting 'useAssay=NULL' assumes that the first assay is the one to count
      singleCellTK::summarizeSCE(inSCE = vals$original,
                                 useAssay = NULL,
                                 sampleVariableName = sampleVar)
  }, striped = TRUE, border = TRUE, align = "c", spacing = "l")

  output$afterFiltering <- renderTable({
      req(vals$counts)
      if ("Sample" %in% names(colData(vals$counts))) {
        sampleVar <- "Sample"
      } else if ("sample" %in% names(colData(vals$counts))) {
        sampleVar <- "sample"
      } else {
        sampleVar <- NULL
      }
      # Setting 'useAssay=NULL' assumes that the first assay is the one to count
      singleCellTK::summarizeSCE(inSCE = vals$counts,
                                 useAssay = NULL,
                                 sampleVariableName = sampleVar)
  }, striped = TRUE, border = TRUE, align = "c", spacing = "l")

  #Render summary table
  output$summarycontents <- renderTable({
      req(vals$counts)
      if ("Sample" %in% names(colData(vals$counts))) {
        sampleVar <- "Sample"
      } else if ("sample" %in% names(colData(vals$counts))) {
        sampleVar <- "sample"
      } else {
        sampleVar <- NULL
      }
      # Setting 'useAssay=NULL' assumes that the first assay is the one to count
      singleCellTK::summarizeSCE(inSCE = vals$counts,
                                 useAssay = NULL,
                                 sampleVariableName = sampleVar)
  }, striped = TRUE, border = TRUE, align = "c", spacing = "l")


  #Reset the data to the original uploaded dataset
  observeEvent(input$resetData, {
    if (is.null(vals$original)){
      shinyalert::shinyalert("Error!", "Upload data first.", type = "error")
    }
    else{
      vals$counts <- vals$original
      #updateSelectInput(session, "deletesamplelist",
      #                  choices = colnames(vals$counts))
      vals$gsvaRes <- NULL
      vals$vamRes <- NULL
      vals$vamResults <- NULL
      vals$gsvaResults <- NULL
      vals$vamScore <- NULL
      vals$gsvaScore <- NULL
      vals$gsvaLimma <- NULL
      vals$visplotobject <- NULL
      vals$enrichRes <- NULL
      vals$dimRedPlot <- NULL
      vals$dimRedPlot_geneExp <- NULL
      vals$dendrogram <- NULL
      vals$pcX <- NULL
      vals$pcY <- NULL
      vals$batchRes <- NULL
      #Refresh things for the clustering tab
      updateColDataNames()
      updateFeatureAnnots()
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
      isolate({
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
  })

  #Filter the selected samples
  #  observeEvent(input$runFilterSample, {
  #    withBusyIndicatorServer("runFilterSample", {
  #      filter <- colData(vals$counts)[, input$filteredSample] %in% input$filterSampleChoices
  #      vals$counts <- vals$counts[, filter]
  #      vals$gsvaRes <- NULL
  #      vals$enrichRes <- NULL
  #      vals$visplotobject <- NULL
  #      vals$gsvaLimma <- NULL
  #      vals$dimRedPlot <- NULL
  #      vals$dimRedPlot_geneExp <- NULL
  #      vals$dendrogram <- NULL
  #      vals$pcX <- NULL
  #      vals$pcY <- NULL
  #      updateNumSamples()
  #    })
  #  })

  #  observeEvent(input$filteredFeature, {
  #    output$filterFeatureOptions <- renderUI({
  #      if (input$filteredFeature != "none")({
  #        if (length(unique(rowData(vals$counts)[, input$filteredFeature])) < 100){
  #          L <- vector("list", 3)
  #          L[[1]] <- renderText("Select features to keep")
  #          L[[2]] <- wellPanel(style = "overflow-y:scroll; max-height: 100px",
  #                              list(checkboxGroupInput("filterFeatureChoices",
  #                                                      label = NULL,
  #                                                      choices = unique(rowData(vals$counts)[, input$filteredFeature]))))
  #          L[[3]] <- list(actionButton("runFilterFeature", "Filter"))
  #          return(L)
  #        } else {
  #          L <- list(renderText("Annotation must have fewer than 100 options"))
  #          return(L)
  #        }
  #      }) else {
  #        L <- list()
  #      }
  #    })
  #  })

  #  observeEvent(input$orgOrganism, {
  #    library(input$orgOrganism, character.only = TRUE)
  #    indb <- get(paste(input$orgOrganism))
  #    output$orgConvertColumns <- renderUI({
  #      tagList(
  #        selectInput("orgFromCol", "Select From Annotation:", columns(indb)),
  #        selectInput("orgToCol", "Select To Annotation:", columns(indb))
  #      )
  #    })
  #  })

  #  observeEvent(input$convertGenes, {
  #    req(vals$counts)
  #    withBusyIndicatorServer("convertGenes", {
  #      vals$counts <- convertGeneIDs(inSCE = vals$counts,
  #                                    inSymbol = input$orgFromCol,
  #                                    outSymbol = input$orgToCol,
  #                                    database = input$orgOrganism)
  #      updateGeneNames()
  #      vals$gsvaRes <- NULL
  #      vals$enrichRes <- NULL
  #      vals$visplotobject <- NULL
  #      vals$dimRedPlot <- NULL
  #      vals$dimRedPlot_geneExp <- NULL
  #      vals$dendrogram <- NULL
  #      vals$pcX <- NULL
  #      vals$pcY <- NULL
  #    })
  #  })

  #Filter the selected features
  #  observeEvent(input$runFilterFeature, {
  #    filter <- rowData(vals$counts)[, input$filteredFeature] %in% input$filterFeatureChoices
  #    vals$counts <- vals$counts[filter, ]
  #    updateGeneNames()
  #    vals$gsvaRes <- NULL
  #    vals$enrichRes <- NULL
  #    vals$visplotobject <- NULL
  #    vals$dimRedPlot <- NULL
  #    vals$dimRedPlot_geneExp <- NULL
  #    vals$dendrogram <- NULL
  #    vals$pcX <- NULL
  #    vals$pcY <- NULL
  #  })

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

  # Delete Data ####

  output$reducedDimsList <- renderUI({
    req(vals$counts)
    if (!is.null(vals$counts) &&
        length(names(reducedDims(vals$counts))) > 0){
      panel(heading = "ReducedDims",
            checkboxGroupInput(
              inputId = "checkboxRedDimToRemove",
              label = NULL,
              choices = names(reducedDims(vals$counts))
              )
            )
    }
  })

  output$assaysList <- renderUI({
    req(vals$counts)
    if (!is.null(vals$counts)){
      panel(heading = "Assays",
            checkboxGroupInput(
              inputId = "checkboxAssaysToRemove",
              label = NULL,
              choices = assayNames(vals$counts)
              )
            )
    }
  })

  output$rowDataList <- renderUI({
    req(vals$counts)
    if (!is.null(vals$counts)
        && length(colnames(rowData(vals$counts))) > 0){
      panel(heading = "Row Annotation",
            checkboxGroupInput(
              inputId = "checkboxRowDataToRemove",
              label = NULL,
              choices = colnames(rowData(vals$counts))
            )
      )
    }
  })

  output$colDataList <- renderUI({
    req(vals$counts)
    if (!is.null(vals$counts)
        && length(colnames(colData(vals$counts))) > 0){
      panel(heading = "Column Annotation",
            checkboxGroupInput(
              inputId = "checkboxColDataToRemove",
              label = NULL,
              choices = colnames(colData(vals$counts))
            )
      )
    }
  })

  output$altExpList <- renderUI({
    req(vals$counts)
    if (!is.null(vals$counts)
        && length(altExpNames(vals$counts)) > 0){
      panel(heading = "Subsets",
            checkboxGroupInput(
              inputId = "checkboxAltExpToRemove",
              label = NULL,
              choices = altExpNames(vals$counts)
            )
      )
    }
  })

  observeEvent(input$delRedDim, {
    req(vals$counts)
    if(length(input$checkboxAssaysToRemove) > 0){
      for(i in seq(input$checkboxAssaysToRemove)){
        expData(vals$counts, input$checkboxAssaysToRemove[i]) <- NULL
        vals$counts <- expDeleteDataTag(vals$counts, input$checkboxAssaysToRemove[i])
      }
    }
    if(length(input$checkboxRedDimToRemove) > 0){
      for(i in seq(input$checkboxRedDimToRemove)){
        reducedDim(vals$counts, input$checkboxRedDimToRemove[i]) <- NULL
      }
    }
    if(length(input$checkboxRowDataToRemove) > 0){
      for(i in seq(input$checkboxRowDataToRemove)){
        rowData(vals$counts)[[input$checkboxRowDataToRemove[i]]] <- NULL
      }
    }
    if(length(input$checkboxColDataToRemove) > 0){
      for(i in seq(input$checkboxColDataToRemove)){
        colData(vals$counts)[[input$checkboxColDataToRemove[i]]] <- NULL
      }
    }
    if(length(input$checkboxAltExpToRemove) > 0){
      for(i in seq(input$checkboxAltExpToRemove)){
        altExps(vals$counts)[[input$checkboxAltExpToRemove[i]]] <- NULL
      }
    }
    updateAssayInputs()
    updateReddimInputs()
    updateFeatureAnnots()
    updateColDataNames()
  })

  # Normalization ####

  observeEvent(input$customNormalizeAssayMethodSelect,{
    if(input$customNormalizeAssayMethodSelect == "LogNormalize"
       || input$customNormalizeAssayMethodSelect == "CLR"
       || input$customNormalizeAssayMethodSelect == "SCTransform"
       || input$customNormalizeAssayMethodSelect == "logNormCounts"){
      updateAwesomeCheckbox(
        session = session,
        inputId = "customNormalizeOptionsTransform",
        value = FALSE
      )
    }
  })

  output$normalizationDataTagUI <- renderUI({
    req(vals$counts)
    tag <- ""
    if(input$normalizeAssayMethodSelect != "custom"){
      if(input$normalizeAssayMethodSelect
         %in% c("LogNormalize", "SCTransform", "CLR", "logNormCounts")){
        tag <- "transformed"
      }
      else{
        tag <- "normalized"
      }
      if(input$normalizationScale){
        tag <- "scaled"
      }
    }
    else{
      if(input$customNormalizeOptionsNormalize){
        if(input$customNormalizeAssayMethodSelect
           %in% c("LogNormalize", "SCTransform", "CLR", "logNormCounts")){
          tag <- "transformed"
        }
        else{
          tag <- "normalized"
        }
      }
      if(input$customNormalizeOptionsTransform){
        tag <- "transformed"
      }
      if(input$customNormalizeOptionsScale){
        tag <- "scaled"
      }
    }
    return(tag)
  })

  output$normalizationNormalizeSelectedMethodUI <- renderUI({
    req(vals$counts)
    if(input$normalizeAssayMethodSelect != "custom"){
      h5(input$normalizeAssayMethodSelect)
    }
    else{
      NULL
    }

  })

  observeEvent(input$modifyAssay, {
    req(vals$counts)
    withBusyIndicatorServer("modifyAssay", {
      if (!(input$modifyAssaySelect %in% names(assays(vals$counts)))) {
        stop("Assay does not exist!")
      } else if (input$modifyAssayOutname == "") {
        stop("Assay name cannot be empty!")
      } else if (input$modifyAssayOutname %in% names(assays(vals$counts))) {
        stop("Assay name already exists! Use another assay name!")
      } else if(is.na(input$trimUpperValueAssay)
                || is.na(input$trimLowerValueAssay)){
        stop("Upper or lower trim value cannot be empty!")
      } else {
        checkedOptions <- c(input$customNormalizeOptionsNormalize,
                            input$customNormalizeOptionsTransform,
                            input$customNormalizeOptionsPsuedocounts,
                            input$customNormalizeOptionsScale,
                            input$customNormalizeOptionsTrim)

        if(!any(checkedOptions)){
          stop("Must select at least one option!")
        }

        #Setting initial parameters
        normalizeMethod <- NULL
        transformMethod <- NULL
        pseudocountsBefore <- NULL
        pseudocountsAfter <- NULL
        doScale <- input$customNormalizeOptionsScale
        trimOptions <- NULL

        if(input$customNormalizeOptionsNormalize)
          normalizeMethod <- input$customNormalizeAssayMethodSelect
        if(input$customNormalizeOptionsTransform)
          transformMethod <- input$customNormalizeTransformOptions
        if(input$customNormalizePseudoOptionsBefore)
          pseudocountsBefore <- input$customNormalizePseudoValueBefore
        if(input$customNormalizePseudoOptionsAfter)
          pseudocountsAfter <- input$customNormalizePseudoValueAfter
        if(input$customNormalizeOptionsTrim)
          trimOptions <- c(input$trimUpperValueAssay, input$trimLowerValueAssay)

        outAssayName <- input$modifyAssayOutname
        useAssay <- input$modifyAssaySelect

        args <- list(
          inSCE = vals$counts,
          useAssay = useAssay,
          outAssayName = outAssayName,
          normalizationMethod = normalizeMethod,
          scale = doScale,
          transformation = transformMethod,
          pseudocountsBeforeNorm = pseudocountsBefore,
          pseudocountsBeforeTransform = pseudocountsAfter,
          trim = trimOptions
        )

        vals$counts <- do.call("runNormalization", args)

        # Show downstream analysis options
        callModule(module = nonLinearWorkflow, id = "nlw-nbc", parent = session, dr = TRUE, fs = TRUE)
      }
    })
  })

  observeEvent(input$normalizeAssay, {
    req(vals$counts)
    withBusyIndicatorServer("normalizeAssay", {
      if(!(input$normalizeAssaySelect %in% expDataNames(vals$counts))){
        stop("Selected assay does not exist!")
      }
      else if(input$normalizeAssayOutname == ""){
        stop("Assay Name cannot be empty!")
      }
      else if(input$normalizeAssayOutname %in% expDataNames(vals$counts)){
        stop("Your selected Assay Name already exists! Try another Assay Name!")
      }
      else if(input$normalizeAssaySelect == ""){
        stop("Please select an assay before proceeding with normalization!")
      }
      else if(is.na(as.numeric(input$normalizationScaleFactor))){
        stop("Scaling factor must be a numeric non-empty value!")
      }
      else{
        #Setting initial parameters
        normalizeMethod <- input$normalizeAssayMethodSelect
        doScale <- input$normalizationScale
        trimOptions <- NULL
        scaleFactor <- input$normalizationScaleFactor

        if(doScale && input$normalizationTrim)
          trimOptions <- c(input$normalizationTrimUpper, input$normalizationTrimLower)

        outAssayName <- input$normalizeAssayOutname
        useAssay <- input$normalizeAssaySelect

        args <- list(
          inSCE = vals$counts,
          useAssay = useAssay,
          outAssayName = outAssayName,
          normalizationMethod = normalizeMethod,
          scale = doScale,
          seuratScaleFactor = scaleFactor,
          trim = trimOptions
        )

        vals$counts <- do.call("runNormalization", args)

        # Show downstream analysis options
        callModule(module = nonLinearWorkflow, id = "nlw-nbc", parent = session, dr = TRUE, fs = TRUE)
      }
    })
  })

  observeEvent(input$normalizeAssayMethodSelect, {
    if(input$normalizeAssayMethodSelect == "LogNormalize") {
      updateTextInput(session = session, inputId = "normalizeAssayOutname", value = "SeuratLogNormalize")
    } else if(input$normalizeAssayMethodSelect == "CLR"){
      updateTextInput(session = session, inputId = "normalizeAssayOutname", value = "SeuratCLR")
    } else if(input$normalizeAssayMethodSelect == "RC"){
      updateTextInput(session = session, inputId = "normalizeAssayOutname", value = "SeuratRC")
    } else if(input$normalizeAssayMethodSelect == "CPM"){
      updateTextInput(session = session, inputId = "normalizeAssayOutname", value = "ScaterCPMCounts")
    } else if(input$normalizeAssayMethodSelect == "logNormCounts"){
      updateTextInput(session = session, inputId = "normalizeAssayOutname", value = "ScaterLogNormCounts")
    } else if(input$normalizeAssayMethodSelect == "SCTransform"){
      updateTextInput(session = session, inputId = "normalizeAssayOutname", value = "SeuratSCTransform")
    }
  })


  #  output$colDataDataFrame <- DT::renderDataTable({
  #    if (!is.null(vals$counts)){
  #      data.frame(colData(vals$counts))
  #    }
  #  }, options = list(scrollX = TRUE, scrollY = "40vh", pageLength = 30))

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
        HTML(paste("clicking on selected radio button option modifies the selected condition's data in the backend."))
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

  output$annotModifyUIHelpText <- renderUI({
    if (!is.null(vals$counts)){
      if (input$annotModifyChoice != "none"){
        HTML(paste(tags$h5("Note: Clicking on selected radio button option modifies the selected condition's data in the backend.")))
      }
    }
  })

  #-----------------------------------------------------------------------------
  # Page 3: dimRed ####
  #-----------------------------------------------------------------------------

  output$dimRedNameUI <- renderUI({
      defaultText <- paste(input$dimRedAssaySelect, input$dimRedPlotMethod,
                           sep = '_')
    textInput('dimRedNameInput', "reducedDim Name:", defaultText)
  })

  output$dimRedNameUI_tsneUmap <- renderUI({
      defaultText <- paste(input$dimRedAssaySelect_tsneUmap, input$dimRedPlotMethod_tsneUmap,
                           sep = '_')
    textInput('dimRedNameInput_tsneUmap', "reducedDim Name:", defaultText)
  })

  observeEvent(input$updateHeatmap_dimRed, {
    req(vals$counts)
    if (!is.null(input$picker_dimheatmap_components_dimRed)) {
      if(vals$runDimred$dimRedAssaySelect %in% assayNames(vals$counts)){
        output$plot_heatmap_dimRed <- renderPlot({
          isolate({
            singleCellTK:::.plotHeatmapMulti(
              plots = vals$counts@metadata$seurat$heatmap_dimRed,
              components = input$picker_dimheatmap_components_dimRed,
              nCol = input$slider_dimheatmap_dimRed)
          })
        })
      }
      else if(vals$runDimred$dimRedAssaySelect %in% expDataNames(vals$counts)){
        output$plot_heatmap_dimRed <- renderPlot({
          isolate({
            singleCellTK:::.plotHeatmapMulti(
              plots = altExps(vals$counts)[[vals$runDimred$dimRedAssaySelect]]@metadata$seurat$heatmap_dimRed,
              components = input$picker_dimheatmap_components_dimRed,
              nCol = input$slider_dimheatmap_dimRed)
          })
        })
      }
    }
    session$sendCustomMessage("close_dropDownDimRedHeatmap", "")
  })

  observeEvent(input$closeDropDownDimRedHeatmap, {
    session$sendCustomMessage("close_dropDownDimRedHeatmap", "")
  })

  observeEvent(input$runDimred, {
    if (!is.null(vals$counts)){
      withBusyIndicatorServer("runDimred", {
        vals$runDimred$dimRedAssaySelect <- input$dimRedAssaySelect
        if (vals$runDimred$dimRedAssaySelect %in% altExpNames(vals$counts)) {
          dimRedUseAltExp <- vals$runDimred$dimRedAssaySelect
        } else {
          dimRedUseAltExp <- NULL
        }
        if (input$dimRedNameInput == ""){
          shinyalert::shinyalert("Error", "enter a reducedDim name", type = "error")
        } #check for named entered and if its a duplicate
        else if (!is.null(input$dimRedNameInput)){
          dimrednamesave <- gsub(" ", "_", input$dimRedNameInput)
          if (input$dimRedNameInput %in% reducedDimNames(vals$counts)){
            stop("A reducedDim with name '", dimrednamesave, "' is already stored in the object. Please specify a different name for this reducedDim.")
          } else {
            vals$counts <- runDimReduce(
              inSCE = vals$counts,
              useAssay = vals$runDimred$dimRedAssaySelect,
              useAltExp = dimRedUseAltExp,
              method = input$dimRedPlotMethod,
              nComponents = input$dimRedNumberDims,
              reducedDimName = dimrednamesave)
            #vals$counts <- runDimensionalityReduction(
            #  inSCE = vals$counts,
            #  useAssay = vals$runDimred$dimRedAssaySelect,
            #  reducedDimName = dimrednamesave,
            #  method = input$dimRedPlotMethod,
            #  nComponents = input$dimRedNumberDims
            #)
            updateReddimInputs()
            # Show downstream analysis options
            callModule(module = nonLinearWorkflow, id = "nlw-dr", parent = session, cl = TRUE, cv = TRUE)
          }
        }
      })
    }
    dimrednamesave <- gsub(" ", "_", input$dimRedNameInput)
    if(input$dimRedPlotMethod == "scaterPCA"){
      redDim <- reducedDim(vals$counts, dimrednamesave)
      new_pca <- CreateDimReducObject(
        embeddings = redDim,
        assay = "RNA",
        loadings = attr(redDim, "rotation"),
        stdev = as.numeric(attr(redDim, "percentVar")),
        key = "PC_")
    }

    removeTab(inputId = "dimRedPCAICA_plotTabset", target = "Component Plot")
    removeTab(inputId = "dimRedPCAICA_plotTabset", target = "Elbow Plot")
    removeTab(inputId = "dimRedPCAICA_plotTabset", target = "Heatmap Plot")
    removeTab(inputId = "dimRedPCAICA_plotTabset", target = "JackStraw Plot")

    shinyjs::show(selector = ".dimRedPCAICA_plotTabset_class")

    if(input$computeElbowPlot
       && input$dimRedPlotMethod != "seuratICA"){
      appendTab(
        inputId = "dimRedPCAICA_plotTabset",
        tabPanel(
          title = "Elbow Plot",
          panel(
            #heading = "Elbow Plot",
            plotlyOutput(outputId = "plotDimRed_elbow")
          )
        ), select = TRUE
      )
      if (input$dimRedPlotMethod == "seuratPCA"){
        withProgress(message = "Generating Elbow Plot", max = 1, value = 1, {
          if(vals$runDimred$dimRedAssaySelect %in% assayNames(vals$counts)){
            output$plotDimRed_elbow <- renderPlotly({
              seuratElbowPlot(inSCE = vals$counts, )
            })
          } else if(vals$runDimred$dimRedAssaySelect %in% expDataNames(vals$counts)){
            output$plotDimRed_elbow <- renderPlotly({
              seuratElbowPlot(inSCE = altExps(vals$counts)[[vals$runDimred$dimRedAssaySelect]])
            })
          }
        })
      } else {
        withProgress(message = "Generating Elbow Plot", max = 1, value = 1, {
          if(input$dimRedAssaySelect %in% assayNames(vals$counts)){
            output$plotDimRed_elbow <- renderPlotly({
              seuratElbowPlot(inSCE = vals$counts,
                              externalReduction = new_pca)
            })
          } else if(input$dimRedAssaySelect %in% expDataNames(vals$counts)){
            output$plotDimRed_elbow <- renderPlotly({
              seuratElbowPlot(inSCE = altExps(vals$counts)[[vals$runDimred$dimRedAssaySelect]],
                              externalReduction = new_pca)
            })
          }
        })
      }
    }

    if(input$computeHeatmapPlot){
      appendTab(
        inputId = "dimRedPCAICA_plotTabset",
        tabPanel(
          title = "Heatmap Plot",
          tags$script("Shiny.addCustomMessageHandler('close_dropDownDimRedHeatmap', function(x){
                  $('html').click();
                });"),
          panel(
            fluidRow(
              column(4, dropdown(
                fluidRow(actionBttn(inputId = "closeDropDownDimRedHeatmap", label = NULL, style = "simple", color = "danger", icon = icon("times"), size = "xs"), align = "right"),
                      selectizeInput(inputId = "picker_dimheatmap_components_dimRed",
                                     label = "Select principal components to plot:",
                                     choices = c(),
                                     multiple = TRUE),
                      numericInput(
                        inputId = "slider_dimheatmap_dimRed",
                        label = "Number of columns for the plot: ",
                        min = 1,
                        max = 4,
                        value = 2
                      ),
                  actionBttn(
                    inputId = "updateHeatmap_dimRed",
                    label = "Update",
                    style = "bordered",
                    color = "primary",
                    size = "sm"
                  ),
                inputId = "dropDownDimRedHeatmap",
                icon = icon("cog"),
                status = "primary",
                circle = FALSE,
                inline = TRUE
              )),
              column(6, fluidRow(h6("Heatmaps of the top features correlated with each component"), align = "center"))
            ),
            hr(),
            br(),
              shinyjqui::jqui_resizable(
                plotOutput(outputId = "plot_heatmap_dimRed"),
                options = list(maxWidth = 700)
              )
          )
        )
      )
      if (input$dimRedPlotMethod == "seuratPCA") {
        withProgress(message = "Generating Heatmaps", max = 1, value = 1, {
          if(input$dimRedAssaySelect %in% assayNames(vals$counts)){
            vals$counts@metadata$seurat$heatmap_dimRed <- singleCellTK::computeHeatmap(
              inSCE = vals$counts,
              useAssay = input$dimRedAssaySelect,
              dims = 1:input$dimRedNumberDims,
              nfeatures = input$dimRedNFeaturesHeatmap,
              reduction = "pca"
            )
            output$plot_heatmap_dimRed <- renderPlot({
              singleCellTK:::.plotHeatmapMulti(vals$counts@metadata$seurat$heatmap_dimRed)
            })
          }
          else if(vals$runDimred$dimRedAssaySelect %in% expDataNames(vals$counts)){
            altExps(vals$counts)[[vals$runDimred$dimRedAssaySelect]]@metadata$seurat$heatmap_dimRed <- singleCellTK::computeHeatmap(
              inSCE = altExps(vals$counts)[[vals$runDimred$dimRedAssaySelect]],
              useAssay = vals$runDimred$dimRedAssaySelect,
              dims = 1:input$dimRedNumberDims,
              nfeatures = input$dimRedNFeaturesHeatmap,
              reduction = "pca"
            )
            output$plot_heatmap_dimRed <- renderPlot({
              singleCellTK:::.plotHeatmapMulti(altExps(vals$counts)[[vals$runDimred$dimRedAssaySelect]]@metadata$seurat$heatmap_dimRed)
            })
          }
        })
      }
      else if(input$dimRedPlotMethod == "seuratICA"){
        withProgress(message = "Generating Heatmaps", max = 1, value = 1, {
          if(vals$runDimred$dimRedAssaySelect %in% assayNames(vals$counts)){
            vals$counts@metadata$seurat$heatmap_dimRed <- singleCellTK::computeHeatmap(
              inSCE = vals$counts,
              useAssay = input$dimRedAssaySelect,
              dims = 1:input$dimRedNumberDims,
              nfeatures = input$dimRedNFeaturesHeatmap,
              reduction = "ica"
            )
            output$plot_heatmap_dimRed <- renderPlot({
              singleCellTK:::.plotHeatmapMulti(vals$counts@metadata$seurat$heatmap_dimRed)
            })
          }
          else if(vals$runDimred$dimRedAssaySelect %in% expDataNames(vals$counts)){
            altExps(vals$counts)[[vals$runDimred$dimRedAssaySelect]]@metadata$seurat$heatmap_dimRed <- singleCellTK::computeHeatmap(
              inSCE = altExps(vals$counts)[[vals$runDimred$dimRedAssaySelect]],
              useAssay = vals$runDimred$dimRedAssaySelect,
              dims = 1:input$dimRedNumberDims,
              nfeatures = input$dimRedNFeaturesHeatmap,
              reduction = "ica"
            )
            output$plot_heatmap_dimRed <- renderPlot({
              singleCellTK:::.plotHeatmapMulti(altExps(vals$counts)[[vals$runDimred$dimRedAssaySelect]]@metadata$seurat$heatmap_dimRed)
            })
          }
        })
      }
      else{
        withProgress(message = "Generating Heatmaps", max = 1, value = 1, {
          if(input$dimRedAssaySelect %in% assayNames(vals$counts)){
            vals$counts@metadata$seurat$heatmap_dimRed <- singleCellTK::computeHeatmap(
              inSCE = vals$counts,
              useAssay = input$dimRedAssaySelect,
              dims = 1:input$dimRedNumberDims,
              nfeatures = input$dimRedNFeaturesHeatmap,
              externalReduction = new_pca
            )
            output$plot_heatmap_dimRed <- renderPlot({
              singleCellTK:::.plotHeatmapMulti(vals$counts@metadata$seurat$heatmap_dimRed)
            })
          }
          else if(input$dimRedAssaySelect %in% expDataNames(vals$counts)){
            altExps(vals$counts)[[input$dimRedAssaySelect]]@metadata$seurat$heatmap_dimRed <- singleCellTK::computeHeatmap(
              inSCE = altExps(vals$counts)[[input$dimRedAssaySelect]],
              useAssay = input$dimRedAssaySelect,
              dims = 1:input$dimRedNumberDims,
              nfeatures = input$dimRedNFeaturesHeatmap,
              externalReduction = new_pca
            )
            output$plot_heatmap_dimRed <- renderPlot({
              singleCellTK:::.plotHeatmapMulti(altExps(vals$counts)[[vals$runDimred$dimRedAssaySelect]]@metadata$seurat$heatmap_dimRed)
            })
          }
        })
      }

      if(input$dimRedPlotMethod == "seuratICA"){
        updateSelectizeInput(session = session, inputId = "picker_dimheatmap_components_dimRed", choices = rep(paste0("IC",seq(as.numeric(input$dimRedNumberDims)))))
      }
      else{
        updateSelectizeInput(session = session, inputId = "picker_dimheatmap_components_dimRed", choices = rep(paste0("PC",seq(as.numeric(input$dimRedNumberDims)))))
      }
    }

    appendTab(inputId = "dimRedPCAICA_plotTabset", tabPanel(title = "Component Plot",
                                                            panel(
                                                              tags$script("Shiny.addCustomMessageHandler('close_dropDownDimRedComponentPlot', function(x){$('html').click();});"),
                                                              fluidRow(
                                                                column(4, dropdown(
                                                                  fluidRow(
                                                                    column(12,
                                                                           fluidRow(actionBttn(inputId = "closeDropDownDimRedComponentPlot", label = NULL, style = "simple", color = "danger", icon = icon("times"), size = "xs"), align = "right"),
                                                                           selectizeInput(
                                                                             inputId = "plotDimRed_pca_selectRedDim",
                                                                             label = "Select reducedDim:",
                                                                             choices = reducedDimNames(vals$counts)
                                                                           ),
                                                                           numericInput(inputId = "plotDimRed_pca_dimX", label = "Select component for X-axis:", value = 1),
                                                                           numericInput(inputId = "plotDimRed_pca_dimY", label = "Select component for Y-axis:", value = 2),
                                                                           actionBttn(
                                                                             inputId = "updateRedDimPlot_pca",
                                                                             label = "Update",
                                                                             style = "bordered",
                                                                             color = "primary",
                                                                             size = "sm"
                                                                           )
                                                                    )
                                                                  ),
                                                                  inputId = "dropDownDimRedComponentPlot",
                                                                  icon = icon("cog"),
                                                                  status = "primary",
                                                                  circle = FALSE,
                                                                  inline = TRUE
                                                                )),
                                                                column(6, fluidRow(h6("Scatterplot of cells on selected components from a dimensionality reduction"), align = "center"))
                                                              ),
                                                              hr(),
                                                              br(),
                                                                  plotlyOutput(outputId = "plotDimRed_pca")
                                                            )
    ))

    withProgress(message = "Plotting PCA/ICA", max = 1, value = 1, {
        output$plotDimRed_pca <- renderPlotly({
          plotly::ggplotly(
            plotDimRed(
              inSCE = vals$counts,
              useReduction = dimrednamesave,
              xAxisLabel = paste0(input$dimRedPlotMethod, "_1"),
              yAxisLabel = paste0(input$dimRedPlotMethod, "_2"))
          )
        })
    })

        if(input$computeJackstrawPlot
           && input$dimRedPlotMethod != "seuratICA"){
          appendTab(inputId = "dimRedPCAICA_plotTabset", tabPanel(title = "JackStraw Plot",
                                                                  panel(heading = "JackStraw Plot",
                                                                        shinyjqui::jqui_resizable(plotOutput(outputId = "plot_jackstraw_dimRed"))
                                                                  )
          ))

          if (input$dimRedPlotMethod == "seuratPCA"){
            withProgress(message = "Generating JackStraw Plot", max = 1, value = 1, {
              if(vals$runDimred$dimRedAssaySelect %in% assayNames(vals$counts)){
                vals$counts <- seuratComputeJackStraw(inSCE = vals$counts,
                                                      useAssay = input$dimRedAssaySelect,
                                                      dims = input$dimRedNumberDims)
                output$plot_jackstraw_dimRed <- renderPlot({
                  seuratJackStrawPlot(inSCE = vals$counts, dims = input$dimRedNumberDims)
                })
              }
              else if(vals$runDimred$dimRedAssaySelect %in% expDataNames(vals$counts)){
                altExps(vals$counts)[[vals$runDimred$dimRedAssaySelect]] <- seuratComputeJackStraw(inSCE = altExps(vals$counts)[[vals$runDimred$dimRedAssaySelect]],
                                                      useAssay = vals$runDimred$dimRedAssaySelect,
                                                      dims = input$dimRedNumberDims)
                output$plot_jackstraw_dimRed <- renderPlot({
                  seuratJackStrawPlot(inSCE = altExps(vals$counts)[[vals$runDimred$dimRedAssaySelect]], dims = input$dimRedNumberDims)
                })
              }
            })
          }
          else{
            withProgress(message = "Generating JackStraw Plot", max = 1, value = 1, {
              if(input$dimRedAssaySelect %in% assayNames(vals$counts)){
                vals$counts <- seuratComputeJackStraw(inSCE = vals$counts,
                                                      useAssay = input$dimRedAssaySelect,
                                                      dims = input$dimRedNumberDims,
                                                      externalReduction = new_pca)
                output$plot_jackstraw_dimRed <- renderPlot({
                  seuratJackStrawPlot(inSCE = vals$counts,
                                      dims = input$dimRedNumberDims)
                })
              }
              else if(input$dimRedAssaySelect %in% expDataNames(vals$counts)){
                altExps(vals$counts)[[input$dimRedAssaySelect]] <- seuratComputeJackStraw(inSCE = altExps(vals$counts)[[input$dimRedAssaySelect]],
                                                      useAssay = input$dimRedAssaySelect,
                                                      dims = input$dimRedNumberDims,
                                                      externalReduction = new_pca)
                output$plot_jackstraw_dimRed <- renderPlot({
                  seuratJackStrawPlot(inSCE = altExps(vals$counts)[[vals$runDimred$dimRedAssaySelect]],
                                      dims = input$dimRedNumberDims)
                })
              }
            })
          }
        }
  })

  observeEvent(input$updateRedDimPlot_pca,{
    req(vals$counts)
      output$plotDimRed_pca <- renderPlotly({
        isolate({
          plotly::ggplotly(
            plotDimRed(
              inSCE = vals$counts,
              useReduction = input$plotDimRed_pca_selectRedDim,
              xDim = input$plotDimRed_pca_dimX,
              yDim = input$plotDimRed_pca_dimY,
              xAxisLabel = paste0(input$dimRedPlotMethod, "_", input$plotDimRed_pca_dimX),
              yAxisLabel = paste0(input$dimRedPlotMethod, "_", input$plotDimRed_pca_dimY))
          )
        })
      })
      session$sendCustomMessage("close_dropDownDimRedComponentPlot", "")
  })

  observeEvent(input$closeDropDownDimRedComponentPlot, {
    req(vals$counts)
    session$sendCustomMessage("close_dropDownDimRedComponentPlot", "")
  })

  observeEvent(input$dimRedAssaySelect_tsneUmap, {
    req(vals$counts)
    if (!is.null(input$dimRedAssaySelect_tsneUmap)) {
      if (input$dimRedAssaySelect_tsneUmap %in% reducedDimNames(vals$counts)) {
        shinyjs::disable("reductionMethodUMAPTSNEDimRed")
      } else {
        shinyjs::enable("reductionMethodUMAPTSNEDimRed")
      }
    } else {
      shinyjs::enable("reductionMethodUMAPTSNEDimRed")
    }
  })
  observeEvent(input$runDimred_tsneUmap, {
    if (!is.null(vals$counts)){
      withBusyIndicatorServer("runDimred_tsneUmap", {
        vals$runDimred$dimRedAssaySelect_tsneUmap <- input$dimRedAssaySelect_tsneUmap
        if (vals$runDimred$dimRedAssaySelect_tsneUmap %in% reducedDimNames(vals$counts)) {
          embedUseAssay <- NULL
          embedUseRedDim <- vals$runDimred$dimRedAssaySelect_tsneUmap
          embedUseAltExp <- NULL
        } else if (vals$runDimred$dimRedAssaySelect_tsneUmap %in% altExpNames(vals$counts)) {
          embedUseAssay <- vals$runDimred$dimRedAssaySelect_tsneUmap
          embedUseRedDim <- NULL
          embedUseAltExp <- vals$runDimred$dimRedAssaySelect_tsneUmap
        } else if (vals$runDimred$dimRedAssaySelect_tsneUmap %in% assayNames(vals$counts)) {
          embedUseAssay <- vals$runDimred$dimRedAssaySelect_tsneUmap
          embedUseRedDim <- NULL
          embedUseAltExp <- NULL
        }
        if (input$dimRedNameInput_tsneUmap == ""){
          shinyalert::shinyalert("Error", "enter a reducedDim name", type = "error")
        } #check for named entered and if its a duplicate
        else if (!is.null(input$dimRedNameInput_tsneUmap)){
          if (input$dimRedNameInput_tsneUmap %in% names(reducedDims(vals$counts))){
            stop("A reducedDim with name '", input$dimRedNameInput_tsneUmap, "' is already stored in the object. Please specify a different name for this reducedDim.")
          } else {
            dimrednamesave <- gsub(" ", "_", input$dimRedNameInput_tsneUmap)
            if (input$dimRedPlotMethod_tsneUmap == "rTSNE"){
              vals$counts <- runDimReduce(
                inSCE = vals$counts,
                useAssay = embedUseAssay,
                useReducedDim = embedUseRedDim,
                useAltExp = embedUseAltExp,
                method = "rTSNE",
                reducedDimName = dimrednamesave,
                perplexity = input$perplexityTSNE,
                nIterations = input$iterTSNE
              )
              #vals$counts <- runDimensionalityReduction(
              #  inSCE = vals$counts,
              #  useAssay = input$dimRedAssaySelect_tsneUmap,
              #  reducedDimName = dimrednamesave,
              #  method = input$dimRedPlotMethod_tsneUmap,
              #  perplexity = input$perplexityTSNE,
              #  nIterations = input$iterTSNE
              #)
            } else if(input$dimRedPlotMethod_tsneUmap == "seuratTSNE"){
              if (!is.null(embedUseRedDim)) {
                vals$counts <- runDimReduce(
                  inSCE = vals$counts,
                  useAssay = embedUseAssay,
                  useReducedDim = embedUseRedDim,
                  useAltExp = embedUseAltExp,
                  method = "seuratTSNE",
                  reducedDimName = dimrednamesave,
                  dims = input$dimRedNumberDims_tsneUmap,
                  perplexity = input$perplexityTSNE
                )
              } else {
                vals$counts <- runDimReduce(
                  inSCE = vals$counts,
                  useAssay = embedUseAssay,
                  useReducedDim = embedUseRedDim,
                  useAltExp = embedUseAltExp,
                  method = "seuratTSNE",
                  reducedDimName = dimrednamesave,
                  dims = input$dimRedNumberDims_tsneUmap,
                  perplexity = input$perplexityTSNE,
                  useReduction = input$reductionMethodUMAPTSNEDimRed
                )
              }
              #vals$counts <- runDimensionalityReduction(
              #  inSCE = vals$counts,
              #  useAssay = input$dimRedAssaySelect_tsneUmap,
              #  reducedDimName = dimrednamesave,
              #  method = input$dimRedPlotMethod_tsneUmap,
              #  nComponents = input$dimRedNumberDims_tsneUmap,
              #  perplexity = input$perplexityTSNE,
              #  useReduction = input$reductionMethodUMAPTSNEDimRed
              #)
            } else if(input$dimRedPlotMethod_tsneUmap == "seuratUMAP"){
              if (!is.null(embedUseRedDim)) {
                vals$counts <- runDimReduce(
                  inSCE = vals$counts,
                  useAssay = embedUseAssay,
                  useReducedDim = embedUseRedDim,
                  useAltExp = embedUseAltExp,
                  method = "seuratUMAP",
                  reducedDimName = dimrednamesave,
                  dims = input$dimRedNumberDims_tsneUmap,
                  minDist = input$minDistUMAPDimRed,
                  nNeighbors = input$nNeighboursUMAPDimRed,
                  spread = input$spreadUMAPDimRed
                )
              } else {
                vals$counts <- runDimReduce(
                  inSCE = vals$counts,
                  useAssay = embedUseAssay,
                  useReducedDim = embedUseRedDim,
                  useAltExp = embedUseAltExp,
                  method = "seuratUMAP",
                  reducedDimName = dimrednamesave,
                  dims = input$dimRedNumberDims_tsneUmap,
                  minDist = input$minDistUMAPDimRed,
                  nNeighbors = input$nNeighboursUMAPDimRed,
                  spread = input$spreadUMAPDimRed,
                  useReduction = input$reductionMethodUMAPTSNEDimRed
                )
              }
              #vals$counts <- runDimensionalityReduction(
              #  inSCE = vals$counts,
              #  useAssay = input$dimRedAssaySelect_tsneUmap,
              #  reducedDimName = dimrednamesave,
              #  method = input$dimRedPlotMethod_tsneUmap,
              #  nComponents = input$dimRedNumberDims_tsneUmap,
              #  minDist = input$minDistUMAPDimRed,
              #  nNeighbors = input$nNeighboursUMAPDimRed,
              #  spread = input$spreadUMAPDimRed,
              #  useReduction = input$reductionMethodUMAPTSNEDimRed
              #)
            }
            else {
              if (is.na(input$alphaUMAP)) {
                stop("Learning rate (alpha) must be a numeric non-empty value!")
              }
              vals$counts <- runDimReduce(
                inSCE = vals$counts,
                useAssay = embedUseAssay,
                useReducedDim = embedUseRedDim,
                useAltExp = embedUseAltExp,
                method = "scaterUMAP",
                reducedDimName = dimrednamesave,
                nNeighbors = input$neighborsUMAP,
                nIterations = input$iterUMAP,
                minDist = input$mindistUMAP,
                alpha = input$alphaUMAP,
                spread = input$spreadUMAP
              )
              #vals$counts <- runDimensionalityReduction(
              #  inSCE = vals$counts,
              #  useAssay = input$dimRedAssaySelect_tsneUmap,
              #  reducedDimName = dimrednamesave,
              #  method = input$dimRedPlotMethod_tsneUmap,
              #  nNeighbors = input$neighborsUMAP,
              #  nIterations = input$iterUMAP,
              #  minDist = input$mindistUMAP,
              #  alpha = input$alphaUMAP
              #)
            }
            updateReddimInputs()
            # Show downstream analysis options
            callModule(module = nonLinearWorkflow, id = "nlw-dr", parent = session, cl = TRUE, cv = TRUE)
          }
        }
      })
    }

    redDimName <- gsub(" ", "_", input$dimRedNameInput_tsneUmap)

    updateSelectizeInput(session, "selectRedDimPlot_tsneUmap",
                         choices = c(reducedDimNames(vals$counts)),
                         selected = redDimName,
                         server = TRUE)

    withProgress(message = "Plotting tSNE/UMAP", max = 1, value = 1, {
        output$plotDimRed_tsneUmap <- renderPlotly({
          plotly::ggplotly(plotDimRed(
            inSCE = vals$counts,
            useReduction = redDimName,
            xAxisLabel = paste0(input$dimRedPlotMethod_tsneUmap,"_1"),
            yAxisLabel = paste0(input$dimRedPlotMethod_tsneUmap,"_2")
          ))
        })
    })
  })

  observeEvent(input$updateRedDimPlot_tsneUmap,{
    req(vals$counts)
    output$plotDimRed_tsneUmap <- renderPlotly({
      isolate({
        plotly::ggplotly(plotDimRed(
          inSCE = vals$counts,
          useReduction = input$selectRedDimPlot_tsneUmap,
          xAxisLabel = paste0(input$selectRedDimPlot_tsneUmap,"_1"),
          yAxisLabel = paste0(input$selectRedDimPlot_tsneUmap,"_2")
        ))
      })
    })

    session$sendCustomMessage("close_dropDownDimRedEmbedding", "")
  })

  observeEvent(input$closeDropDownDimRedEmbedding,{
    req(vals$counts)
    session$sendCustomMessage("close_dropDownDimRedEmbedding", "")
  })

  #-----------------------------------------------------------------------------
  # Page 3: Clustering ####
  #-----------------------------------------------------------------------------

  output$clustNameUI <- renderUI({
    if(input$clustAlgo %in% seq(6)){
      # Scran SNN
      textInput("clustName", "Name of Clustering Result:",
                "scran_snn_cluster")
    } else if(input$clustAlgo %in% seq(7, 9)){
      # K-Means
      textInput("clustName", "Name of Clustering Result:",
                "kmeans_cluster")
    } else if(input$clustAlgo %in% seq(10, 12)){
      algoList <- list('10' = "louvain",
                       '11' = "multilevel", '12' = "SLM")
      algo <- algoList[[as.character(input$clustAlgo)]]
      disabled(
        textInput("clustName", "Name of Clustering Result:",
                  paste0("Seurat", "_", algo, "_",
                         "Resolution", input$clustSeuratRes))
      )
    }
  })

  clustResults <- reactiveValues(names = NULL)

  getTypeByMat <- function(inSCE, matName) {
    if (matName %in% assayNames(inSCE)) {
      return("assay")
    } else if (matName %in% altExpNames(inSCE)) {
      return("altExp")
    } else if (matName %in% reducedDimNames(inSCE)) {
      return("reducedDim")
    } else {
      for (i in altExpNames(inSCE)) {
        if (matName %in% reducedDimNames(altExp(inSCE, i))) {
          return(c("reducedDim", i))
        }
      }
      return()
    }
  }

  observeEvent(input$clustRun, {
    if (is.null(vals$counts)){
      shinyalert::shinyalert("Error!", "Upload data first.", type = "error")
    } else if (input$clustName == "") {
      shinyalert::shinyalert("Error!", "Cluster name required.", type = "error")
    } else {
      withBusyIndicatorServer("clustRun", {
        saveClusterName = gsub(" ", "_", input$clustName)
        if(input$clustAlgo %in% seq(6)){
          # Scran SNN
          if(is.na(input$clustScranSNNK)){
            stop("K must be a numeric non-empty value!")
          }
          if(is.na(input$clustScranSNNd)){
            stop("Number of components must be a numeric non-empty value!")
          }
          algoList <- list('1' = "walktrap", '2' = "louvain", '3' = "infomap",
                           '4' = "fastGreedy", '5' = "labelProp",
                           '6' = "leadingEigen")
          algo <- algoList[[as.character(input$clustAlgo)]]
          params = list(inSCE = vals$counts,
                        clusterName = saveClusterName,
                        k = input$clustScranSNNK,
                        weightType = input$clustScranSNNType,
                        algorithm = algo)
          matType <- getTypeByMat(vals$counts, input$clustScranSNNMat)
          if (is.null(matType)) {
            return()
          } else if (length(matType) == 1) {
            if (matType == "assay") {
              params$useAssay = input$clustScranSNNMat
              params$nComp = input$clustScranSNNd
              plotReddim <- NULL
            } else if (matType == "reducedDim") {
              params$useReducedDim = input$clustScranSNNMat
              updateSelectInput(session, "clustVisReddim",
                                selected = input$clustScranSNNMat)
              plotReddim <- input$clustScranSNNMat
            } else if (matType == "altExp") {
              params$useAltExp = input$clustScranSNNMat
              params$altExpAssay = input$clustScranSNNMat
              params$nComp = input$clustScranSNNd
              plotReddim <- NULL
            }
          } else if (length(matType) == 2 &&
                     matType[1] == "reducedDim") {
            # Using reddims saved in altExp
            params$useAltExp = matType[2]
            params$altExpRedDim = input$clustScranSNNMat
            updateSelectInput(session, "clustVisReddim",
                              selected = input$clustScranSNNMat)
          }
          vals$counts <- do.call(runScranSNN, params)
        } else if (input$clustAlgo %in% seq(7, 9)) {
          # K-Means
          if(input$clustKMeansReddim == ""){
            stop("Must select a reducedDim! If none available, compute one in the Dimensionality Reduction tab.")
          }
          if(is.na(input$clustKMeansN)){
            stop("Number of clusters/centers must be a numeric non-empty value!")
          }
          if(is.na(input$clustKMeansNIter)){
            stop("Max number of iterations must be a numeric non-empty value!")
          }
          if(is.na(input$clustKMeansNStart)){
            stop("Number of random sets must be a numeric non-empty value!")
          }
          algoList <- list('7' = "Hartigan-Wong",
                           '8' = "Lloyd", '9' = "MacQueen")
          algo <- algoList[[as.character(input$clustAlgo)]]
          vals$counts <- runKMeans(inSCE = vals$counts,
                                   useReducedDim = input$clustKMeansReddim,
                                   nCenters = input$clustKMeansN,
                                   nIter = input$clustKMeansNIter,
                                   nStart = input$clustKMeansNStart,
                                   algorithm = algo,
                                   clusterName = saveClusterName)
          updateSelectInput(session, "clustVisReddim",
                            selected = input$clustKMeansReddim)
          plotReddim <- input$clustKMeansReddim
        } else if (input$clustAlgo %in% seq(10, 12)) {
          # Seurat
          if(input$clustSeuratReddim == ""){
            stop("Must select a reducedDim! If none available, compute one in the Dimensionality Reduction tab.")
          }
          if(is.na(input$clustSeuratDims)){
            stop("Number of dimensions must be a numeric non-empty value!")
          }
          if(is.na(input$clustSeuratRes)){
            stop("Resolution must be a numeric non-empty value!")
          }
          reddim <- reducedDim(vals$counts, input$clustSeuratReddim)
          rownames(reddim) <- gsub("_", "-", rownames(reddim))
          if ("percentVar" %in% names(attributes(reddim))) {
            stdev <- as.numeric(attr(reddim, "percentVar"))
            new_pca <- CreateDimReducObject(embeddings = reddim, assay = "RNA",
                                            stdev = stdev, key = "PC_")
          } else {
            new_pca <- CreateDimReducObject(embeddings = reddim, assay = "RNA",
                                            key = "PC_")
          }
          if (input$clustSeuratDims > ncol(reddim)) {
            warning("More dimensions specified in dims than have been computed")
            dims <- ncol(reddim)
          } else {
            dims <- input$clustSeuratDims
          }
          useAssay <- assayNames(vals$counts)[1]
          algoList <- list('10' = "louvain",
                           '11' = "multilevel", '12' = "SLM")
          algo <- algoList[[as.character(input$clustAlgo)]]
          vals$counts <- seuratFindClusters(inSCE = vals$counts,
                                            useAssay = useAssay,
                                            useReduction = "pca",
                                            externalReduction = new_pca,
                                            dims = dims,
                                            algorithm = algo,
                                            groupSingletons = input$clustSeuratGrpSgltn,
                                            resolution = input$clustSeuratRes)
          updateSelectInput(session, "clustVisReddim",
                            selected = input$clustSeuratReddim)
          plotReddim <- input$clustSeuratReddim
        }
        updateColDataNames()
        clustResults$names <- c(clustResults$names, saveClusterName)
        updateSelectInput(session, "clustVisRes", choices = clustResults$names)
        if (!is.null(plotReddim)) {
          output$clustVisPlot <- renderPlotly({
            isolate({
              plotSCEDimReduceColData(inSCE = vals$counts,
                                      colorBy = saveClusterName,
                                      conditionClass = "factor",
                                      reducedDimName = plotReddim,
                                      labelClusters = TRUE,
                                      dim1 = 1, dim2 = 2,
                                      legendTitle = saveClusterName)
            })
          })
        }
        # Show downstream analysis options
        callModule(module = nonLinearWorkflow, id = "nlw-cl", parent = session, de = TRUE, pa = TRUE, cv = TRUE)
      })
    }
  })
  
  observeEvent(input$closeDropDownClust, {
    session$sendCustomMessage("close_dropDownClust", "")
  })

  observeEvent(input$clustPlot, {
    if (is.null(vals$counts)){
      shinyalert::shinyalert("Error!", "Upload data first.", type = "error")
    } else {
      choice <- NULL
      if (input$clustVisChoicesType == 1) {
        # Use result
        if (is.null(input$clustVisRes) ||
            input$clustVisRes == "") {
          shinyalert::shinyalert("Error!", "Select the clusters to plot",
                                 type = "error")
        }
        choice <- input$clustVisRes
      } else if (input$clustVisChoicesType == 2) {
        # Use colData
        if (is.null(input$clustVisCol) ||
            input$clustVisCol == "") {
          shinyalert::shinyalert("Error!", "Select the clusters to plot",
                                 type = "error")
        }
        choice <- input$clustVisCol
      }
      if (is.null(input$clustVisReddim) || input$clustVisReddim == "") {
        shinyalert::shinyalert("Error!",
          "No reduction selected. Select one or run dimension reduction first",
          type = "error")
      }
      if (!is.null(choice) && choice != "" &&
          !is.null(input$clustVisReddim) && input$clustVisReddim != "") {
        output$clustVisPlot <- renderPlotly({
          isolate({
            plotSCEDimReduceColData(inSCE = vals$counts,
                                    colorBy = choice,
                                    conditionClass = "factor",
                                    reducedDimName = input$clustVisReddim,
                                    labelClusters = TRUE,
                                    dim1 = 1, dim2 = 2,
                                    legendTitle = choice)
          })
        })
      }
      session$sendCustomMessage("close_dropDownClust", "")
    }
  })

  #-----------------------------------------------------------------------------
  # Page 3.2: Celda ####
  #-----------------------------------------------------------------------------

  observeEvent(input$navbar, {
    if(!is.null(vals$counts)){
      if(input$navbar == "CeldaWorkflow"){
        updateSelectInput(session, "celdaassayselect", choices = c(names(assays(vals$counts))))
      }
    }
  })

  modsplit <- reactiveVal()
  cellsplit <- reactiveVal(NULL)

  observeEvent(input$celdamodsplit, {
    removeTab(inputId = "celdaModsplitTabset", target = "Perplexity Plot")
    removeTab(inputId = "celdaModsplitTabset", target = "Perplexity Difference Plot")
    appendTab(inputId = "celdaModsplitTabset", tabPanel(title = "Rate of perplexity change",
                                                        panel(heading = "RPC Plot",
                                                              plotlyOutput(outputId = "plot_modsplit_perpdiff", height = "auto")
                                                        )
    ), select = TRUE)
    appendTab(inputId = "celdaModsplitTabset", tabPanel(title = "Perplexity Plot",
      panel(heading = "Perplexity Plot",
        plotlyOutput(outputId = "plot_modsplit_perp", height = "auto")
      )
    ))
    
    withBusyIndicatorServer("celdamodsplit",{
      if (input$celdafeatureselect == "None"){
        vals$counts <- selectFeatures(vals$counts, minCount = input$celdarowcountsmin,
                                      minCell = input$celdacolcountsmin, useAssay = input$celdaassayselect)
      }else if(input$celdafeatureselect == "SeuratFindHVG"){
        vals$counts <- seuratNormalizeData(vals$counts, useAssay = input$celdaassayselect)
        vals$counts <- seuratFindHVG(vals$counts, useAssay = "seuratNormData",
                                     hvgMethod = input$celdaseurathvgmethod, hvgNumber = input$celdafeaturenum)
        
        g <- getTopHVG(vals$counts, method = input$celdaseurathvgmethod, n = input$celdafeaturenum)
        altExp(vals$counts, "featureSubset") <- vals$counts[g, ]
        
        vals$counts <- selectFeatures(vals$counts[g, ], minCount = input$celdarowcountsmin,
                                     minCell = input$celdacolcountsmin, useAssay = input$celdaassayselect, altExpName = "featureSubset")
      }else if(input$celdafeatureselect == "Scran_modelGeneVar"){
        if (!("ScaterLogNormCounts" %in% names(assays(vals$counts)))){
          vals$counts <- scater::logNormCounts(vals$counts, name = "ScaterLogNormCounts",
                                               exprs_values = input$celdaassayselect)
        }
        vals$counts <- scranModelGeneVar(vals$counts, assayName = "ScaterLogNormCounts")
        g <- getTopHVG(vals$counts, method = "modelGeneVar", n = input$celdafeaturenum)
        altExp(vals$counts, "featureSubset") <- vals$counts[g, ]
        
        vals$counts <- selectFeatures(vals$counts[g, ], minCount = input$celdarowcountsmin,
                                      minCell = input$celdacolcountsmin, useAssay = input$celdaassayselect, altExpName = "featureSubset")
      }
      counts(altExp(vals$counts)) <- as.matrix(counts(altExp(vals$counts)))
      updateNumericInput(session, "celdaLselect", min = input$celdaLinit, max = input$celdaLmax, value = input$celdaLinit)
      modsplit(recursiveSplitModule(vals$counts, useAssay = input$celdaassayselect, altExpName = "featureSubset",  initialL = input$celdaLinit, maxL = input$celdaLmax))
      output$plot_modsplit_perpdiff <- renderPlotly({plotRPC(modsplit(), sep = 10)})
      output$plot_modsplit_perp <- renderPlotly({plotGridSearchPerplexity(modsplit())})
      
    })

    shinyjs::enable(
      selector = ".celda_modsplit_plots a[data-value='Perplexity Plot']")
    shinyjs::enable(
      selector = ".celda_modsplit_plots a[data-value='Perplexity Diff Plot']")
    shinyjs::show(selector = ".celda_modsplit_plots")
    showNotification("Module splitting complete.")
    shinyjs::show(id = "celdaLselect")
    shinyjs::show(id = "celdaLbtn")
  })

  observeEvent(input$celdaLbtn, {
    vals$counts <- subsetCeldaList(modsplit(), params = list(L = input$celdaLselect))
    showNotification("Number of Feature Modules Selected.")
    updateCollapse(session = session, "CeldaUI", style = list("Identify Number of Feature Modules" = "success"))
    shinyjs::enable(selector = "div[value='Identify Number of Cell Clusters']")
  })

  output$celdaKplots <- renderUI({
    if (!is.null(vals$counts)){
      if (!is.null(cellsplit())){
        clusterlist <- runParams(cellsplit())$K
        plot_output_list <- lapply(runParams(cellsplit())$K, function(i){
          plotname <- paste0("Cluster", i)
          tabPanel(title = sprintf("Cluster %s", i),
                   panel(heading = sprintf("Cluster %s", i),
                    plotlyOutput(plotname)
                   )
          )
        })
        myTabs <- lapply(clusterlist, tabPanel)
        do.call(tabsetPanel, plot_output_list)
      }
    }
  })

  observeEvent(input$celdacellsplit, {
    withBusyIndicatorServer("celdacellsplit", {
      cellsplit(recursiveSplitCell(vals$counts, useAssay = input$celdaassayselect, initialK = input$celdaKinit, maxK = input$celdaKmax,
                                        yInit = celdaModules(vals$counts)))
      temp_umap <- celdaUmap(vals$counts)
      output$plot_cellsplit_perpdiff <- renderPlotly({plotRPC(cellsplit(), sep = 10)})
      output$plot_cellsplit_perp <- renderPlotly({plotGridSearchPerplexity(cellsplit())})
      
      for (i in runParams(cellsplit())$K){
        local({
          my_i <- i
          plotname <- paste0("Cluster", my_i)
          celdamod <- subsetCeldaList(cellsplit(), params = list(K = my_i))
          output[[plotname]] <- renderPlotly(plotDimReduceCluster(celdamod,
                                                                  dim1= reducedDim(altExp(temp_umap), "celda_UMAP")[, 1],
                                                                  dim2 = reducedDim(altExp(temp_umap), "celda_UMAP")[, 2],
                                                                  labelClusters = TRUE))
        })
      }
    })
    shinyjs::show(selector = ".celda_cellsplit_plots")
    showNotification("Cell Clustering Complete.")
    updateNumericInput(session, "celdaKselect", min = input$celdaKinit, max = input$celdaKmax, value = input$celdaKinit)
    shinyjs::show(id = "celdaKselect")
    shinyjs::show(id = "celdaKbtn")
  })

  observeEvent(input$celdaKbtn, {
    vals$counts <- subsetCeldaList(cellsplit(), params = list(K = input$celdaKselect))
    showNotification("Number of Cell Clusters Selected.")
    updateCollapse(session = session, "CeldaUI", style = list("Identify Number of Cell Clusters" = "success"))
    shinyjs::enable(
      selector = "div[value='Visualization']")
    updateNumericInput(session, "celdamodheatmapnum", min = 1, max = input$celdaLselect, value = 1)
    # Show downstream analysis options
    callModule(module = nonLinearWorkflow, id = "nlw-celda", parent = session, de = TRUE, pa = TRUE)
    
  })

  output$celdaheatmapplt <- renderPlot({plot(celdaHeatmap(vals$counts))})
  output$celdaprobmapplt <- renderPlot({celdaProbabilityMap(vals$counts)})

  observeEvent(input$CeldaUmap, {
    withBusyIndicatorServer("CeldaUmap", {
      vals$counts <- celdaUmap(vals$counts,
                               useAssay = input$celdaassayselect,
                               maxCells = input$celdaUMAPmaxCells,
                               minClusterSize = input$celdaUMAPminClusterSize,
                               seed = input$celdaUMAPSeed,
                               minDist = input$celdaUMAPmindist,
                               spread = input$celdaUMAPspread,
                               nNeighbors = input$celdaUMAPnn)
      output$celdaumapplot <- renderPlotly({plotDimReduceCluster(vals$counts, reducedDimName = "celda_UMAP", xlab = "UMAP_1",
                                                                 ylab = "UMAP_2", labelClusters = TRUE)})
      
    })
    showNotification("Umap complete.")
    colData(vals$counts)$celda_clusters <- celdaClusters(vals$counts)
    updateColDataNames()
    shinyjs::enable("CeldaTsne")
  })

  observeEvent(input$CeldaTsne, {
    withBusyIndicatorServer("CeldaTsne", {
      vals$counts <- celdaTsne(vals$counts,
                               useAssay = input$celdaassayselect,
                               maxCells = input$celdatSNEmaxCells,
                               minClusterSize = input$celdatSNEminClusterSize,
                               perplexity = input$celdatSNEPerplexity,
                               maxIter = input$celdatSNEmaxIter,
                               seed = input$celdatSNESeed)
      output$celdatsneplot <- renderPlotly({plotDimReduceCluster(vals$counts, reducedDimName = "celda_tSNE", xlab = "tSNE_1",
                                                                 ylab = "tSNE_2", labelClusters = TRUE)})
    })
    showNotification("Tsne complete.")
  })

  observeEvent(input$celdamodheatmapbtn,{
    output$celdamodheatmapplt <- renderPlot({moduleHeatmap(vals$counts, topCells= input$celdamodheatmaptopcells, featureModule = input$celdamodheatmapnum)})
    output$celdamodprobplt <- renderPlot({plotDimReduceModule(vals$counts, modules =  input$celdamodheatmapnum, reducedDimName = "celda_UMAP")})
    showNotification("Module heatmap complete.")
  })

  observe({
    if(!is.null(vals$counts)){
      #If data is uploaded in data tab, enable first tab i.e. Normalization tab in Seurat workflow
      shinyjs::enable(
        selector = "div[value='Identify Number of Feature Modules']")
    }else{
      #If no data uploaded in data tab, disabled all tabs and plots.

      #Disable tabs
      shinyjs::disable(
        selector = "div[value='Identify Number of Feature Modules']")
      shinyjs::disable(
        selector = "div[value='Identify Number of Cell Clusters']")
      shinyjs::disable(
        selector = "div[value='Visualization']")

      #Disable plots inside Modsplit subtab
      shinyjs::disable(
        selector = ".celda_modsplit_plots a[data-value='Perplexity Plot']")
      shinyjs::disable(
        selector = ".celda_modsplit_plots a[data-value='Perplexity Diff Plot']")
    }
  })

  #-----------------------------------------------------------------------------
  # Page 3.3: Cell Viewer ####
  #-----------------------------------------------------------------------------
  #-+-+-+-+-+-For Functional Panel collapse##############
  shinyjs::onclick("cv_button1", shinyjs::toggle(id = "cv_collapse1",
                                                 anim = TRUE), add = TRUE)
  shinyjs::onclick("cv_button2", shinyjs::toggle(id = "cv_collapse2",
                                                 anim = TRUE), add = TRUE)
  shinyjs::onclick("cv_button3", shinyjs::toggle(id = "cv_collapse3",
                                                 anim = TRUE), add = TRUE)
  shinyjs::addClass(id = "cv_button1", class = "btn-block")
  shinyjs::addClass(id = "cv_button2", class = "btn-block")
  shinyjs::addClass(id = "cv_button3", class = "btn-block")
  colorbrewer_list <- rownames(RColorBrewer::brewer.pal.info)
  color_table <- RColorBrewer::brewer.pal.info %>% data.frame()
  color_seqdiv <- rownames(color_table[which(color_table$category == "div"
                                             |color_table$category == "seq"),])

  #-+-+-+-+-+-For Input Observe##############
  observeEvent(input$navbar,{
    if (input$navbar == "CellViewer"){
      # is there an error or not
      if (is.null(vals$counts)){
        shinyalert::shinyalert("Error!", "Upload data first.", type = "error")
      }else{
        gene_list <- rownames(vals$counts)
        annotation_list <- names(colData(vals$counts))
        annotation_list2 <- list()
        for (i in 1:length(annotation_list)){
          if(!all.is.numeric(vals$counts[[annotation_list[i]]])){
            annotation_list2$Categorical <- c(annotation_list2$Categorical, annotation_list[i])
          }else{
            annotation_list2$Numeric <- c(annotation_list2$Numeric, annotation_list[i])
          }
        }
        annotation_list <- annotation_list2
        rm(annotation_list2)
        updateSelectizeInput(session, "GeneSelect_Assays_Xaxis",
                          choices = c(gene_list), server = TRUE)
        updateSelectInput(session, "AnnotationSelect_Xaxis",
                          choices = c(annotation_list))
        updateSelectizeInput(session, "GeneSelect_Assays_Yaxis",
                          choices = c(gene_list), server = TRUE)
        updateSelectInput(session, "AnnotationSelect_Yaxis",
                          choices = c(annotation_list))
        updateSelectizeInput(session, "GeneSelect_Assays_Colorby",
                          choices = c(gene_list), server = TRUE)
        updateSelectInput(session, "AnnotationSelect_Colorby",
                          choices = c(annotation_list))
        updateSelectizeInput(session, "adjustgroupby", label = NULL, choices = c("None", annotation_list))
        updateSelectizeInput(session,"adjustbrewer", label = "Color Palettes:",
                             choices = c("RdYlBu",color_seqdiv))
      }
    }

    # if(input$navbar == "Feature Selection & Dimensionality Reduction"){
    #   gene_list <- rownames(vals$counts)
    #   updateSelectizeInput(session, "scatterFSGenes",
    #                        choices = c(gene_list),
    #                        server = TRUE)
    # }
  })

  hide_TypeSelect <- reactiveVal("hide")
  hide_bins <- reactiveVal()

  observeEvent(input$viewertabs, {
    if(!is.null(vals$counts)) {
      if(!is.null(reducedDims(vals$counts))) {
        approach_list <- names(reducedDims(vals$counts))
        if (input$viewertabs != "Scatter Plot") {
          updateSelectInput(session, "QuickAccess",
                            choices = c("Custom"))
          shinyjs::delay(5,shinyjs::disable("QuickAccess"))

          updateSelectInput(session, "TypeSelect_Xaxis",
                            choices = c("None", "Cell Annotation"))
          updateSelectInput(session, "TypeSelect_Yaxis",
                            choices = c("Expression Assays", "Cell Annotation"))
          updateSelectInput(session, "TypeSelect_Colorby",
                            selected = "Single Color")
          updateSelectInput(session, "adjustgroupby",
                            selected = "None")
          updatePrettyToggle(session, "checkColorbinning",
                             value = FALSE)
          hide_TypeSelect("hide")
          shinyjs::delay(5,shinyjs::disable("TypeSelect_Colorby"))
          shinyjs::delay(5,shinyjs::disable("adjustgroupby"))

          shinyjs::delay(5, shinyjs::disable("adjustlegendtitle"))
          shinyjs::delay(5, shinyjs::disable("adjustlegendtitlesize"))
          shinyjs::delay(5, shinyjs::disable("adjustlegendsize"))
        } else {
          updateSelectInput(session, "QuickAccess",
                            choices = c("", approach_list, "Custom"))
          shinyjs::delay(5,shinyjs::enable("QuickAccess"))

          updateSelectInput(session, "TypeSelect_Xaxis",
                            choices = c("Reduced Dimensions", "Expression Assays", "Cell Annotation"))
          updateSelectInput(session, "TypeSelect_Yaxis",
                            choices = c("Reduced Dimensions", "Expression Assays", "Cell Annotation"))
          updateSelectInput(session, "TypeSelect_Colorby",
                            selected = "Single Color")
          updateSelectInput(session, "adjustgroupby",
                            selected = "None")
          updatePrettyToggle(session, "checkColorbinning",
                             value = FALSE)
          hide_TypeSelect("hide")
          shinyjs::delay(5,shinyjs::enable("TypeSelect_Colorby"))
          shinyjs::delay(5,shinyjs::enable("adjustgroupby"))
          shinyjs::delay(5, shinyjs::enable("adjustlegendtitle"))
          if (!is.null(input$adjustgridlines) &
              isFALSE(input$adjustgridlines)) {
            shinyjs::delay(5, shinyjs::enable("adjustlegendtitlesize"))
            shinyjs::delay(5, shinyjs::enable("adjustlegendsize"))
          }
        }

        if (input$viewertabs != "Bar Plot") {
          shinyjs::delay(5, shinyjs::enable("adjustalpha"))
          shinyjs::delay(5, shinyjs::enable("adjustsize"))
        } else {
          shinyjs::delay(5, shinyjs::disable("adjustalpha"))
          shinyjs::delay(5, shinyjs::disable("adjustsize"))
        }
      }
    }
  })

  observeEvent(input$adjustgridlines, {
    req(vals$counts)
    if (!is.null(input$adjustgridlines)) {
      if (isTRUE(input$adjustgridlines)) {
        shinyjs::delay(5, shinyjs::disable("adjustlegendtitlesize"))
        shinyjs::delay(5, shinyjs::disable("adjustlegendsize"))
        shinyjs::delay(5, shinyjs::disable("adjustaxissize"))
        shinyjs::delay(5, shinyjs::disable("adjustaxislabelsize"))
      } else {
        if (input$viewertabs == "Scatter Plot") {
          shinyjs::delay(5, shinyjs::enable("adjustlegendtitlesize"))
          shinyjs::delay(5, shinyjs::enable("adjustlegendsize"))
        }
        shinyjs::delay(5, shinyjs::enable("adjustaxissize"))
        shinyjs::delay(5, shinyjs::enable("adjustaxislabelsize"))
      }
    }
  })


  #-+-+-+-+-+-For Advanced Input Observe##############
  ###ApproachSelect to DimensionSelect X-Axis
  observeEvent(input$ApproachSelect_Xaxis, {
    if (!is.null(vals$counts)){
      len <- length(SingleCellExperiment::reducedDims(vals$counts))
      if (!is.null(input$ApproachSelect_Xaxis) & len > 0){
        Df <- data.frame(SingleCellExperiment::reducedDim(vals$counts,input$ApproachSelect_Xaxis))
        xs <- colnames(Df)
        updateSelectInput(session, "ColumnSelect_Xaxis", choices = c(xs))
        rm(Df)
      }
    }
  })
  ###ApproachSelect to DimensionSelect Y-Axis
  observeEvent(input$ApproachSelect_Yaxis, {
    if (!is.null(vals$counts)){
      len <- length(SingleCellExperiment::reducedDims(vals$counts))
      if (!is.null(input$ApproachSelect_Yaxis) & len > 0){
        Df2 <- data.frame(SingleCellExperiment::reducedDim(vals$counts,input$ApproachSelect_Yaxis))
        xs2 <- colnames(Df2)
        xs2 <- sort(xs2, decreasing = TRUE)
        updateSelectInput(session, "ColumnSelect_Yaxis", choices = c(xs2))
        rm(Df2)
      }
    }
  })
  ###ApproachSelect to DimensionSelect Colorby
  observeEvent(input$ApproachSelect_Colorby, {
    if (!is.null(vals$counts)){
      len <- length(SingleCellExperiment::reducedDims(vals$counts))
      if (!is.null(input$ApproachSelect_Colorby) & len > 0){
        Df3 <- data.frame(SingleCellExperiment::reducedDim(vals$counts,input$ApproachSelect_Colorby))
        xs3 <- colnames(Df3)
        prefix <- input$ApproachSelect_Colorby
        suffix <- seq(1:length(xs3))
        columns <- paste(prefix, suffix, sep = "_")
        updateSelectInput(session, "ColumnSelect_Colorby", choices = c(columns))
        rm(Df3)
      }
    }
  })

  #-+-+-+-+-+-Observe Color by###################################################
  ###Observe Radio Button Select Value Type
  # input$AnnotationSelect_Colorby,

  observe({
    # All inputs to listen for
    input$TypeSelect_Colorby
    input$AnnotationSelect_Colorby

    #Reduced Dimensions
    input$ApproachSelect_Colorby
    input$ColumnSelect_Colorby

    #Expression Assay
    input$AdvancedMethodSelect_Colorby
    input$GeneSelect_Assays_Colorby

    if(input$TypeSelect_Colorby == 'Cell Annotation'){
      ###If Cell Annotation is not numeric
      if(!is.numeric(colData(vals$counts)@listData[[input$AnnotationSelect_Colorby]])){
        updateRadioButtons(session, "SelectColorType", "Categorical or Continuous",
                           choices = c("Categorical", "Continuous"),
                           selected = "Categorical")
        hide_TypeSelect("hide")
      }else if(is.integer(colData(vals$counts)@listData[[input$AnnotationSelect_Colorby]])
               &length(levels(as.factor(colData(vals$counts)@listData[[input$AnnotationSelect_Colorby]])))<=25){
        updateRadioButtons(session, "SelectColorType", "Categorical or Continuous",
                           choices = c("Categorical", "Continuous"),
                           selected = "Categorical")
        hide_TypeSelect("show")
      }else{
        updateRadioButtons(session, "SelectColorType", "Categorical or Continuous",
                           choices = c("Categorical", "Continuous"),
                           selected = "Continuous")
        hide_TypeSelect("hide")
      }
    } else if(input$TypeSelect_Colorby == 'Reduced Dimensions'){
      updateRadioButtons(session, "SelectColorType", "Categorical or Continuous",
                         choices = c("Categorical", "Continuous"),
                         selected = "Continuous")
      hide_TypeSelect("hide")
    } else if(input$TypeSelect_Colorby == "Expression Assays"){
      updateRadioButtons(session, "SelectColorType", "Categorical or Continuous",
                         choices = c("Categorical", "Continuous"),
                         selected = "Continuous")
      hide_TypeSelect("hide")
    } else {
      # single color
      updateRadioButtons(session, "SelectColorType", "Categorical or Continuous",
                         choices = c("Categorical", "Continuous"),
                         selected = "Categorical")
      hide_TypeSelect("hide")
    }
  })

  observeEvent(input$SelectColorType,{
    if(input$SelectColorType == "Categorical"){
      hide_bins("hide")
    }else{
      hide_bins("show")
    }
  })

  output$hide_typebtns <- renderText({
    hide_TypeSelect()
  })

  outputOptions(output, "hide_typebtns", suspendWhenHidden = FALSE)

  output$hide_bins <- renderText({
    hide_bins()
  })

  outputOptions(output, "hide_bins", suspendWhenHidden = FALSE)

  numColors <- NULL
  colorLabels <- NULL
  ### Observe input changes that should trigger categorical color generator
  observeEvent(c(input$SelectColorType, input$TypeSelect_Colorby, input$AnnotationSelect_Colorby, input$colorTheme), {
    if (input$TypeSelect_Colorby == "Single Color") {
      shinyjs::hide("categoricalColorConditional")
      shinyjs::hide("continuousColorConditional")
    }
    if (input$SelectColorType == "Categorical" && input$TypeSelect_Colorby != "Single Color") {
      if(input$TypeSelect_Colorby == "Cell Annotation") {
        req(vals$counts, input$AnnotationSelect_Colorby)
        labels = sort(unique(SingleCellExperiment::colData(vals$counts)[, input$AnnotationSelect_Colorby]))
        if (length(labels) <=25 ) {
          colorLabels <<- labels
          numColors <<- length(labels)
          defaultColors <- discreteColorPalette(numColors, input$colorTheme)
          output$categoricalColorUI <- renderUI({
            lapply(1:numColors, function(i){
              colourInput(inputId=paste0(i, "_color"), label=labels[i], value=defaultColors[i], showColour="background")
            })
          })
          shinyjs::show("categoricalColorConditional")
          shinyjs::hide("continuousColorConditional")
        }
      }
    } else if (input$SelectColorType == "Continuous" && input$TypeSelect_Colorby != "Single Color") {
      shinyjs::hide("categoricalColorConditional")
      shinyjs::show("continuousColorConditional")
    }
  })

  #-+-+-+-+-+-cellviewer prepare step1: choose data. (next steps included)###########################################################
  #cellviewer <- eventReactive(input$runCellViewer,{
  observeEvent(input$runCellViewer,{
    colors <- c()
    if (!is.null(numColors) && input$SelectColorType == 'Categorical') {
      for (i in 1: numColors) {
        colors[i] <- input[[ paste0(i,"_color")]]
      }
      names(colors) = colorLabels
    }
    #-+-+-+-+-+-cellviewer prepare3 : prepare Axis Label Name#####################
    ###Xaxis label name
    if (!is.null(input$adjustxlab) &
        input$adjustxlab != "") {
      xname <- input$adjustxlab
    } else {
      if (input$QuickAccess != "Custom" &
          input$QuickAccess != "") {
        # reddim selected
        xname <- paste0(input$QuickAccess, 1)
      } else if (input$TypeSelect_Xaxis == 'Reduced Dimensions') {
        xname <- paste0(input$ApproachSelect_Xaxis, "_",
                        substr(input$ColumnSelect_Xaxis,
                               str_length(input$ColumnSelect_Xaxis),
                               str_length(input$ColumnSelect_Xaxis)))
      } else if (input$TypeSelect_Xaxis == 'Expression Assays') {
        xname <- input$GeneSelect_Assays_Xaxis
      } else if (input$TypeSelect_Xaxis == "Cell Annotation") {
        xname <- input$AnnotationSelect_Xaxis
      } else {
        xname <- ""
      }
    }
    xname <- gsub("-", "_", xname)
    ###Yaxis label name
    if (!is.null(input$adjustylab) &
        input$adjustylab != "") {
      yname <- input$adjustylab
    } else {
      if (input$QuickAccess != "Custom" &
          input$QuickAccess != "") {
        # reddim selected
        yname <- paste0(input$QuickAccess, 2)
      } else if (input$TypeSelect_Yaxis == 'Reduced Dimensions') {
        yname <- paste0(input$ApproachSelect_Yaxis, "_",
                        substr(input$ColumnSelect_Yaxis,
                               str_length(input$ColumnSelect_Yaxis),
                               str_length(input$ColumnSelect_Yaxis)))
      } else if (input$TypeSelect_Yaxis == 'Expression Assays') {
        yname <- input$GeneSelect_Assays_Yaxis
      } else {
        yname <- input$AnnotationSelect_Yaxis
      }
    }
    yname <- gsub("-", "_", yname)
    ###Legend name
    if (input$TypeSelect_Colorby != 'Pick a Color') {
      if (input$TypeSelect_Colorby == 'Reduced Dimensions' && input$adjustlegendtitle == "") {
        legendname <- paste0(input$ApproachSelect_Colorby,"_",substr(input$ColumnSelect_Colorby,
                                                                     str_length(input$ColumnSelect_Colorby),str_length(input$ColumnSelect_Colorby)))
      } else if (input$TypeSelect_Colorby == 'Expression Assays' && input$adjustlegendtitle == "") {
        legendname <- input$GeneSelect_Assays_Colorby
      } else if (input$adjustlegendtitle == "") {
        legendname <- input$AnnotationSelect_Colorby
      } else {
        legendname <- input$adjustlegendtitle
      }
    }
    legendname <- gsub("-", "_", legendname)
    #-+-+-+-+-+-cellviewer prepare4 : choose group by and create plotly function###################
    pltVars <- list()
    if(input$viewertabs == "Violin/Box Plot" || input$viewertabs == "Bar Plot"){
      if(input$TypeSelect_Xaxis == "None"){
        pltVars$groupby <- NULL
      }else if(input$TypeSelect_Xaxis == "Expression Assays"){
        pltVars$groupby <- input$GeneSelect_Assays_Xaxis
      }else if(input$TypeSelect_Xaxis == "Cell Annotation"){
        pltVars$groupby <- input$AnnotationSelect_Xaxis
      }
    }else if(input$adjustgroupby != "None"){
      pltVars$groupby <- input$adjustgroupby
    }else{
      pltVars$groupby <- NULL
    }
    if (input$checkColorbinning == TRUE && input$SelectColorType == "Continuous"){
      pltVars$bin <- input$adjustColorbinning
    }else{
      pltVars$bin <- NULL
    }
    if (input$SelectColorType == "Categorical"){
      pltVars$class <- "factor"
    }else{
      pltVars$class <- "numeric"
    }

    if(input$adjustgridlines == TRUE){
      pltVars$defTheme <- FALSE
    }else{
      pltVars$defTheme <- TRUE
    }

    if(input$viewertabs == "Scatter Plot"){
      #### Prepare Custom plotting matrix axis ####
      if (input$QuickAccess == "Custom") {
        # X Axis
        message("CellViewer: Custom plotting mode, making up the axis")
        if (input$TypeSelect_Xaxis == "Expression Assays") {
          message("X axis: Using expression of ", input$GeneSelect_Assays_Xaxis,
                  " from ", input$AdvancedMethodSelect_Xaxis)
          xvec <- expData(vals$counts, input$AdvancedMethodSelect_Xaxis)[input$GeneSelect_Assays_Xaxis,]
        } else if (input$TypeSelect_Xaxis == "Reduced Dimensions") {
          message("X axis: Using dimension reduction ", input$ColumnSelect_Xaxis,
                  " from ", input$ApproachSelect_Xaxis)
          xvec <- reducedDim(vals$counts, input$ApproachSelect_Xaxis)[,input$ColumnSelect_Xaxis]
        } else if (input$TypeSelect_Xaxis == 'Cell Annotation') {
          message("X axis: Using cell annotation ",
                  input$AnnotationSelect_Xaxis)
          xvec <- vals$counts[[input$AnnotationSelect_Xaxis]]
        }
        # Y Axis
        if (input$TypeSelect_Yaxis == "Expression Assays") {
          message("Y axis: Using expression of ", input$GeneSelect_Assays_Yaxis,
                  " from ", input$AdvancedMethodSelect_Yaxis)
          yvec <- expData(vals$counts, input$AdvancedMethodSelect_Yaxis)[input$GeneSelect_Assays_Yaxis,]
        } else if (input$TypeSelect_Yaxis == "Reduced Dimensions") {
          message("Y axis: Using dimension reduction ", input$ColumnSelect_Yaxis,
                  " from ", input$ApproachSelect_Yaxis)
          yvec <- reducedDim(vals$counts, input$ApproachSelect_Yaxis)[,input$ColumnSelect_Yaxis]
        } else if (input$TypeSelect_Yaxis == 'Cell Annotation') {
          message("Y axis: Using cell annotation ",
                  input$AnnotationSelect_Yaxis)
          yvec <- vals$counts[[input$AnnotationSelect_Yaxis]]
        }
        # Merge and insert to reducedDim(sce, "Custom")
        customMat <- matrix(c(xvec, yvec), nrow = length(xvec))
        colnames(customMat) <- c(xname, yname)
        rownames(customMat) <- names(xvec)
        reducedDim(vals$counts, "Custom") <- customMat
      }
      if(input$TypeSelect_Colorby == "Single Color"){
        a <- plotSCEScatter(vals$counts, reducedDimName = input$QuickAccess,
                            xlab = xname, ylab = yname, title = input$adjusttitle, groupBy = pltVars$groupby,
                            transparency = input$adjustalpha, dotSize = input$adjustsize, combinePlot = "none",
                            axisSize = input$adjustaxissize, axisLabelSize = input$adjustaxislabelsize,
                            legendSize = input$adjustlegendsize, legendTitleSize = input$adjustlegendtitlesize,
                            conditionClass = pltVars$class, defaultTheme = as.logical(pltVars$defTheme))
      }else if(input$TypeSelect_Colorby == "Expression Assays"){
        a <- plotSCEDimReduceFeatures(vals$counts, feature = input$GeneSelect_Assays_Colorby,
                                      reducedDimName = input$QuickAccess, useAssay = input$AdvancedMethodSelect_Colorby,
                                      xlab = xname, ylab = yname, legendTitle = legendname, title = input$adjusttitle,
                                      groupBy = pltVars$groupby, bin = pltVars$bin, transparency = input$adjustalpha,
                                      colorLow = input$lowColor, colorMid = input$midColor, colorHigh = input$highColor,
                                      dotSize = input$adjustsize, combinePlot = "none", axisSize = input$adjustaxissize,
                                      axisLabelSize = input$adjustaxislabelsize, legendSize = input$adjustlegendsize,
                                      legendTitleSize = input$adjustlegendtitlesize)
      } else if (input$TypeSelect_Colorby == "Cell Annotation") {
        a <- plotSCEDimReduceColData(vals$counts,reducedDimName = input$QuickAccess,xlab = xname,ylab = yname,
                                     colorBy = input$AnnotationSelect_Colorby,groupBy = pltVars$groupby,legendTitle = legendname,
                                     title = input$adjusttitle,bin = pltVars$bin,transparency = input$adjustalpha,colorScale = colors,
                                     colorLow = input$lowColor, colorMid = input$midColor, colorHigh = input$highColor, dotSize = input$adjustsize,
                                     combinePlot = "none",axisSize = input$adjustaxissize,axisLabelSize = input$adjustaxislabelsize,
                                     legendSize = input$adjustlegendsize,legendTitleSize = input$adjustlegendtitlesize,conditionClass = pltVars$class)
      }else if(input$TypeSelect_Colorby == "Reduced Dimensions"){
        a <- plotSCEScatter(vals$counts, reducedDimName = input$QuickAccess, slot = "reducedDims",
                            annotation = input$ColumnSelect_Colorby, transparency = input$adjustalpha,
                            colorLow = input$lowColor, colorMid = input$midColor, colorHigh = input$highColor,
                            groupBy = pltVars$groupby, title = input$adjusttitle, legendTitle = legendname,
                            xlab = xname, ylab = yname, dotSize = input$adjustsize, bin = pltVars$bin,
                            combinePlot = "none", axisSize = input$adjustaxissize, axisLabelSize = input$adjustaxislabelsize,
                            legendSize = input$adjustlegendsize, legendTitleSize = input$adjustlegendtitlesize)
      }
    }else if(input$viewertabs == "Bar Plot"){
      if(input$TypeSelect_Yaxis == "Expression Assays"){
        a <- plotSCEBarAssayData(vals$counts, title = input$adjusttitle, xlab = xname, ylab = yname,
                                 useAssay = input$AdvancedMethodSelect_Yaxis, groupBy = pltVars$groupby,
                                 feature = input$GeneSelect_Assays_Yaxis,
                                 combinePlot = "none", axisSize = input$adjustaxissize,
                                 axisLabelSize = input$adjustaxislabelsize, defaultTheme = as.logical(pltVars$defTheme))
      }else if(input$TypeSelect_Yaxis == "Cell Annotation"){
        a <- plotSCEBarColData(vals$counts, title = input$adjusttitle, xlab = xname, ylab = yname,
                               coldata = input$AnnotationSelect_Yaxis, groupBy = pltVars$groupby,
                               combinePlot = "none",
                               axisSize = input$adjustaxissize, axisLabelSize = input$adjustaxislabelsize,
                               defaultTheme = as.logical(pltVars$defTheme))
      }
    }else if(input$viewertabs == "Violin/Box Plot"){
      if(isTRUE(input$vlnboxcheck)){
        vln <- TRUE
        bx <- FALSE
      }else if(isFALSE(input$vlnboxcheck)){
        vln <- FALSE
        bx <- TRUE
      }
      if(input$TypeSelect_Yaxis == "Expression Assays"){
        a <- plotSCEViolinAssayData(vals$counts, violin = vln, box = bx, xlab = xname, ylab = yname,
                                    useAssay = input$AdvancedMethodSelect_Yaxis, title = input$adjusttitle,
                                    feature = input$GeneSelect_Assays_Yaxis, groupBy = pltVars$groupby,
                                    transparency = input$adjustalpha, dotSize = input$adjustsize, combinePlot = "none",
                                    axisSize = input$adjustaxissize, axisLabelSize = input$adjustaxislabelsize,
                                    defaultTheme = as.logical(pltVars$defTheme))
      }else if(input$TypeSelect_Yaxis == "Cell Annotation"){
        a <- plotSCEViolinColData(vals$counts, title = input$adjusttitle, xlab = xname, ylab = yname,
                                  coldata = input$AnnotationSelect_Yaxis, violin = vln, box = bx,
                                  groupBy = pltVars$groupby, transparency = input$adjustalpha,
                                  dotSize = input$adjustsize, combinePlot = "none", axisSize = input$adjustaxissize,
                                  axisLabelSize = input$adjustaxislabelsize, defaultTheme = as.logical(pltVars$defTheme))
      }
    }
    if (input$TypeSelect_Colorby == "Single Color"){
      a$layers[[1]]$aes_params$colour <- input$Col
    }
    if (isTRUE(input$adjustgridlines)){
      a <- a + ggplot2::theme_bw()
    }
    a <- plotly::ggplotly(a)
    output$scatter <- renderPlotly({
      plotly::subplot(plotlist = a, titleX = TRUE, titleY = TRUE)
    })
  })

  #-----------------------------------------------------------------------------
  # Page 3.4: Heatmap ####
  #-----------------------------------------------------------------------------

  hmTemp <- reactiveValues(
    sce = NULL,
    cellIndex = NULL,
    geneIndex = NULL,
    colDataName = NULL,
    rowDataName = NULL,
    colSplitBy = NULL,
    rowSplitBy = NULL,
    cellTableCol = NULL,
    geneTableCol = NULL,
    colColorPresets = list(),
    rowColorPresets = list()
  )

  observeEvent(vals$counts, {
    if(!is.null(vals$counts)){
      hmTemp$sce <- vals$counts
    }
  })

  # Heatmap: Import Analysis ####
  output$hmImpDEGUI <- renderUI({
    if(!is.null(vals$counts)){
      if("diffExp" %in% names(metadata(vals$counts))){
        analysis <- names(metadata(vals$counts)$diffExp)
        selectInput('hmImpDEG', "Import results from analysis:",
                    analysis)
      } else {
        p("Differential expression analysis not performed yet.")
      }
    }
  })

  observeEvent(input$hmImportRun, {
    if(!is.null(vals$counts)){
      if(!is.null(input$hmImport) &&
         input$hmImport == "Differential Expression"){
        if(!is.null(input$hmImpDEG)){
          result <- metadata(vals$counts)$diffExp[[input$hmImpDEG]]
          useAssay <- result$useAssay
          updateSelectInput(session, "hmAssay", selected = useAssay)
          method <- result$method
          # Cell side
          addColData <- data.frame(row.names = colnames(vals$counts))
          idx <- rep(NA, ncol(vals$counts))
          idx[result$select$ix1] <- result$groupNames[1]
          idx[result$select$ix2] <- result$groupNames[2]
          hmTemp$cellIndex <- which(!is.na(idx))
          conditionColName <- paste(method, input$hmImpDEG, "condition",
                                    sep = '_')
          addColData[[conditionColName]] <- factor(idx,
                                                   levels = result$groupNames)
          colData(hmTemp$sce) <- cbind(colData(hmTemp$sce), addColData)
          hmTemp$cellTableCol <- conditionColName
          hmTemp$colDataName <- conditionColName
          hmTemp$colSplitBy <- conditionColName
          hmTemp$colColorPresets[[conditionColName]] <- c('red', 'cyan',
                                                          'white')
          names(hmTemp$colColorPresets[[conditionColName]]) <-
            c(result$groupNames, "NA")
          # Gene side
          addRowData <- data.frame(row.names = rownames(vals$counts))
          deg <- result$result
          deg <- deg[stats::complete.cases(deg),]
          logFCColName <- paste(method, input$hmImpDEG, "Log2FC",
                                sep = '_')
          FDRColName <- paste(method, input$hmImpDEG, "FDR",
                              sep = '_')
          addRowData[deg$Gene, logFCColName] <- deg$Log2_FC
          addRowData[deg$Gene, FDRColName] <- deg$FDR
          regColName <- paste(method, input$hmImpDEG, "regulation",
                              sep = '_')
          degUp <- deg[deg$Log2_FC > 0,]
          degDown <- deg[deg$Log2_FC < 0,]
          addRowData[degUp$Gene, regColName] <- "up"
          addRowData[degDown$Gene, regColName] <- "down"
          addRowData[[regColName]] <- factor(addRowData[[regColName]],
                                             levels = c('up', 'down'))
          rowData(hmTemp$sce) <- cbind(rowData(hmTemp$sce), addRowData)
          hmTemp$geneTableCol <- c(regColName, logFCColName, FDRColName)
          hmTemp$geneIndex <- which(rownames(vals$counts) %in% deg$Gene)
          hmTemp$rowDataName <- regColName
          hmTemp$rowSplitBy <- regColName
          hmTemp$rowColorPresets[[regColName]] <- c('red', 'cyan', 'white')
          names(hmTemp$rowColorPresets[[regColName]]) <- c('up', 'down', 'NA')
        }
      } else if (!is.null(input$hmImport) &&
                 input$hmImport == "Find Marker"){
        markerTable <- metadata(vals$counts)$findMarker
        if(!is.null(markerTable) &&
           dim(markerTable)[1] > 0){
          markerTable <- markerTable[stats::complete.cases(markerTable),]
          # Cell side
          cluster <- colnames(markerTable)[5]
          hmTemp$cellIndex <- seq_len(ncol(hmTemp$sce))
          hmTemp$colDataName <- cluster
          hmTemp$cellTableCol <- cluster
          hmTemp$colSplitBy <- cluster
          # Gene side
          dup.gene <- unique(markerTable$Gene[duplicated(markerTable$Gene)])
          for(g in dup.gene){
            deg.gix <- markerTable$Gene == g
            deg.gtable <- markerTable[deg.gix,]
            toKeep <- which.max(deg.gtable$Log2_FC)
            toRemove <- which(deg.gix)[-toKeep]
            markerTable <- markerTable[-toRemove,]
          }
          hmTemp$geneIndex <- which(rownames(vals$counts) %in% markerTable$Gene)
          addRowData <- data.frame(row.names = rownames(vals$counts))
          addRowData[markerTable$Gene, "Marker_for_Cluster"] <- markerTable[,5]
          addRowData[markerTable$Gene, "findMarker_Log2FC"] <-
            markerTable$Log2_FC
          addRowData[markerTable$Gene, "findMarker_FDR"] <- markerTable$FDR
          rowData(hmTemp$sce) <- cbind(rowData(hmTemp$sce), addRowData)
          hmTemp$geneTableCol <- c("Marker_for_Cluster",
                                   "findMarker_Log2FC",
                                   "findMarker_FDR")
          hmTemp$rowDataName <- "Marker_for_Cluster"
          hmTemp$rowSplitBy <- "Marker_for_Cluster"
          hmTemp$rowColorPresets$Marker_for_Cluster <-
            hmAnnAllColors$col[[cluster]]
        }
      }
    }
  })
  # Heatmap: Subsetting Cells####
  output$hmCellColUI <- renderUI({
    if(!is.null(vals$counts)){
      selectInput(
        'hmCellCol',
        "Columns to display",
        names(colData(hmTemp$sce)), multiple = TRUE, width = '550px',
        selected = hmTemp$cellTableCol)
    }
  })

  output$hmCellColTable <- DT::renderDataTable({
    if(!is.null(vals$counts)){
      df <- as.data.frame(colData(hmTemp$sce))
      rowNameCol <- data.frame(Row_Names = colnames(vals$counts))
      df <- cbind(rowNameCol, df)
      rownames(df) <- NULL
      DT::datatable(
        df,
        filter = 'top', options = list(stateSave = TRUE, scrollX = TRUE)
      )
    }
  }, server = TRUE)

  hmCellColTable_proxy <- DT::dataTableProxy("hmCellColTable")

  observeEvent(input$hmCellCol, {
    colNames <- c('Row_Names', names(colData(hmTemp$sce)))
    showIdx <- which(colNames %in% input$hmCellCol)
    showIdx <- c(1, showIdx)
    DT::showCols(hmCellColTable_proxy, showIdx, reset = TRUE)
  })

  observeEvent(input$hmCellColTable_state, {
    DT::selectRows(hmCellColTable_proxy, hmTemp$cellIndex)
  })

  observeEvent(input$hmCellColTable_rows_selected, {
    hmTemp$cellIndex <- input$hmCellColTable_rows_selected
  })

  observeEvent(input$hmCellColTable_addAll, {
    DT::selectRows(hmCellColTable_proxy,
                   sort(unique(c(input$hmCellColTable_rows_selected,
                                 input$hmCellColTable_rows_all))))
    hmTemp$cellIndex <- sort(unique(c(input$hmCellColTable_rows_selected,
                                      input$hmCellColTable_rows_all)))
  })

  observeEvent(input$hmCellColTable_clear, {
    DT::selectRows(hmCellColTable_proxy, NULL)
    hmTemp$cellIndex <- NULL
  })

  observeEvent(input$hmCellCol, {
    hmTemp$cellTableCol <- input$hmCellCol
  })

  output$hmCellNEnteredUI <- renderUI({
    inputList <- str_trim(scan(text = input$hmCellText,
                               sep='\n', what = 'character', quiet = TRUE))
    uniqInput <- unique(inputList)
    nInput <- length(uniqInput)
    if(!is.null(vals$counts) && nInput > 0){
      if(!is.null(input$hmCellTextBy) && input$hmCellTextBy == 'Row Names'){
        BY <- NULL
      } else {
        BY <- input$hmCellTextBy
      }
      matched <- retrieveSCEIndex(vals$counts, uniqInput, 'cell',
                                  by = BY, exactMatch = input$hmCellTextEM,
                                  firstMatch = input$hmCellTextFM)
      nMatched <- length(matched)
    } else {
      nMatched <- 0
    }

    p(paste0(nInput, " unique input, ", nMatched, "matched."))
  })

  observeEvent(input$hmCellAddFromText, {
    if(!is.null(vals$counts)){
      inputList <- str_trim(scan(text = input$hmCellText,
                                 sep='\n', what = 'character', quiet = TRUE))
      uniqInput <- unique(inputList)
      if(length(uniqInput) > 0){
        if(!is.null(input$hmCellTextBy) && input$hmCellTextBy == 'Row Names'){
          BY <- NULL
        } else {
          BY <- input$hmCellTextBy
        }
        newIdx <- retrieveSCEIndex(vals$counts, uniqInput, 'cell',
                                   by = BY, exactMatch = input$hmCellTextEM,
                                   firstMatch = input$hmCellTextFM)
        DT::selectRows(hmCellColTable_proxy,
                       sort(unique(c(input$hmCellColTable_rows_selected,
                                     newIdx))))
      }
    }
  })

  output$hmCellSumUI <- renderUI({
    nCell <- length(hmTemp$cellIndex)
    if(nCell == 0){
      p("No cells selected, going to use them all", style = 'margin-top: 5px;')
    } else {
      p(paste0("Totally ", nCell, " cells selected."),
        style = 'margin-top: 5px;')
    }
  })

  # Heatmap: Subsetting Genes ####
  output$hmGeneColUI <- renderUI({
    if(!is.null(vals$counts)){
      selectInput(
        'hmGeneCol',
        "Columns to display",
        names(rowData(hmTemp$sce)), multiple = TRUE, width = '550px',
        selected = hmTemp$geneTableCol)
    }
  })

  output$hmGeneColTable <- DT::renderDataTable({
    if(!is.null(vals$counts)){
      df <- as.data.frame(rowData(hmTemp$sce))
      rowNameCol <- data.frame(Row_Names = rownames(vals$counts))
      df <- cbind(rowNameCol, df)
      rownames(df) <- NULL
      DT::datatable(
        df,
        filter = 'top', options = list(stateSave = TRUE, scrollX = TRUE)
      )
    }
  }, server = TRUE)

  hmGeneColTable_proxy <- DT::dataTableProxy("hmGeneColTable")

  observeEvent(input$hmGeneCol, {
    colNames <- c('Row_Names', names(rowData(hmTemp$sce)))
    showIdx <- which(colNames %in% input$hmGeneCol)
    showIdx <- c(1, showIdx)
    DT::showCols(hmGeneColTable_proxy, showIdx, reset = TRUE)
  })

  observeEvent(input$hmGeneColTable_state, {
    DT::selectRows(hmGeneColTable_proxy, hmTemp$geneIndex)
  })

  observeEvent(input$hmGeneColTable_rows_selected, {
    hmTemp$geneIndex <- input$hmGeneColTable_rows_selected
  })

  observeEvent(input$hmGeneColTable_addAll, {
    DT::selectRows(hmGeneColTable_proxy,
                   sort(unique(c(input$hmGeneColTable_rows_selected,
                                 input$hmGeneColTable_rows_all))))
    hmTemp$geneIndex <- sort(unique(c(input$hmGeneColTable_rows_selected,
                                      input$hmGeneColTable_rows_all)))
  })

  observeEvent(input$hmGeneColTable_clear, {
    DT::selectRows(hmGeneColTable_proxy, NULL)
    hmTemp$geneIndex <- NULL
  })

  observeEvent(input$hmGeneCol, {
    hmTemp$geneTableCol <- input$hmGeneCol
  })

  output$hmGeneNEnteredUI <- renderUI({
    inputList <- str_trim(scan(text = input$hmGeneText,
                               sep='\n', what = 'character', quiet = TRUE))
    uniqInput <- unique(inputList)
    nInput <- length(uniqInput)
    if(!is.null(vals$counts) && nInput > 0){
      if(!is.null(input$hmGeneTextBy) && input$hmGeneTextBy == 'Row Names'){
        BY <- NULL
      } else {
        BY <- input$hmGeneTextBy
      }
      matched <- retrieveSCEIndex(vals$counts, uniqInput, 'gene',
                                  by = BY, exactMatch = input$hmGeneTextEM,
                                  firstMatch = input$hmGeneTextFM)
      nMatched <- length(matched)
    } else {
      nMatched <- 0
    }

    p(paste0(nInput, " unique input, ", nMatched, "matched."))
  })

  observeEvent(input$hmGeneAddFromText, {
    if(!is.null(vals$counts)){
      inputList <- str_trim(scan(text = input$hmGeneText,
                                 sep='\n', what = 'character', quiet = TRUE))
      uniqInput <- unique(inputList)
      if(length(uniqInput) > 0){
        if(!is.null(input$hmGeneTextBy) && input$hmGeneTextBy == 'Row Names'){
          BY <- NULL
        } else {
          BY <- input$hmGeneTextBy
        }
        newIdx <- retrieveSCEIndex(vals$counts, uniqInput, 'gene',
                                   by = BY, exactMatch = input$hmGeneTextEM,
                                   firstMatch = input$hmGeneTextFM)
        DT::selectRows(hmGeneColTable_proxy,
                       sort(unique(c(input$hmGeneColTable_rows_selected,
                                     newIdx))))
      }
    }
  })

  output$hmGeneSumUI <- renderUI({
    nGene <- length(hmTemp$geneIndex)
    if(nGene == 0){
      p("No features selected, going to use them all",
        style = 'margin-top: 5px;')
    } else {
      p(paste0("Totally ", nGene, " features selected."),
        style = 'margin-top: 5px;')
    }
  })

  # Heatmap: Annotation color assignment ####

  output$hmCellAnnUI <- renderUI({
    if(!is.null(vals$counts)){
      classes <- names(colData(hmTemp$sce))
      selectInput('hmCellAnn', 'Add cell annotation', classes,
                  multiple = TRUE, selected = hmTemp$colDataName)
    }
  })

  output$hmGeneAnnUI <- renderUI({
    if(!is.null(vals$counts)){
      classes <- names(rowData(hmTemp$sce))
      selectInput('hmGeneAnn', 'Add feature annotation', classes,
                  multiple = TRUE, selected = hmTemp$rowDataName)
    }
  })

  observeEvent(input$hmCellAnn, {
    hmTemp$colDataName <- input$hmCellAnn
  })
  observeEvent(input$hmGeneAnn, {
    hmTemp$rowDataName <- input$hmGeneAnn
  })

  hmAnnAllColors <- reactiveValues(
    col = NULL,
    row = NULL
  )
  observeEvent(vals$counts, {
    if(!is.null(vals$counts)){
      hmAnnAllColors$col <- singleCellTK:::dataAnnotationColor(hmTemp$sce, 'col')
      hmAnnAllColors$row <- singleCellTK:::dataAnnotationColor(hmTemp$sce, 'row')
    }
  })

  generateAnnColAssUI <- function(colname, axis){
    if(axis == "row"){
      data <- as.vector(rowData(hmTemp$sce)[[colname]])
    } else if(axis == 'col'){
      data <- as.vector(colData(hmTemp$sce)[[colname]])
    }
    nUniq <- length(as.vector(unique(data[!is.na(data)])))
    if(colname %in% names(hmTemp[[paste0(axis, "ColorPresets")]])){
      cats = names(hmTemp[[paste0(axis, "ColorPresets")]][[colname]])
      fluidRow(style = "padding-left:20px;",
               h4(colname),
               lapply(seq_along(cats), function(i) {
                 column(
                   width = 3,
                   colourpicker::colourInput(
                     inputId = paste0('hm', axis, colname, cats[i]),
                     label = cats[i],
                     value = hmTemp[[paste0(axis, "ColorPresets")]][[colname]][[cats[i]]]
                   )
                 )
               })
      )
    } else if(nUniq > 12){
      if(is.numeric(data)){
        fluidRow(style = "padding-left:20px;",
                 h4(colname),
                 p(paste0("Numeric annotation with ", nUniq, " unique values detected. Please choose the type of legend.")),
                 radioButtons(
                   inputId = paste0('hm', axis, colname, 'type'),
                   label = NULL,
                   choices = c('Categorical', 'Continuous'),
                   inline = TRUE
                 ),
                 conditionalPanel(
                   condition = paste0("input.hm", axis, colname, "type == 'Categorical'"),
                   p("Since more than 12 unique values detected, discrete colors will be assigned for this class")
                 ),
                 conditionalPanel(
                   condition = paste0("input.hm", axis, colname, "type == 'Continuous'"),
                   p("We generate a gradient color legend for continuous annotation value"),
                   column(
                     width = 6,
                     colourpicker::colourInput(
                       inputId = paste0('hm', axis, colname, 'High'),
                       label = 'High Value'
                     )
                   ),
                   column(
                     width = 6,
                     colourpicker::colourInput(
                       inputId = paste0('hm', axis, colname, 'Low'),
                       label = 'Low Value'
                     )
                   )
                 ),
        )
      } else {
        fluidRow(style = "padding-left:20px;", h4(colname),
                 p(paste0("Totally ", nUniq, " unique values in this class of annotation, which is too many to provide manual selection. Coloring will be provided by default."))
        )
      }

    } else if(nUniq >= 1 && nUniq <= 12){
      cats <- as.character(unique(data))
      fluidRow(style = "padding-left:20px;",
               h4(colname),
               lapply(seq_along(cats), function(i) {
                 if(!is.na(cats[i])){
                   column(
                     width = 3,
                     colourpicker::colourInput(
                       inputId = paste0('hm', axis, colname, cats[i]),
                       label = cats[i],
                       value = hmAnnAllColors[[axis]][[colname]][[cats[i]]]
                     )
                   )
                 } else {
                   column(
                     width = 3,
                     colourpicker::colourInput(
                       inputId = paste0('hm', axis, colname, cats[i]),
                       label = "NA",
                       value = #FFFFFF
                     )
                   )
                 }
               })
      )
    } else {
      fluidRow(style = "padding-left:20px;",
               h4(colname),
               p("No effective category found for the class.")
      )
    }
  }

  observeEvent(input$hmCellAnn, {
    if(!is.null(input$hmCellAnn)){
      output$hmCellAnnAssUI <- renderUI({
        panel(
          lapply(input$hmCellAnn, generateAnnColAssUI, axis = 'col')
        )
      })
    }
  })

  observeEvent(input$hmGeneAnn, {
    if(!is.null(input$hmGeneAnn)){
      output$hmGeneAnnAssUI <- renderUI({
        panel(
          lapply(input$hmGeneAnn, generateAnnColAssUI, axis = 'row')
        )
      })
    }
  })

  observe({
    for (i in names(hmTemp$colColorPresets)){
      if (i %in% hmTemp$colDataName){
        for (j in names(hmTemp$colColorPresets[[i]])){
          if(!is.null(input[[paste0('hmcol', i, j)]])){
            hmTemp$colColorPresets[[i]][[j]] <- input[[paste0('hmcol', i, j)]]
          }
        }
      }
    }
  })
  observe({
    for (i in names(hmTemp$rowColorPresets)){
      if (i %in% hmTemp$rowDataName){
        for (j in names(hmTemp$rowColorPresets[[i]])){
          if(!is.null(input[[paste0('hmrow', i, j)]])){
            hmTemp$rowColorPresets[[i]][[j]] <- input[[paste0('hmrow', i, j)]]
          }
        }
      }
    }
  })

  # Heatmap: Others ####
  output$hmColSplitUI <- renderUI({
    selectInput(
      'hmColSplit',
      "Split columns (cell) by (Leave this for not splitting)",
      hmTemp$colDataName, multiple = TRUE, selected = hmTemp$colSplitBy
    )
  })

  output$hmRowSplitUI <- renderUI({
    selectInput(
      'hmRowSplit',
      "Split rows (feature) by (Leave this for not splitting)",
      hmTemp$rowDataName, multiple = TRUE, selected = hmTemp$rowSplitBy
    )
  })
  observeEvent(input$hmColSplit, {
    hmTemp$colSplitBy <- input$hmColSplit
  })
  observeEvent(input$hmRowSplit, {
    hmTemp$rowSplitBy <- input$hmRowSplit
  })

  output$hmTrimUI <- renderUI({
    if(!is.null(vals$counts) && !is.null(input$hmAssay)){
      # This might be slow when running with real data
      mat <- as.matrix(assay(vals$counts, input$hmAssay))
      if(isTRUE(input$hmScale)){
        mat <- as.matrix(computeZScore(mat))
      }
      sliderInput("hmTrim",  "Trim", min = floor(min(mat)),
                  max = ceiling(max(mat)), value = c(-2, 2), step = 0.5)
    }
  })

  # Heatmap: Color Scheme ####
  observe({
    # Palette preset coding refers:
    # https://stackoverflow.com/a/52552008/13676674
    vals$hmCSURL <- session$registerDataObj(
      name = 'uniquename1',
      data = vals$hmCSPresets,
      filter = function(data, req) {
        query <- parseQueryString(req$QUERY_STRING)
        palette <- query$palette
        cols <- data[[palette]]
        image <- tempfile()
        tryCatch({
          png(image, width = 75, height = 25, bg = 'transparent')
          par(mar = c(0, 0, 0, 0))
          barplot(rep(1, length(cols)), col = cols, axes = FALSE)
        },finally = dev.off())

        shiny:::httpResponse(
          200, 'image/png',readBin(image, 'raw', file.info(image)[,'size'])
        )
      }
    )

    updateSelectizeInput(
      session, 'hmCSPalette', server = TRUE,
      choices = names(vals$hmCSPresets),
      selected = "RWB",
      options = list(
        render = I(
          sprintf(
            "{
            option: function(item, escape) {
            return '<div><img width=\"75\" height=\"25\" ' +
            'src=\"%s&palette=' + escape(item.value) + '\" />' +
            escape(item.value) + '</div>';
            }
          }",
            vals$hmCSURL
          )
        )
      )
    )
  })

  observeEvent(input$hmCSPalette, {
    if(!input$hmCSPalette == ""){
      lowColor <- vals$hmCSPresets[[input$hmCSPalette]][1]
      colourpicker::updateColourInput(session, 'hmCSLow', value = lowColor)
    }
    if(!input$hmCSPalette == ""){
      mediumColor <- vals$hmCSPresets[[input$hmCSPalette]][2]
      colourpicker::updateColourInput(session, 'hmCSMedium', value = mediumColor)
    }
    if(!input$hmCSPalette == ""){
      highColor <- vals$hmCSPresets[[input$hmCSPalette]][3]
      colourpicker::updateColourInput(session, 'hmCSHigh', value = highColor)
    }
  })

  # Heatmap: Final run ####
  observeEvent(input$plotHeatmap, {
    if (is.null(vals$counts)){
      shinyalert::shinyalert("Error!", "Upload data first.", type = "error")
    } else {
      # Move all plotting process into alert callback, thus auto-re-render can
      # be avoided while tuning parameters.
      shinyalert(
        title = "Confirm",
        text = "Large dataset might take time to rerun. Are you sure with the parameters?",
        type = "warning",
        showCancelButton = TRUE,
        confirmButtonText = "Plot",
        cancelButtonText = "Check Once More",
        callbackR = function(x){
          if(isTRUE(x)){
            withBusyIndicatorServer("plotHeatmap", {
              if(!is.null(hmTemp$colDataName)){
                cellAnnColor <- list()
                for(i in hmTemp$colDataName){
                  uniqs <- as.vector(unique(colData(hmTemp$sce)[[i]]))
                  uniqs[is.na(uniqs)] <- 'NA'
                  if (i %in% names(hmTemp$colColorPresets)) {
                    cellAnnColor[[i]] <- hmTemp$colColorPresets[[i]]
                  } else if (length(uniqs) <= 12) {
                    cellAnnColor[[i]] <- vector()
                    for(j in uniqs){
                      inputId <- paste0('hmcol', i, j)
                      cellAnnColor[[i]] <- c(cellAnnColor[[i]], input[[inputId]])
                    }
                    names(cellAnnColor[[i]]) <- uniqs
                  } else {
                    if(is.numeric(colData(hmTemp$sce)[[i]])){
                      if(input[[paste0('hmcol', i, 'type')]] == 'Continuous'){
                        cFun <- circlize::colorRamp2(
                          c(min(colData(hmTemp$sce)[[i]]),
                            max(colData(hmTemp$sce)[[i]])),
                          c(input[[paste0('hmcol', i, 'Low')]],
                            input[[paste0('hmcol', i, 'High')]])
                        )
                        cellAnnColor[[i]] <- cFun
                      } else {
                        c <- distinctColors(length(uniqs))
                        names(c) <- uniqs
                        cellAnnColor[[i]] <- c
                      }
                    }
                  }
                }
              } else {
                cellAnnColor <- NULL
              }
              if(!is.null(hmTemp$rowDataName)){
                geneAnnColor <- list()
                for(i in hmTemp$rowDataName){
                  uniqs <- as.vector(unique(rowData(hmTemp$sce)[[i]]))
                  if (i %in% names(hmTemp$rowColorPresets)) {
                    geneAnnColor[[i]] <- hmTemp$rowColorPresets[[i]]
                  } else if(length(uniqs) <= 12){
                    geneAnnColor[[i]] <- vector()
                    for(j in uniqs){
                      inputId <- paste0('hmrow', i, j)
                      geneAnnColor[[i]] <- c(geneAnnColor[[i]], input[[inputId]])
                    }
                    names(geneAnnColor[[i]]) <- uniqs
                  } else {
                    if(is.numeric(rowData(hmTemp$sce)[[i]])){
                      if(input[[paste0('hmrow', i, 'type')]] == 'Continuous'){
                        cFun <- circlize::colorRamp2(
                          c(min(rowData(hmTemp$sce)[[i]]),
                            max(rowData(hmTemp$sce)[[i]])),
                          c(input[[paste0('hmrow', i, 'Low')]],
                            input[[paste0('hmrow', i, 'High')]])
                        )
                        geneAnnColor[[i]] <- cFun
                      } else {
                        c <- distinctColors(length(uniqs))
                        names(c) <- uniqs
                        geneAnnColor[[i]] <- c
                      }
                    }
                  }
                }
              } else {
                geneAnnColor <- NULL
              }
              hmAddLabel <- list(cell = FALSE, gene = FALSE)
              if(!is.null(input$hmAddLabel)){
                if("1" %in% input$hmAddLabel){
                  if(input$hmAddCellLabel == "Default cell IDs"){
                    hmAddLabel$cell <- TRUE
                  } else {
                    hmAddLabel$cell <- input$hmAddCellLabel
                  }
                }
                if("2" %in% input$hmAddLabel){
                  if(input$hmAddGeneLabel == "Default feature IDs"){
                    hmAddLabel$gene <- TRUE
                  } else {
                    hmAddLabel$gene <- input$hmAddGeneLabel
                  }
                }
              }
              hmShowDendro <- c(FALSE, FALSE)
              hmShowDendro[as.numeric(input$hmShowDendro)] <- TRUE
              #if(is.null(hmTemp$rowSplitBy)){
              #  hmRowSplit <- NULL
              #} else {
              #  hmRowSplit <- hmTemp$rowSplitBy
              #}
              #if(is.null(hmTemp$colSplitBy)){
              #  hmColSplit <- NULL
              #} else {
              #  hmColSplit <- hmTemp$colSplitBy
              #}
              cs <- circlize::colorRamp2(
                c(input$hmTrim[1], mean(input$hmTrim), input$hmTrim[2]),
                c(input$hmCSLow, input$hmCSMedium, input$hmCSHigh)
              )
              output$Heatmap <- renderPlot({
                isolate({
                  plotSCEHeatmap(
                    inSCE = hmTemp$sce, useAssay = input$hmAssay, colorScheme = cs,
                    featureIndex = hmTemp$geneIndex, cellIndex = hmTemp$cellIndex,
                    rowDataName = hmTemp$rowDataName, colDataName = hmTemp$colDataName,
                    rowSplitBy = hmTemp$rowSplitBy, colSplitBy = hmTemp$colSplitBy,
                    rowLabel = hmAddLabel$gene, colLabel = hmAddLabel$cell,
                    rowDend = hmShowDendro[2], colDend = hmShowDendro[1],
                    scale = input$hmScale, trim = input$hmTrim,
                    width = unit(20, 'cm'), height = unit(20, 'cm'),
                    featureAnnotationColor = geneAnnColor,
                    cellAnnotationColor = cellAnnColor
                  )
                })
              }, height = 800)
            })
          }
        }
      )
    }
  })

  #-----------------------------------------------------------------------------
  # Page 4: Batch Correction ####
  #-----------------------------------------------------------------------------

  observeEvent(input$closeDropDownBC, {
    session$sendCustomMessage("close_dropDownBC", "")
  })
  
  observeEvent(input$batchCorrMethods, {
    if (!is.null(vals$counts) &&
        !is.null(input$batchCorrMethods)) {
      # What type of assays are required, according their docs
      # ComBatSeq - counts
      # BBKNN - filtered, normalized, and scaled
      # fastMNN - log-expression
      # Limma - log-expression
      # MNN - log-expression
      # scanorama - normalized, log1p
      # scMerge - logcounts
      # zinbwave - counts
      bc.recommended <- NULL
      method.log <- c("FastMNN", "Limma", "MNN")
      method.scale <- c("BBKNN")
      method.raw <- c("ZINBWaVE", "ComBatSeq")
      if (is.null(input$batchCorrMethods)) {
        bc.recommended <- "raw"
      } else if (input$batchCorrMethods %in% method.log) {
        bc.recommended <- c("transformed", "normalized")
      } else if (input$batchCorrMethods %in% method.raw) {
        bc.recommended <- "raw"
      } else if (input$batchCorrMethods %in% method.scale) {
        bc.recommended <- "scaled"
      }
      updateSelectInputTag(session, "batchCorrAssay",
                           label = "Select Assay to Correct:",
                           choices = assayNames(vals$counts),
                           recommended = bc.recommended)
    }
  })

  output$batchCheckResUI <- renderUI({
    selectInput("batchCheckCorrName", "Corrected Matrix",
                c(names(vals$batchRes)))
  })

  observeEvent(input$plotBatchCheck, {
    if(!is.null(vals$counts) &&
       !is.null(input$batchCheckCorrName) &&
       input$batchCheckVar != input$batchCheckCond){
      withBusyIndicatorServer("plotBatchCheck", {
        ## Generals
        if(input$batchCheckCond == "None"){
          shapeBy <- NULL
        } else {
          shapeBy <- input$batchCheckCond
        }
        ## Original assay PCA
        oriAssayPCAName <- paste0(input$batchCheckOrigAssay, "_PCA")
        if(!oriAssayPCAName %in% names(reducedDims(vals$counts))){
          vals$counts <- scaterPCA(vals$counts,
                                useAssay = input$batchCheckOrigAssay,
                                reducedDimName = oriAssayPCAName)
          updateReddimInputs()
        }
        resName <- input$batchCheckCorrName
        ## Corrected assay/altExp PCA
        if (vals$batchRes[[resName]] == 'assay'){
          corrAssayPCAName = paste0(resName, "_PCA")
          vals$counts <- scaterPCA(vals$counts, useAssay = resName,
                                reducedDimName = corrAssayPCAName)
          updateReddimInputs()
        } else if (vals$batchRes[[resName]] == 'altExp'){
          ae <- altExp(vals$counts, resName)
          corrAltExpPCAName <- paste0(resName, "_PCA")
          ae <- scaterPCA(ae, useAssay = resName,
                       reducedDimName = corrAltExpPCAName)
          reducedDim(vals$counts, corrAltExpPCAName) <-
            reducedDim(ae, corrAltExpPCAName)
          updateReddimInputs()
        }
        ## Update plots
        output$batchOriVars <- renderPlot({
          isolate({
            plotBatchVariance(inSCE = vals$counts,
                              useAssay = input$batchCheckOrigAssay,
                              batch = input$batchCheckVar,
                              condition = shapeBy)
          })
        })
        output$batchOriPCA <- renderPlot({
          isolate({
            plotSCEDimReduceColData(vals$counts, colorBy = input$batchCheckVar,
                                    shape = shapeBy,
                                    reducedDimName = oriAssayPCAName,
                                    dim1 = 1, dim2 = 2,
                                    title = paste0("Original ",
                                                   input$batchCheckOrigAssay,
                                                   " PCA"))
          })
        })
        output$batchCorrVars <- renderPlot({
          isolate({
            if (vals$batchRes[[resName]] == 'reddim'){
              plotBatchVariance(inSCE = vals$counts, useReddim = resName,
                                batch = input$batchCheckVar,
                                condition = shapeBy)
            } else if (vals$batchRes[[resName]] == 'assay'){
              plotBatchVariance(inSCE = vals$counts, useAssay = resName,
                                batch = input$batchCheckVar,
                                condition = shapeBy)
            } else if (vals$batchRes[[resName]] == 'altExp'){
              plotBatchVariance(inSCE = vals$counts, useAltExp = resName,
                                batch = input$batchCheckVar,
                                condition = shapeBy)
            }
          })
        })
        output$batchCorrReddim <- renderPlot({
          isolate({
            if (vals$batchRes[[resName]] == 'reddim'){
              plotSCEDimReduceColData(vals$counts,
                                      colorBy = input$batchCheckVar,
                                      shape = shapeBy,
                                      reducedDimName = resName,
                                      conditionClass = "character",
                                      dim1 = 1, dim2 = 2,
                                      title = paste0(resName, " corrected"))
            } else if (vals$batchRes[[resName]] == 'assay'){
              plotSCEDimReduceColData(vals$counts,
                                      colorBy = input$batchCheckVar,
                                      shape = shapeBy,
                                      reducedDimName = corrAssayPCAName,
                                      conditionClass = "character",
                                      dim1 = 1, dim2 = 2,
                                      title = paste0(resName, " corrected"))
            } else if (vals$batchRes[[resName]] == 'altExp'){
              plotSCEDimReduceColData(vals$counts,
                                      colorBy = input$batchCheckVar,
                                      shape = shapeBy,
                                      reducedDimName = corrAltExpPCAName,
                                      dim1 = 1, dim2 = 2,
                                      title = paste0(resName, " corrected"))
            }
          })
        })
      })
    }
    session$sendCustomMessage("close_dropDownBC", "")
  })

  observeEvent(input$BBKNNRun, {
    if (is.null(vals$counts)){
      shinyalert::shinyalert("Error!", "Upload data first.", type = "error")
    } else {
      withBusyIndicatorServer("BBKNNRun", {
        saveassayname <- gsub(" ", "_", input$BBKNNSaveReddim)
        vals$counts <- runBBKNN(vals$counts,
                                useAssay = input$batchCorrAssay,
                                batch = input$batchCorrVar,
                                reducedDimName = saveassayname,
                                nComponents = input$BBKNNNComp)
        shinyalert::shinyalert('Success!', 'BBKNN completed.',
                               type = 'success')
        vals$batchRes[[saveassayname]] <- 'reddim'
        updateReddimInputs()
      })
    }
  })

  observeEvent(input$combatRun, {
    if (is.null(vals$counts)){
      shinyalert::shinyalert("Error!", "Upload data first.", type = "error")
    }
    else{
      withBusyIndicatorServer("combatRun", {
        #check for zeros
        if (any(rowSums(assay(vals$counts, input$batchCorrAssay)) == 0)){
          shinyalert::shinyalert("Error!", "Rows with a sum of zero found. Filter data to continue.", type = "error")
        } else {
          saveassayname <- gsub(" ", "_", input$combatSaveAssay)
          if (input$combatKnownCT == "Yes") {
            cov <- input$combatCond
          } else {
            cov <- NULL
          }
          if (input$combatCTBalance == "Yes") {
            useSVA <- FALSE
          } else {
            useSVA <- TRUE
          }
          if (input$combatBioCond == "None") {
            combatBioCond <- NULL
          } else {
            combatBioCond <- input$combatBioCond
          }
          vals$counts <- runComBatSeq(inSCE = vals$counts,
                                      useAssay = input$batchCorrAssay,
                                      batch = input$batchCorrVar,
                                      covariates = cov,
                                      bioCond = combatBioCond,
                                      useSVA = useSVA,
                                      assayName = saveassayname,
                                      shrink = input$combatShrink,
                                      shrinkDisp = input$combatShrinkDisp,
                                      nGene = input$combatNGene                                      )
          vals$batchRes[[saveassayname]] <- 'assay'
          updateAssayInputs()
          shinyalert::shinyalert('Success!', 'ComBatSeq completed.',
                                 type = 'success')
        }
      })
    }
  })

  observeEvent(input$FastMNNRun, {
    if (is.null(vals$counts)){
      shinyalert::shinyalert("Error!", "Upload data first.", type = "error")
    } else {
      withBusyIndicatorServer("FastMNNRun", {
        saveassayname <- gsub(" ", "_", input$FastMNNSaveReddim)
        if(isTRUE(input$FastMNNPcInput)){
          fmnnAssay <- input$FastMNNReddim
        } else {
          fmnnAssay <- input$batchCorrAssay
        }
        vals$counts <- runFastMNN(vals$counts,
                                  useAssay = fmnnAssay,
                                  batch = input$batchCorrVar,
                                  reducedDimName = saveassayname,
                                  pcInput = input$FastMNNPcInput
        )
        shinyalert::shinyalert('Success!', 'FastMNN completed.',
                               type = 'success')
        vals$batchRes[[saveassayname]] <- 'reddim'
        updateReddimInputs()
      })
    }
  })

  # observeEvent(input$HarmonyRun, {
  #   if (is.null(vals$counts)){
  #     shinyalert::shinyalert("Error!", "Upload data first.", type = "error")
  #   } else {
  #     withBusyIndicatorServer("HarmonyRun", {
  #       saveassayname <- gsub(" ", "_", input$HarmonySaveReddim)
  #       if(isTRUE(input$HarmonyPcInput)){
  #         useAssay <- input$HarmonyReddim
  #       } else {
  #         useAssay <- input$batchCorrAssay
  #       }
  #       if(is.na(as.numeric(input$HarmonyTheta))){
  #         stop("Theta value must be numeric.")
  #       } else {
  #         theta <- as.numeric(input$HarmonyTheta)
  #       }
  #       vals$counts <- runHarmony(vals$counts, useAssay = useAssay,
  #                                 pcInput = input$HarmonyPcInput,
  #                                 batch = input$batchCorrVar,
  #                                 reducedDimName = saveassayname,
  #                                 nComponents = input$HarmonyNComp,
  #                                 theta = theta, nIter = input$HarmonyNIter)
  #       shinyalert::shinyalert('Success!', 'Harmony completed.',
  #                              type = 'success')
  #       vals$batchRes[[saveassayname]] <- 'reddim'
  #       updateReddimInputs()
  #     })
  #   }
  # })

  observeEvent(input$limmaRun, {
    if (is.null(vals$counts)){
      shinyalert::shinyalert("Error!", "Upload data first.", type = "error")
    } else {
      withBusyIndicatorServer("limmaRun", {
        saveassayname <- gsub(" ", "_", input$limmaSaveAssay)
        vals$counts <- runLimmaBC(vals$counts,
                                  useAssay = input$batchCorrAssay,
                                  batch = input$batchCorrVar,
                                  assayName = saveassayname)
        shinyalert::shinyalert('Success!', 'Limma completed.',
                               type = 'success')
        vals$batchRes[[saveassayname]] <- 'assay'
        updateAssayInputs()
      })
    }
  })

  observeEvent(input$ligerRun, {
    if (is.null(vals$counts)){
      shinyalert::shinyalert("Error!", "Upload data first.", type = "error")
    }
    else{
      withBusyIndicatorServer("ligerRun", {
        #check for zeros
        if (any(rowSums(assay(vals$counts, input$batchCorrAssay)) == 0)){
          shinyalert::shinyalert("Error!", "Rows with a sum of zero found. Filter data to continue.", type = "error")
        } else {
          saveassayname <- gsub(" ", "_", input$ligerSaveReddim)
          vals$counts <-
            runLIGER(inSCE = vals$counts,
                     useAssay = input$batchCorrAssay,
                     batch = input$batchCorrVar,
                     reducedDimName = saveassayname,
                     nComponents = input$ligerNComp,
                     lambda = input$ligerLambda,
                     resolution = input$ligerResolution)
          shinyalert::shinyalert('Success!', 'LIGER completed.',
                                 type = 'success')
          vals$batchRes[[saveassayname]] <- 'reddim'
          updateReddimInputs()
        }
      })
    }
  })

  observeEvent(input$MNNRun, {
    if (is.null(vals$counts)){
      shinyalert::shinyalert("Error!", "Upload data first.", type = "error")
    } else {
      withBusyIndicatorServer("MNNRun", {
        saveassayname <- gsub(" ", "_", input$MNNSaveAssay)
        vals$counts <- runMNNCorrect(vals$counts,
                                     useAssay = input$batchCorrAssay,
                                     batch = input$batchCorrVar,
                                     k = input$MNNK, sigma = input$MNNSigma,
                                     assayName = saveassayname)
        shinyalert::shinyalert('Success!', 'MNN completed.',
                               type = 'success')
        vals$batchRes[[saveassayname]] <- 'assay'
        updateAssayInputs()
      })
    }
  })

  observeEvent(input$scnrmRun, {
    if (is.null(vals$counts)){
      shinyalert::shinyalert("Error!", "Upload data first.", type = "error")
    } else {
      withBusyIndicatorServer("scnrmRun", {
        saveassayname <- gsub(" ", "_", input$scnrmSaveAssay)
        vals$counts <- runSCANORAMA(vals$counts,
                                    useAssay = input$batchCorrAssay,
                                    batch = input$batchCorrVar,
                                    SIGMA = input$scnrmSIGMA,
                                    ALPHA = input$scnrmALPHA,
                                    KNN = input$scnrmKNN,
                                    assayName = saveassayname)
        shinyalert::shinyalert('Success!', 'SCANORAMA completed.',
                               type = 'success')
        vals$batchRes[[saveassayname]] <- 'assay'
        updateAssayInputs()
      })
    }
  })

  output$scMergeNBatch <- renderUI({
    if(!is.null(vals$counts) &&
       !is.null(input$batchCorrVar)){
      nBatch <- length(unique(SummarizedExperiment::colData(vals$counts)[[input$batchCorrVar]]))
      span(paste0("Please input ", nBatch, " integer(s), separated by ','."))
    }
  })

  observeEvent(input$scMergeRun, {
    if (is.null(vals$counts)){
      shinyalert::shinyalert("Error!", "Upload data first.", type = "error")
    } else {
      withBusyIndicatorServer("scMergeRun", {
        saveassayname <- gsub(" ", "_", input$scMergeSaveAssay)
        if(input$scMergeSEGOpt == 1){
          seg <- NULL
        } else if(input$scMergeSEGOpt == 2){
          data("SEG")
          seg <- SEG[[input$scMergeSEGSpecies]]
        } else {
          seg <- str_trim(scan(text = input$scMergeSEGCustom,
                               sep='\n', what = 'character'))
        }
        if(isTRUE(input$scMergeAutoKmk)){
          kmk <- NULL
        } else {
          kmk <- scan(text = input$scMergeUserKmk, sep=',')
        }
        vals$counts <- runSCMerge(inSCE = vals$counts,
                                  useAssay = input$batchCorrAssay,
                                  batch = input$batchCorrVar,
                                  cellType = input$scMergeCT,
                                  seg = seg, kmeansK = kmk,
                                  assayName = saveassayname
        )
        shinyalert::shinyalert('Success!', 'scMerge completed.',
                               type = 'success')
        vals$batchRes[[saveassayname]] <- 'assay'
        updateAssayInputs()
      })
    }
  })

  output$Srt3IntNAnchUI <- renderUI({
    if(!is.null(vals$counts)){
      ngene <- nrow(vals$counts)
      tagList(

        numericInput('Srt3IntNAnch', "Number of anchors:",
                     value = ngene, min = 30, max = ngene, step = 1),

        numericInput('Srt3IntKWeight', "kWeight:",
                     value = 0, min = 0, step = 1),

        numericInput('Srt3IntKFilter', "kFilter:",
                     value = 0, min = 0, step = 1),

        numericInput('Srt3IntNDims', "Number of Dimensions:",
                     value = 2, min = 2, step = 1)

      )
    }
  })

  observeEvent(input$Srt3IntRun, {
    if (is.null(vals$counts)){
      shinyalert::shinyalert("Error!", "Upload data first.", type = "error")
    } else {
      withBusyIndicatorServer("Srt3IntRun", {
        saveassayname <- gsub(" ", "_", input$Srt3IntSaveAssay)
        vals$counts <- seuratIntegration(
          inSCE = vals$counts,
          batch = input$batchCorrVar,
          newAssayName = saveassayname,
          kAnchor = input$Srt3IntNAnch,
          kWeight = input$Srt3IntKWeight,
          kFilter = input$Srt3IntKFilter,
          ndims = input$Srt3IntNDims)
        vals$batchRes[[saveassayname]] <- 'altExp'
        shinyalert::shinyalert('Success!', 'Seurat3 Integration completed.',
                               type = 'success')
      })
    }
  })

  output$zinbwaveNHvgUI <- renderUI({
    if(!is.null(vals$counts)){
      ngenes <- nrow(vals$counts)
      zwdefault <- min(ngenes, 1000)
      numericInput('zinbwaveNHVG', 'Number of highly variable genes to use:',
                   value = zwdefault, max = ngenes)
    }
  })

  output$zinbwaveEpsUI <- renderUI({
    if(!is.null(vals$counts)){
      ngenes <- nrow(vals$counts)
      zwdefault <- min(ngenes, 1000)
      numericInput('zinbwaveEps', 'Epsilon value:',
                   value = zwdefault, max = ngenes)
    }
  })

  observeEvent(input$zinbwaveRun, {
    if (is.null(vals$counts)){
      shinyalert::shinyalert("Error!", "Upload data first.", type = "error")
    } else {
      withBusyIndicatorServer("zinbwaveRun", {
        saveassayname <- gsub(" ", "_", input$limmaSaveAssay)
        vals$counts <- runZINBWaVE(vals$counts,
                                   useAssay = input$batchCorrAssay,
                                   batch = input$batchCorrVar,
                                   reducedDimName = saveassayname,
                                   epsilon = input$zinbwaveEps,
                                   nHVG = input$zinbwaveNHVG,
                                   nIter = input$zinbwaveNIter,
                                   nComponents = input$zinbwaveNComp
        )
        shinyalert::shinyalert('Success!', 'ZINBWaVE completed.',
                               type = 'success')
        vals$batchRes[[saveassayname]] <- 'reddim'
        updateReddimInputs()
      })
    }
  })

  #-----------------------------------------------------------------------------
  # Page 4.1: Feature Selection ####
  #-----------------------------------------------------------------------------

  observeEvent(input$findHvgButtonFS, {
    withBusyIndicatorServer("findHvgButtonFS", {
      if (!is.null(vals$counts)) {
        tryCatch(vals$counts <- runFeatureSelection(
          inSCE = vals$counts,
          useAssay = input$assaySelectFS_Norm,
          hvgMethod = input$hvgMethodFS
        ), error = function(e) stop("HVG computation failed. Try re-computing with a normalized assay!"))
        vals$hvgCalculated$status <- TRUE
        vals$hvgCalculated$method <- input$hvgMethodFS
        vals$hvgCalculated$assayName <- input$assaySelectFS_Norm
      }
    })
    updateSelectInputTag(session, "hvgSubsetAssay",
                         recommended = c("scaled", "transformed", "normalized"),
                         choices = expDataNames(vals$counts))
  })

  observeEvent(input$updatePlotFS, {
    req(vals$counts)
      if (isTRUE(vals$hvgCalculated$status) &&
          !is.null(vals$hvgCalculated$method)) {
        #checks
        if(is.na(input$hvgNoFeaturesViewFS)){
          stop("Number of features cannot be empty!")
        }
        #processing
        HVGs <- getTopHVG(inSCE = vals$counts,
                          method = input$hvgMethodFS,
                          n = input$hvgNoFeaturesViewFS)
        vals$vfplot <- plotTopHVG(
          inSCE =  vals$counts,
          method = input$hvgMethodFS,
          hvgList = HVGs
        )
        output$plotFS <- renderPlot({
          isolate({
            if (!is.null(vals$vfplot)) {
              vals$vfplot
            }
          })
        })
        output$hvgOutputFS <- renderText({
          isolate({
            HVGs
          })
        })
      }
    session$sendCustomMessage("close_dropDownFS", "")
  })

  observeEvent(input$closeDropDownFS, {
    req(vals$counts)
    session$sendCustomMessage("close_dropDownFS", "")
  })

  observeEvent(input$hvgSubsetRun, {
    withBusyIndicatorServer("hvgSubsetRun", {
      if (isTRUE(vals$hvgCalculated$status) &&
          !is.null(vals$hvgCalculated$method) &&
          !is.null(input$hvgSubsetAssay)) {
        if(is.na(input$hvgNumberSelect)){
          stop("Number of HVG cannot be empty!")
        }
        if(input$hvgAltExpName == ""){
          stop("Name of the subset cannot be empty!")
        }
        if (input$hvgAltExpName %in% altExpNames(vals$counts)) {
          shinyalert(
            "Warning",
            "Entered subset name is already there.",
            "warning",
            showCancelButton = TRUE,
            confirmButtonText = "Overwrite",
            callbackR = function(x){
              if (isTRUE(x)) {
                HVGs <- getTopHVG(inSCE = vals$counts,
                                  method = input$hvgMethodFS,
                                  n = input$hvgNumberSelect)

                #make sure no NA's are introduced in HVGs
                HVGs <- stats::na.omit(HVGs)
                #tempAssay <- expData(vals$counts, vals$hvgCalculated$assayName)[HVGs,]
                tempAssay <- expData(vals$counts, input$hvgSubsetAssay)[HVGs,]
                expData(vals$counts, input$hvgAltExpName, tag = "hvg", altExp = TRUE) <- tempAssay
                updateAssayInputs()
              }
            })
        } else {
          HVGs <- getTopHVG(inSCE = vals$counts,
                            method = input$hvgMethodFS,
                            n = input$hvgNumberSelect)

          #make sure no NA's are introduced in HVGs
          HVGs <- stats::na.omit(HVGs)
          #tempAssay <- expData(vals$counts, vals$hvgCalculated$assayName)[HVGs,]
          tempAssay <- expData(vals$counts, input$hvgSubsetAssay)[HVGs,]
          expData(vals$counts, input$hvgAltExpName, tag = "hvg", altExp = TRUE) <- tempAssay
          updateAssayInputs()

          # added this
          HVGs <- getTopHVG(inSCE = vals$counts,
                            method = input$hvgMethodFS,
                            n = 100)
          vals$vfplot <- plotTopHVG(
            inSCE =  vals$counts,
            method = input$hvgMethodFS,
            hvgList = HVGs
          )
          output$plotFS <- renderPlot({
            isolate({
              if (!is.null(vals$vfplot)) {
                vals$vfplot
              }
            })
          })
          output$hvgOutputFS <- renderText({
            isolate({
              HVGs
            })
          })
        }
        # Show downstream analysis options
        callModule(module = nonLinearWorkflow, id = "nlw-fs", parent = session, dr = TRUE, cl = TRUE)
      } else {
        shinyalert::shinyalert(
          "Error",
          text = "Please compute the variance before the subsetting!",
          type = "error"
        )
      }
    })
  })

  # observeEvent(input$scatterFSRun,{
  #   useAssay <- "tophat_counts"
  #   xname <- input$scatterFSGenes[1]
  #   yname <- input$scatterFSGenes[2]
  #   xvec <- expData(vals$counts, useAssay)[xname,]
  #   yvec <- expData(vals$counts, useAssay)[yname,]
  #   customMat <- matrix(c(xvec, yvec), nrow = length(xvec))
  #   colnames(customMat) <- c(xname, yname)
  #   rownames(customMat) <- names(xvec)
  #   reducedDim(vals$counts, "Custom") <- customMat
  #   reducedDimName <- "Custom"
  #   colorLow <- "#FFFFFF"
  #   colorMid <- "#666666"
  #   colorHigh <- "#0000FF"
  #
  #
  #   a <- plotSCEDimReduceFeatures(vals$counts, feature = xname,
  #                                 reducedDimName = reducedDimName, useAssay = useAssay,
  #                                 xlab = xname, ylab = yname, transparency = 1,
  #                                 colorLow = colorLow, colorMid = colorMid, colorHigh = colorHigh,
  #                                 combinePlot = "none")
  #
  #
  #   a <- a + ggplot2::theme_bw()
  #   a <- plotly::ggplotly(a)
  #
  #   output$scatterFS <- renderPlotly({
  #     plotly::subplot(plotlist = a, titleX = TRUE, titleY = TRUE)
  #   })
  # })

  #-----------------------------------------------------------------------------
  # Page 5.1: Differential Expression ####
  #-----------------------------------------------------------------------------
  observeEvent(input$deMethod, {
    if (!is.null(vals$counts)) {
      if (is.null(input$deMethod)) {
        updateSelectInputTag(session, "deAssay", tags = c("raw", "transformed", "uncategorized", "normalized", "scaled"), recommended = c("transformed", "normalized"))
      } else if (input$deMethod == "DESeq2") {
        updateSelectInputTag(session, "deAssay", tags = c("raw", "transformed", "uncategorized", "normalized", "scaled"), recommended = c("raw"))
      } else {
        updateSelectInputTag(session, "deAssay", tags = c("raw", "transformed", "uncategorized", "normalized", "scaled"), recommended = c("transformed", "normalized"))
      }
    }
  })

  ## DE - Thresholding Vis ####
  observeEvent(input$deViewThresh, {
    if (!is.null(vals$counts) &&
        !is.null(input$deAssay)) {
      shinyjs::showElement(id= "deThreshpanel")
      withProgress(message = "Plotting thresholding...", max = 1, value = 1, {
        withBusyIndicatorServer("deViewThresh", {
          # MAST style sanity check for whether logged or not
          x <- expData(vals$counts, input$deAssay)
          if (!all(floor(x) == x, na.rm = TRUE) & max(x, na.rm = TRUE) <
              100) {
            output$deSanityWarnThresh <- renderText("")
            isLogged <- TRUE
          } else {
            output$deSanityWarnThresh <- renderText("Selected assay seems not logged (MAST style sanity check). Forcing to plot by automatically applying log-transformation. ")
            isLogged <- FALSE
          }
          suppressMessages({
            thres.grob <- plotMASTThresholdGenes(inSCE = vals$counts,
                                                 useAssay = input$deAssay,
                                                 check_sanity = FALSE,
                                                 isLogged = isLogged,
                                                 doPlot = FALSE)
          })
          nSub <- tail(strsplit(thres.grob$childrenOrder, split = '-'),
                       n = 1)[[1]][3]
          plotHeight <- ceiling(as.numeric(nSub) / 4) * 240

          output$deThreshPlotDiv <- renderUI({
            div(
              style = paste0("height: ", plotHeight, "px;"),
              plotOutput("deThreshplot"))
          })
          output$deThreshplot <- renderPlot({
            grid.draw(thres.grob)
          }, height = plotHeight)
          updateActionButton(session, "deViewThresh", "Refresh")
        })
      })
    }

  })

  observeEvent(input$deHideThresh, {
    shinyjs::hideElement(id= "deThreshpanel")
    updateActionButton(session, "deViewThresh", "View Thresholding")
  })

  ## DE - condition determination method1 ####
  output$deC1G1UI <- renderUI({
    if(!is.null(vals$counts) &
       !input$deC1Class == "None"){
      classCol <- colData(vals$counts)[[input$deC1Class]]
      classChoices <- sort(as.vector(unique(classCol)))
      selectInput(inputId = "deC1G1", label = "Select Condition(s)",
                  choices = classChoices, multiple = TRUE)
    } else {
      selectInput(inputId = "deC1G1", label = "Select Condition(s)",
                  choices = NULL, multiple = TRUE)
    }
  })

  output$deC1G2UI <- renderUI({
    if(!is.null(vals$counts) &
       !input$deC1Class == "None"){
      classCol <- colData(vals$counts)[[input$deC1Class]]
      classChoices <- sort(as.vector(unique(classCol)))
      selectInput(inputId = "deC1G2", label = "Select Condition(s)",
                  choices = classChoices, multiple = TRUE)
    } else {
      selectInput(inputId = "deC1G2", label = "Select Condition(s)",
                  choices = NULL, multiple = TRUE)
    }
  })

  output$deC1G1CellCheckUI <- renderUI({
    if(!is.null(input$deC1G1) &
       length(input$deC1G1) > 0){
      g1Idx <- colData(vals$counts)[[input$deC1Class]] %in% input$deC1G1
      g1Cells <- colnames(vals$counts)[g1Idx]
      g1CellsText <- paste(g1Cells, collapse = "\n")
      textAreaInput("deC1G1CellCheck", "Cells selected:", g1CellsText,
                    height = '100px', placeholder = "Nothing selected")
    } else {
      textAreaInput("deC1G1CellCheck", "Cells selected:", NULL,
                    height = '100px', placeholder = "Nothing selected")
    }
  })

  output$deC1G2CellCheckUI <- renderUI({
    if(!is.null(input$deC1G2) &
       length(input$deC1G2) > 0){
      g2Idx <- colData(vals$counts)[[input$deC1Class]] %in% input$deC1G2
      g2Cells <- colnames(vals$counts)[g2Idx]
      g2CellsText <- paste(g2Cells, collapse = "\n")
      textAreaInput("deC1G2CellCheck", "Cells selected:", g2CellsText,
                    height = '100px',
                    placeholder = "Leave unselected for all the others.")
    } else {
      textAreaInput("deC1G2CellCheck", "Cells selected:", NULL,
                    height = '100px',
                    placeholder = "Leave unselected for all the others.")
    }
  })

  output$deC1G1NCell <- renderUI({
    if(!is.null(input$deC1G1CellCheck)){
      if(!input$deC1G1CellCheck == ""){
        cellList <- str_trim(scan(text = input$deC1G1CellCheck,
                                  sep='\n', what = 'character', quiet = TRUE))
        cellList <- unique(cellList)
        nCell <- length(which(cellList %in% colnames(vals$counts)))
      } else {
        nCell <- 0
      }
    } else {
      nCell <- 0
    }
    msg <- paste0("Totally ", nCell, " cell(s) selected.")
    span(msg, style = 'margin-left:10px')
  })

  output$deC1G2NCell <- renderUI({
    if(!is.null(input$deC1G2CellCheck)){
      if(!input$deC1G2CellCheck == ""){
        cellList <- str_trim(scan(text = input$deC1G2CellCheck,
                                  sep='\n', what = 'character', quiet = TRUE))
        cellList <- unique(cellList)
        nCell <- length(which(cellList %in% colnames(vals$counts)))
      } else {
        nCell <- 0
      }
    } else {
      nCell <- 0
    }
    msg <- paste0("Totally ", nCell, " cell(s) selected.")
    span(msg, style = 'margin-left:10px')
  })
  ## DE - condition determination method2 ####
  ## condition 1 table operation vvvv
  output$deC2G1Table <- DT::renderDataTable({
    if(!is.null(vals$counts)){
      df <- lapply(colData(vals$counts),
                   function(i){
                     if(is.character(i) && !length(unique(i)) == length(i)){
                       return(as.factor(i))
                     } else if(is.integer(i) &&
                               !length(unique(i)) == length(i)){
                       return(as.factor(i))
                     } else {
                       return(i)
                     }
                   })
      df <- data.frame(df, row.names = colnames(vals$counts))
      DT::datatable(df, filter = "top", options = list(scrollX = TRUE))
    }
  }, server = TRUE)
  deC2G1Table_proxy <- DT::dataTableProxy("deC2G1Table")

  observeEvent(input$deC2G1Col, {
    colNames <- names(colData(vals$counts))
    showIdx <- which(colNames %in% input$deC2G1Col)
    DT::showCols(deC2G1Table_proxy, showIdx, reset = TRUE)
  })

  observeEvent(input$deC2G1Table_addAll, {
    DT::selectRows(deC2G1Table_proxy,
                   sort(unique(c(input$deC2G1Table_rows_selected,
                                 input$deC2G1Table_rows_all))))
  })

  observeEvent(input$deC2G1Table_clear, {
    DT::selectRows(deC2G1Table_proxy, NULL)
  })

  output$deC2G1info <- renderUI({
    nCell <- length(input$deC2G1Table_rows_selected)
    p(paste0("Totally ", nCell, " cells selected for ", input$deG1Name))
  })
  ## condition 1 table operation ^^^^
  ## condition 2 table operation vvvv
  output$deC2G2Table <- DT::renderDataTable({
    if(!is.null(vals$counts)){
      df <- lapply(colData(vals$counts),
                   function(i){
                     if(is.character(i) && !length(unique(i)) == length(i)){
                       return(as.factor(i))
                     } else if(is.integer(i) &&
                               !length(unique(i)) == length(i)){
                       return(as.factor(i))
                     } else {
                       return(i)
                     }
                   })
      df <- data.frame(df, row.names = colnames(vals$counts))
      DT::datatable(df, filter = "top", options = list(scrollX = TRUE))
    }
  }, server = TRUE)
  deC2G2Table_proxy <- DT::dataTableProxy("deC2G2Table")

  observeEvent(input$deC2G2Col, {
    colNames <- names(colData(vals$counts))
    showIdx <- which(colNames %in% input$deC2G2Col)
    DT::showCols(deC2G2Table_proxy, showIdx, reset = TRUE)
  })

  observeEvent(input$deC2G2Table_addAll, {
    DT::selectRows(deC2G2Table_proxy,
                   sort(unique(c(input$deC2G2Table_rows_selected,
                                 input$deC2G2Table_rows_all))))
  })

  observeEvent(input$deC2G2Table_clear, {
    DT::selectRows(deC2G2Table_proxy, NULL)
  })

  output$deC2G2info <- renderUI({
    nCell <- length(input$deC2G2Table_rows_selected)
    p(paste0("Totally ", nCell, " cells selected for ", input$deG2Name))
  })
  ## condition 2 table operation ^^^^
  ## DE - condition determination method3 ####
  output$deC3G1NCell <- renderUI({
    if(!is.null(input$deC3G1Cell)){
      if(!input$deC3G1Cell == ""){
        cellList <- str_trim(scan(text = input$deC3G1Cell,
                                  sep='\n', what = 'character', quiet = TRUE))
        cellList <- unique(cellList)
        nCell <- length(which(cellList %in% colnames(vals$counts)))
      } else {
        nCell <- 0
      }
    } else {
      nCell <- 0
    }
    msg <- paste0("Totally ", nCell, " valid cell name(s) entered.")
    span(msg, style = 'margin-left:10px')
  })

  output$deC3G2NCell <- renderUI({
    if(!is.null(input$deC3G2Cell)){
      if(!input$deC3G2Cell == ""){
        cellList <- str_trim(scan(text = input$deC3G2Cell,
                                  sep='\n', what = 'character', quiet = TRUE))
        cellList <- unique(cellList)
        nCell <- length(which(cellList %in% colnames(vals$counts)))
      } else {
        nCell <- 0
      }
    } else {
      nCell <- 0
    }
    msg <- paste0("Totally ", nCell, " valid cell name(s) entered.")
    span(msg, style = 'margin-left:10px')
  })
  # DE run analysis ####

  runDEfromShiny <- function(overwrite){
    withBusyIndicatorServer("runDE", {
      if(input$deCondMethod == 1){
        vals$counts <- runDEAnalysis(method = input$deMethod,
                                     inSCE = vals$counts,
                                     useAssay = input$deAssay,
                                     class = input$deC1Class,
                                     classGroup1 = input$deC1G1,
                                     classGroup2 = input$deC1G2,
                                     groupName1 = input$deG1Name,
                                     groupName2 = input$deG2Name,
                                     analysisName = input$deAnalysisName,
                                     covariates = input$deCovar,
                                     log2fcThreshold = input$deFCThresh,
                                     fdrThreshold = input$deFDRThresh,
                                     onlyPos = input$mastPosOnly,
                                     overwrite = overwrite)
      } else if(input$deCondMethod == 2){
        vals$counts <- runDEAnalysis(method = input$deMethod,
                                     inSCE = vals$counts,
                                     useAssay = input$deAssay,
                                     index1 = input$deC2G1Table_rows_selected,
                                     index2 = input$deC2G2Table_rows_selected,
                                     groupName1 = input$deG1Name,
                                     groupName2 = input$deG2Name,
                                     analysisName = input$deAnalysisName,
                                     covariates = input$deCovar,
                                     log2fcThreshold = input$deFCThresh,
                                     fdrThreshold = input$deFDRThresh,
                                     onlyPos = input$dePosOnly,
                                     overwrite = overwrite)
      } else {
        g1CellList <- str_trim(scan(text = input$deC3G1Cell,
                                    sep='\n', what = 'character', quiet = TRUE))
        g1CellList <- sort(unique(g1CellList))
        g2CellList <- str_trim(scan(text = input$deC3G2Cell,
                                    sep='\n', what = 'character', quiet = TRUE))
        g2CellList <- sort(unique(g2CellList))
        vals$counts <- runDEAnalysis(method = input$deMethod,
                                     inSCE = vals$counts,
                                     useAssay = input$deAssay,
                                     index1 = g1CellList,
                                     index2 = g2CellList,
                                     groupName1 = input$deG1Name,
                                     groupName2 = input$deG2Name,
                                     analysisName = input$deAnalysisName,
                                     covariates = input$deCovar,
                                     log2fcThreshold = input$deFCThresh,
                                     fdrThreshold = input$deFDRThresh,
                                     onlyPos = input$dePosOnly,
                                     overwrite = overwrite)
      }
      shinyalert::shinyalert(
        "Success",
        text = "Differential expression analysis completed.",
        type = "success"
      )
      res <- names(metadata(vals$counts)$diffExp)
      updateSelectInput(session, "deResSel", choices = res,
                        selected = input$deAnalysisName)
      colSplitBy <- "condition"
      rowSplitBy <- "regulation"

      x <- expData(vals$counts, input$deAssay)
      if (!all(floor(x) == x, na.rm = TRUE) & max(x, na.rm = TRUE) <
          100) {
        isLogged <- TRUE
      } else {
        isLogged <- FALSE
        updateCheckboxGroupInput(session, "deHMDoLog", selected = TRUE)
      }

      output$deHeatmap <- renderPlot({
        isolate({
          plotDEGHeatmap(inSCE = vals$counts,
                         useResult = input$deAnalysisName,
                         onlyPos = input$dePosOnly,
                         log2fcThreshold = input$deFCThresh,
                         fdrThreshold = input$deFDRThresh,
                         colSplitBy = colSplitBy,
                         rowSplitBy = rowSplitBy,
                         doLog = !isLogged)
        })
      })

      output$deViolinPlot <- renderPlot({
        isolate({
          plotDEGViolin(inSCE = vals$counts, useResult = input$deAnalysisName,
                        nrow = input$deVioNRow, ncol = input$deVioNCol, labelBy = NULL,
                        check_sanity = FALSE, isLogged = isLogged)
        })
      })

      output$deRegPlot <- renderPlot({
        isolate({
          plotDEGRegression(inSCE = vals$counts,
                            useResult = input$deAnalysisName,
                            nrow = input$deRegNRow,
                            ncol = input$deRegNCol,
                            labelBy = NULL,
                            check_sanity = FALSE,
                            isLogged = isLogged)
        })
      })

    })
  }

  observeEvent(input$runDE, {
    if (is.null(vals$counts)){
      shinyalert("Error!", "Upload data first.", type = "error")
    } else if(input$deAnalysisName == "" ||
              input$deG1Name == "" ||
              input$deG2Name == ""){
      shinyalert("Error!",
                 "The name of the two conditions and the whole analysis have to be specified!",
                 type = "error")
    } else {
      allRes <- names(metadata(vals$counts)$diffExp)
      if(input$deAnalysisName %in% allRes){
        shinyalert(
          "Warning",
          "Entered differential experiment analysis name is already there.",
          "warning", showCancelButton = TRUE,
          confirmButtonText = "Overwrite",
          callbackR = function(x){if(isTRUE(x)){runDEfromShiny(x)}})
      } else {
        runDEfromShiny(FALSE)
      }
      # Show downstream analysis options
      callModule(module = nonLinearWorkflow, id = "nlw-de", parent = session, pa = TRUE, cv = TRUE)
    }
  })

  # DE: Result visualize ####

  # Data table
  output$deResult <- DT::renderDataTable({
    if(!is.null(input$deResSel) &&
       !is.null(vals$counts)){
      metadata(vals$counts)$diffExp[[input$deResSel]]$result
    }
  }, filter = 'top')

  observeEvent(input$deResSel, {
    if (is.null(input$deResSel) ||
        input$deResSel == "") {
      shinyjs::disable("deDownload")
    } else {
      shinyjs::enable("deDownload")
    }
  })

  output$deDownload <- downloadHandler(
    filename = function() {
      paste0("deResult_", input$deResSel, ".csv")
    },
    content = function(file) {
      fullTable <- metadata(vals$counts)$diffExp[[input$deResSel]]$result
      filteredTable <- fullTable[input$deResult_rows_all,]
      filteredTable <- filteredTable[rowSums(is.na(filteredTable)) != ncol(filteredTable), ]
      utils::write.csv(filteredTable, file, row.names = FALSE, )
    }
  )
  # Violin plot
  output$deVioTotalUI <- renderUI({
    topN <- input$deVioNRow * input$deVioNCol
    p(as.character(topN))
  })
  
  observeEvent(input$closeDropDownDeViolin, {
    session$sendCustomMessage("close_dropDownDeViolin", "")
  })

  observeEvent(input$dePlotVio, {
    if(!is.null(input$deResSel) &&
       !input$deResSel == "" &&
       !is.null(vals$counts)){
      useAssay <- metadata(vals$counts)$diffExp[[input$deResSel]]$useAssay
      if(input$deVioLabel == "Default ID"){
        labelBy = NULL
      } else {
        labelBy = input$deVioLabel
      }
      # MAST style sanity check for whether logged or not
      x <- expData(vals$counts, useAssay)
      if (!all(floor(x) == x, na.rm = TRUE) & max(x, na.rm = TRUE) <
          100) {
        output$deSanityWarnViolin <- renderText("")
        isLogged <- TRUE
      } else {
        output$deSanityWarnViolin <- renderText("Selected assay seems not logged (MAST style sanity check). Forcing to plot by automatically applying log-transformation. ")
        isLogged <- FALSE
      }
      output$deViolinPlot <- renderPlot({
        isolate({
          plotDEGViolin(inSCE = vals$counts, useResult = input$deResSel,
                        #threshP = input$deVioUseThresh,
                        nrow = input$deVioNRow, ncol = input$deVioNCol, labelBy = labelBy,
                        check_sanity = FALSE, isLogged = isLogged)
        })
      })
      session$sendCustomMessage("close_dropDownDeViolin", "")
    }
  })
  # Linear Regression Plot
  output$deRegTotalUI <- renderUI({
    topN <- input$deRegNRow * input$deRegNCol
    p(as.character(topN))
  })
  
  observeEvent(input$closeDropDownDeReg, {
    session$sendCustomMessage("close_dropDownDeReg", "")
  })

  observeEvent(input$dePlotReg, {
    if(!is.null(input$deResSel) &&
       !input$deResSel == "" &&
       !is.null(vals$counts)){
      useAssay <- metadata(vals$counts)$diffExp[[input$deResSel]]$useAssay
      if(input$deRegLabel == "Default ID"){
        labelBy = NULL
      } else {
        labelBy = input$deRegLabel
      }
      # MAST style sanity check for whether logged or not
      x <- expData(vals$counts, useAssay)
      if (!all(floor(x) == x, na.rm = TRUE) & max(x, na.rm = TRUE) <
          100) {
        output$deSanityWarnReg <- renderText("")
        isLogged <- TRUE
      } else {
        output$deSanityWarnReg <- renderText("Selected assay seems not logged (MAST style sanity check). Forcing to plot by automatically applying log-transformation. ")
        isLogged <- FALSE
      }
      output$deRegPlot <- renderPlot({
        isolate({
          plotDEGRegression(inSCE = vals$counts,
                            useResult = input$deResSel,
                            #threshP = input$deVioUseThresh,
                            nrow = input$deRegNRow,
                            ncol = input$deRegNCol,
                            labelBy = labelBy,
                            check_sanity = FALSE,
                            isLogged = isLogged)
        })
      })
      session$sendCustomMessage("close_dropDownDeReg", "")
    }
  })

  # Heatmap
  output$deHMSplitColUI <- renderUI({
    otherAvail <- input$deHMcolData
    selectInput("deHMSplitCol", "Split columns by", multiple = TRUE,
                choices = c('condition', otherAvail),
                selected = 'condition')
  })
  output$deHMSplitRowUI <- renderUI({
    otherAvail <- input$deHMrowData
    selectInput("deHMSplitRow", "Split rows by", multiple = TRUE,
                choices = c('regulation', otherAvail),
                selected = 'regulation')
  })

  observeEvent(input$closeDropDownDeHM, {
    session$sendCustomMessage("close_dropDownDeHM", "")
  })
  
  observeEvent(input$dePlotHM, {
    if(!is.null(input$deResSel) &&
       !input$deResSel == ""){
      output$deHeatmap <- renderPlot({
        isolate({
          plotDEGHeatmap(inSCE = sce <- vals$counts,
                         useResult = input$deResSel,
                         doLog = input$deHMDoLog,
                         onlyPos = input$deHMPosOnly,
                         log2fcThreshold = input$deHMFC,
                         fdrThreshold = input$deHMFDR,
                         rowDataName = input$deHMrowData,
                         colDataName = input$deHMcolData,
                         colSplitBy = input$deHMSplitCol,
                         rowSplitBy = input$deHMSplitRow)
        })
      })
      session$sendCustomMessage("close_dropDownDeHM", "")
    }
  })

  #-----------------------------------------------------------------------------
  # Page 5.2: Find Marker ####
  #-----------------------------------------------------------------------------
  observeEvent(input$fmMethod, {
    if (!is.null(vals$counts)) {
      if (is.null(input$fmMethod)) {
        updateSelectInputTag(session, "fmAssay", recommended = c("transformed", "normalized"))
      } else if (input$fmMethod == "DESeq2") {
        updateSelectInputTag(session, "fmAssay", recommended = c("raw"))
      } else {
        updateSelectInputTag(session, "fmAssay", recommended = c("transformed", "normalized"))
      }
    }
  })
  # findMarker RUN ####
  observeEvent(input$runFM, {
    if (is.null(vals$counts)){
      shinyalert::shinyalert("Error!", "Upload data first.", type = "error")
    } else {
      withBusyIndicatorServer("runFM", {
        if(is.na(input$fmLogFC)){
          stop("Log2FC must be a numeric non-empty value!")
        }
        if(is.na(input$fmFDR)){
          stop("FDR must be a numeric non-empty value!")
        }
        vals$counts <- findMarkerDiffExp(inSCE = vals$counts,
                                         method = input$fmMethod,
                                         useAssay = input$fmAssay,
                                         cluster = input$fmCluster,
                                         covariates = input$fmCovar,
                                         log2fcThreshold = input$fmLogFC,
                                         fdrThreshold = input$fmFDR,
                                         minClustExprPerc = input$fmMinClustExprPerc,
                                         maxCtrlExprPerc = input$fmMaxCtrlExprPerc,
                                         minMeanExpr = input$fmMinMeanExpr)
        shinyalert::shinyalert("Success", "Find Marker completed.",
                               "success")
        updateFMPlot()
      })
    }
  })
  # findMarker ResultTable ####
  output$fmResTable <- DT::renderDataTable({
    if(!is.null(vals$counts) &&
       'findMarker' %in% names(metadata(vals$counts))){
      fullTable <- metadata(vals$counts)$findMarker
      fullTable[,5] <- as.factor(fullTable[,5])
      fullTable
    }
  }, filter = "top", options = list(scrollX = TRUE))

  observe({
    if (!is.null(vals$counts) &&
        !is.null(metadata(vals$counts)$findMarker)) {
      shinyjs::enable("fmDownload")
    } else {
      shinyjs::disable("fmDownload")
    }
  })

  output$fmDownload <- downloadHandler(
    filename = function() {
      paste0("findMarkerResult_", input$fmCluster, ".csv")
    },
    content = function(file) {
      fullTable <- metadata(vals$counts)$findMarker
      filteredTable <- fullTable[input$fmResTable_rows_all,]
      utils::write.csv(filteredTable, file, row.names = FALSE)
    }
  )

  # findMarker Heatmap ####
  observeEvent(input$fmShowHMSetting, {
    if (isTRUE(vals$fmHMshowHide)) {
      shinyjs::hide("fmHMsettings")
      updateActionButton(session, "fmShowHMSetting", label = "Show Settings")
      vals$fmHMshowHide <- FALSE
    } else {
      shinyjs::show("fmHMsettings")
      updateActionButton(session, "fmShowHMSetting", label = "Hide Settings")
      vals$fmHMshowHide <- TRUE
    }
  })


  observeEvent(input$fmUseTopN, {
    if (!isTRUE(input$fmUseTopN)) {
      shinyjs::disable("fmTopN")
    } else {
      shinyjs::enable("fmTopN")
    }
  })
  
  observeEvent(input$closeDropDownFM, {
    session$sendCustomMessage("close_dropDownFM", "")
  })

  observeEvent(input$plotFM, {
    updateFMPlot()
    session$sendCustomMessage("close_dropDownFM", "")
  })

  updateFMPlot <- function() {
    if(!is.null(vals$counts) &&
       'findMarker' %in% names(metadata(vals$counts))){
      withBusyIndicatorServer("plotFM", {
        withProgress(message = "Updating marker heatmap...", max = 1, value = 1, {
          if(isTRUE(input$fmUseTopN)
             && is.na(input$fmTopN)){
            stop("Top N marker must be a numeric non-empty value")
          }
          if(is.na(input$fmHMFC)){
            stop("Log2FC must be a numeric non-empty value!")
          }
          if(is.na(input$fmHMFDR)){
            stop("FDR must be a numeric non-empty value!")
          }
          if(!isTRUE(input$fmUseTopN)) {
            topN <- NULL
          } else {
            topN <- input$fmTopN
          }
          # Take value before rendering plot, so that the plot doesn't auto
          # re-render while we tweak the parameter
          output$fmHeatmap <- renderPlot({
            isolate({
              plotMarkerDiffExp(inSCE = vals$counts,
                                orderBy = input$fmHMOrder,
                                log2fcThreshold = input$fmHMFC,
                                topN = topN,
                                fdrThreshold = input$fmHMFDR,
                                decreasing = input$fmHMdec,
                                rowDataName = input$fmHMrowData,
                                colDataName = input$fmHMcolData,
                                minClustExprPerc = input$fmHMMinClustExprPerc,
                                maxCtrlExprPerc = input$fmHMMaxCtrlExprPerc,
                                minMeanExpr = input$fmHMMinMeanExpr,
                                rowLabel = TRUE)
            })
          })
        })
      })
    }
  }


  #-----------------------------------------------------------------------------
  # Page 6: Pathway Activity Analysis
  #-----------------------------------------------------------------------------

  #colData for grouping the data (optional for user)
  observeEvent(input$pathway, {
    if(!is.null(vals$counts)){
      updateSelectInput(session, "pathwayPlotVar", choices = colnames(colData(vals$counts)))
    }
  })

  #select geneset collection name for pathway analysis
  output$selectPathwayGeneLists <- renderUI({
    if (!is.null(vals$counts)){
      if (!is.null(metadata(vals$counts)$sctk$genesets)) {
        #newGSchoices <- sctkListGeneSetCollections(vals$original)
        newGSchoices <- sctkListGeneSetCollections(vals$counts)
        selectizeInput("PathwayGeneLists", "Select Geneset Collection(s):",
                       choices = newGSchoices, multiple = FALSE)

      }
    } else {
      HTML("<h5><span style='color:red'>Must upload data first!</span></h5></br>")
    }
  })

  #Run algorithm
  observeEvent(input$pathwayRun, {
    if (is.null(vals$counts)){
      shinyalert::shinyalert("Error!", "Upload data first.", type = "error")
    } else {
      withBusyIndicatorServer("pathwayRun", {
        #checks

        if(is.null(input$PathwayGeneLists)){
          stop("Must select atleast one Gene List! Gene Lists can be uploaded/selected from 'Import Gene Sets tab'")
        }
        #update metadata of vals$counts
        #metadata(vals$counts)$sctk <- metadata(vals$original)$sctk

        if(input$pathway == "VAM"){

          vals$vamRes <- runVAM(inSCE = vals$counts,
                                useAssay = input$vamAssay,
                                geneSetCollectionName = input$PathwayGeneLists,
                                center = input$vamCenterParameter,
                                gamma = input$vamGammaParameter)


          #vals$vamCdf <- SingleCellExperiment::reducedDim(vals$vamRes, paste0(paste0("VAM_", input$PathwayGeneLists, "_"), "CDF"))
          vals$vamResults <- SingleCellExperiment::reducedDim(vals$vamRes)
          vals$vamScore <- paste0(paste0("VAM_", input$PathwayGeneLists, "_"), "CDF")
          vals$dimreduced <- reducedDims(vals$vamRes)


        }
        else if (input$pathway == "GSVA"){

          vals$gsvaRes <- runGSVA(inSCE = vals$counts,
                                  useAssay = input$vamAssay,
                                  geneSetCollectionName = input$PathwayGeneLists)

          vals$gsvaResults <- SingleCellExperiment::reducedDim(vals$gsvaRes)
          vals$gsvaScore <- paste0(paste0("GSVA_", input$PathwayGeneLists, "_"), "Scores")
          vals$dimreduced <- append(vals$dimreduced, reducedDims(vals$gsvaRes))
        }

      })
    }
  })

  #select geneset for plotting
  output$selectGeneSets <- renderUI({
    if(input$pathway == "VAM"){
      if (!(is.null(vals$vamRes))){
        selectizeInput("GeneSets", "Select Geneset:",
                       choices = colnames(reducedDim(vals$vamRes, vals$vamScore)), multiple = FALSE)
      }
    }
    else if(input$pathway == "GSVA"){

      if (!(is.null(vals$gsvaRes))){
        selectizeInput("GeneSets", "Select Geneset:",
                       choices = colnames(reducedDim(vals$gsvaRes, vals$gsvaScore)), multiple = FALSE)
      }
    }
  })

  #select from reducedDimNames
  output$selectReduceDim <- renderUI({
    if (!(is.null(vals$counts))){
      selectizeInput("reducedDimNames", "Select Score matrix which you want to plot:",
                     choices = names((vals$dimreduced)), multiple = FALSE)
    }
  })


  #plot results with default values intitially
  observeEvent(input$pathwayRun, {
    output$pathwayPlot <- renderPlot({
      if (input$pathway == "VAM"){
        if (!(is.null(vals$vamRes))){
          plotSCEViolin(inSCE = vals$vamRes, slotName = "reducedDims", itemName = vals$vamScore, dimension = colnames(reducedDim(vals$vamRes, vals$vamScore))[[1]], xlab = "sample", ylab = colnames(reducedDim(vals$vamRes, vals$vamScore))[[1]])
          #groupby can be used for indicating colnames (optional)
        }

      }
      else if (input$pathway == "GSVA"){
        if (!(is.null(vals$gsvaRes))){
          plotSCEViolin(inSCE = vals$gsvaRes, slotName = "reducedDims", itemName = vals$gsvaScore, dimension = colnames(reducedDim(vals$gsvaRes, vals$gsvaScore))[[1]], xlab = "Sample", ylab = colnames(reducedDim(vals$gsvaRes, vals$gsvaScore))[[1]])
          #groupby can be used for indicating colnames (optional)

        }

      }


    })
  })




 #plot results
  observeEvent(input$Plot, {
    output$pathwayPlot <- renderPlot({
      isolate({
      if (input$pathway == "VAM"){
        if (!(is.null(vals$vamRes))){
          plotSCEViolin(inSCE = vals$vamRes, slotName = "reducedDims", itemName = input$reducedDimNames, dimension = input$GeneSets, xlab = "sample", ylab = input$GeneSets, groupBy = input$pathwayPlotVar, violin = input$violinplot, boxplot = input$boxplot, summary = input$summary)
          #groupby can be used for indicating colnames (optional)
        }

      }
      else if (input$pathway == "GSVA"){
        if (!(is.null(vals$gsvaRes))){
          plotSCEViolin(inSCE = vals$gsvaRes, slotName = "reducedDims", itemName = input$reducedDimNames, dimension = input$GeneSets, xlab = "Sample", ylab = input$GeneSets, groupBy = input$pathwayPlotVar, violin = input$violinplot, boxplot = input$boxplot, summary = input$summary)
          #groupby can be used for indicating colnames (optional)

        }

      }

    })
     })
    session$sendCustomMessage("close_dropDownPathway", "")
  })

  observeEvent(input$closeDropDownPathway,{
    session$sendCustomMessage("close_dropDownPathway", "")
  })

  #disable downloadPathway button if the pathway data doesn't exist
  isVamResult <- reactive(is.null(vals$vamResults))
  isGsvaResult <- reactive(is.null(vals$gsvaResults))
  observe({
    if (isVamResult() && isGsvaResult()) {
      shinyjs::disable("downloadPathway")
    } else {
      shinyjs::enable("downloadPathway")
    }
  })

  #download pathway results
  output$downloadPathway <- downloadHandler(
    filename = function() {
      paste("Pathway_results-", Sys.Date(), ".csv", sep = "")
    },
    content = function(file) {
      if(input$pathway == "VAM"){
        utils::write.csv(vals$vamResults, file)
      }
      else if (input$pathway == "GSVA"){
        utils::write.csv(vals$gsvaResults, file)
      }
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
    isoDbs <- isolate(dbs())
    nTabs <- length(isoDbs)
    #create tabPanel with datatable in it
    myTabs <- lapply(seq_len((nTabs)), function(i) {
      tabPanel(paste0(isoDbs[i]),
               DT::dataTableOutput(paste0(isoDbs[i]))
      )
    })
    do.call(tabsetPanel, myTabs)
  })

  #create datatables
  observe({
    req(vals$enrichRes)
    isoDbs <- isolate(dbs())
    enrResults <- vals$enrichRes[, c(1:10)] %>%
      mutate(Database_selected =
               paste0("<a href='", vals$enrichRes[, 11],
                      "' target='_blank'>",
                      vals$enrichRes[, 1], "</a>"))
    lapply(seq_len(length(isoDbs)), function(i){
      output[[paste0(isoDbs[i])]] <- DT::renderDataTable({
        DT::datatable({
          enr <- enrResults[which(vals$enrichRes[, 1] %in% isoDbs[i]), ]
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
        if(is.na(input$minCount)){
          stop("Minimum readcount must be a non-empty numeric value!")
        }
        if(is.na(input$minCells)){
          stop("Minimum number of cells must be a non-empty numeric value!")
        }
        if(is.na(input$iterations)){
          stop("Number of bootstrap iterations must be a non-empty numeric value!")
        }
        vals$subDepth <- downSampleDepth(originalData = vals$counts,
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
        if(is.na(input$minCellNum)
           || is.na(input$maxCellNum)
           || is.na(input$iterations)
           || is.na(input$totalReads)
           || is.na(input$minCount)
           || is.na(input$minCells)
           || is.na(input$depthResolution)){
          stop("One or more parameter values are empty!")
        }
        if (input$useReadCount){
          vals$subCells <- downSampleCells(originalData = vals$counts,
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
          vals$subCells <- downSampleCells(originalData = vals$counts,
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
        if(is.na(input$numCellsSnap)
           || is.na(input$numReadsSnap)
           || is.na(input$iterationsSnap)){
          stop("One or more parameter values are empty!")
        }
        vals$snapshot <- iterateSimulations(originalData = vals$counts,
                                            useAssay = input$snapshotAssay,
                                            realLabels = input$selectSnapshotCondition,
                                            totalReads = input$numReadsSnap,
                                            cells = input$numCellsSnap,
                                            iterations = input$iterationsSnap)
        vals$effectSizes <- calcEffectSizes(countMatrix = expData(vals$counts, input$snapshotAssay), condition = colData(vals$counts)[, input$selectSnapshotCondition])
        output$Snaplot <- renderPlot({
          plot(apply(vals$snapshot, 1, function(x){sum(x <= 0.05) / length(x)}) ~ vals$effectSizes,
               xlab = "Cohen's d effect size", ylab = "Detection power", lwd = 4, main = "Power to detect diffex by effect size")
        })
      })
    }
  })

  #-----------------------------------------------------------------------------
  # Page 8: Seurat Workflow
  #-----------------------------------------------------------------------------

  #Perform normalization
  observeEvent(input$normalize_button, {
    req(vals$counts)
    withProgress(message = "Normalizing", max = 1, value = 1, {
      vals$counts <- seuratNormalizeData(inSCE = vals$counts,
                                         useAssay = input$seuratSelectNormalizationAssay,
                                         normAssayName = "seuratNormData",
                                         normalizationMethod = input$normalization_method,
                                         scaleFactor = as.numeric(input$scale_factor))
      metadata(vals$counts)$sctk$seuratUseAssay <- input$seuratSelectNormalizationAssay
      # updateAssayInputs()
      vals$counts <- singleCellTK:::.seuratInvalidate(inSCE = vals$counts)
    })
    updateCollapse(session = session, "SeuratUI", style = list("Normalize Data" = "success"))
    shinyjs::enable(selector = "div[value='Scale Data']")
    S4Vectors::metadata(vals$counts)$seuratMarkers <- NULL
    shinyjs::hide(
      selector = "div[value='Downstream Analysis']")
    showNotification("Normalization Complete")
  })

  #Perform scaling
  observeEvent(input$scale_button, {
    req(vals$counts)
    withProgress(message = "Scaling", max = 1, value = 1, {
      vals$counts <- seuratScaleData(inSCE = vals$counts,
                                     useAssay = "seuratNormData",
                                     scaledAssayName = "seuratScaledData",
                                     #model = input$model.use,
                                     scale = input$do.scale,
                                     center = input$do.center,
                                     scaleMax = input$scale.max)
      # updateAssayInputs()
      vals$counts <- singleCellTK:::.seuratInvalidate(inSCE = vals$counts, scaleData = FALSE)
    })
    updateCollapse(session = session, "SeuratUI", style = list("Scale Data" = "success"))
    shinyjs::enable(selector = "div[value='Highly Variable Genes']")
    S4Vectors::metadata(vals$counts)$seuratMarkers <- NULL
    shinyjs::hide(
      selector = "div[value='Downstream Analysis']")
    showNotification("Scale Complete")
  })

  #Find HVG
  observeEvent(input$find_hvg_button, {
    req(vals$counts)
    withProgress(message = "Finding highly variable genes", max = 1, value = 1, {
      if(input$hvg_method == "vst"){
        vals$counts <- seuratFindHVG(inSCE = vals$counts,
                                     useAssay = metadata(vals$counts)$sctk$seuratUseAssay,
                                     hvgMethod = input$hvg_method,
                                     hvgNumber = as.numeric(input$hvg_no_features))
      }
      else{
        vals$counts <- seuratFindHVG(inSCE = vals$counts,
                                     useAssay = "seuratNormData",
                                     hvgMethod = input$hvg_method,
                                     hvgNumber = as.numeric(input$hvg_no_features))
      }
      vals$counts <- singleCellTK:::.seuratInvalidate(inSCE = vals$counts, scaleData = FALSE, varFeatures = FALSE)
    })
    withProgress(message = "Plotting HVG", max = 1, value = 1, {
      output$plot_hvg <- renderPlotly({
        isolate({
          plotly::ggplotly(seuratPlotHVG(vals$counts, input$hvg_no_features_view))
        })
      })
    })
    updateCollapse(session = session, "SeuratUI", style = list("Highly Variable Genes" = "success"))
    shinyjs::enable(selector = "div[value='Dimensionality Reduction']")
    S4Vectors::metadata(vals$counts)$seuratMarkers <- NULL
    shinyjs::hide(
      selector = "div[value='Downstream Analysis']")
    showNotification("Find HVG Complete")
  })

  #Display highly variable genes
  output$hvg_output <- renderText({
    if (!is.null(vals$counts)) {
      if (!is.null(vals$counts@metadata$seurat$obj)) {
        if (length(slot(vals$counts@metadata$seurat$obj, "assays")[["RNA"]]@var.features) > 0) {
          isolate({
            singleCellTK:::.seuratGetVariableFeatures(vals$counts, input$hvg_no_features_view)
          })
        }
      }
    }
  })

  #Run PCA
  observeEvent(input$run_pca_button, {
    req(vals$counts)

    #remove tabs if not generated
    removeTab(inputId = "seuratPCAPlotTabset", target = "PCA Plot")
    removeTab(inputId = "seuratPCAPlotTabset", target = "Elbow Plot")
    removeTab(inputId = "seuratPCAPlotTabset", target = "JackStraw Plot")
    removeTab(inputId = "seuratPCAPlotTabset", target = "Heatmap Plot")

    withProgress(message = "Running PCA", max = 1, value = 1, {
      vals$counts <- seuratPCA(inSCE = vals$counts,
                               useAssay = "seuratScaledData",
                               reducedDimName = "seuratPCA",
                               nPCs = input$pca_no_components)

      vals$counts@metadata$seurat$count_pc <- dim(convertSCEToSeurat(vals$counts)[["pca"]])[2]
      vals$counts <- singleCellTK:::.seuratInvalidate(inSCE = vals$counts, scaleData = FALSE, varFeatures = FALSE, PCA = FALSE, ICA = FALSE)
    })

    appendTab(inputId = "seuratPCAPlotTabset", tabPanel(title = "PCA Plot",
                                                        panel(heading = "PCA Plot",
                                                              plotlyOutput(outputId = "plot_pca")
                                                        )
    ), select = TRUE)

    withProgress(message = "Plotting PCA", max = 1, value = 1, {
      output$plot_pca <- renderPlotly({
        isolate({
          plotly::ggplotly(seuratReductionPlot(inSCE = vals$counts,
                                               useReduction = "pca",
                                               showLegend = FALSE))
        })
      })
    })
    if (input$pca_compute_elbow) {
      appendTab(inputId = "seuratPCAPlotTabset", tabPanel(title = "Elbow Plot",
                                                          panel(
                                                            heading = "Elbow Plot",
                                                                plotlyOutput(outputId = "plot_elbow_pca")
                                                          )
      ))

      withProgress(message = "Generating Elbow Plot", max = 1, value = 1, {
        updateNumericInput(session = session, inputId = "pca_significant_pc_counter", value = singleCellTK:::.computeSignificantPC(vals$counts))
        output$plot_elbow_pca <- renderPlotly({
          isolate({
            seuratElbowPlot(inSCE = vals$counts,
                            significantPC = singleCellTK:::.computeSignificantPC(vals$counts))
          })
        })
        output$pca_significant_pc_output <- renderText({
          isolate({
            paste("<p>Number of significant components suggested by ElbowPlot: <span style='color:red'>", singleCellTK:::.computeSignificantPC(vals$counts)," </span> </p> <hr>")
          })
        })
      })
    }
    if (input$pca_compute_jackstraw) {
      appendTab(inputId = "seuratPCAPlotTabset", tabPanel(title = "JackStraw Plot",
                                                          panel(heading = "JackStraw Plot",
                                                                plotlyOutput(outputId = "plot_jackstraw_pca")
                                                          )
      ))

      withProgress(message = "Generating JackStraw Plot", max = 1, value = 1, {
        vals$counts <- seuratComputeJackStraw(inSCE = vals$counts,
                                              useAssay = "seuratScaledData",
                                              dims = input$pca_no_components)
        output$plot_jackstraw_pca <- renderPlotly({
          isolate({
            plotly::ggplotly(seuratJackStrawPlot(inSCE = vals$counts,
                                                 dims = input$pca_no_components))
          })
        })
      })
    }
    if (input$pca_compute_heatmap) {
      appendTab(inputId = "seuratPCAPlotTabset", tabPanel(title = "Heatmap Plot",
                                                          panel(heading = "Heatmap Plot",
                                                                panel(heading = "Plot Options",
                                                                      fluidRow(
                                                                        column(4, dropdown(
                                                                          fluidRow(
                                                                            column(12,
                                                                                   fluidRow(actionBttn(inputId = "closeDropDownSeuratHM", label = NULL, style = "simple", color = "danger", icon = icon("times"), size = "xs"), align = "right"),
                                                                                   fluidRow(
                                                                                     column(6,
                                                                                            pickerInput(inputId = "picker_dimheatmap_components_pca", label = "Select principal components to plot:", choices = c(), options = list(`actions-box` = TRUE, size = 10, `selected-text-format` = "count > 3"), multiple = TRUE)
                                                                                     ),
                                                                                     column(6,
                                                                                            sliderInput(inputId = "slider_dimheatmap_pca", label = "Number of columns for the plot: ", min = 1, max = 4, value = 2)
                                                                                     )
                                                                                   ),
                                                                                   actionBttn(
                                                                                     inputId = "plot_heatmap_pca_button",
                                                                                     label = "Update",
                                                                                     style = "bordered",
                                                                                     color = "primary",
                                                                                     size = "sm"
                                                                                   )
                                                                            )
                                                                          ),
                                                                          inputId = "dropDownSeuratHM",
                                                                          icon = icon("cog"),
                                                                          status = "primary",
                                                                          circle = FALSE,
                                                                          inline = TRUE
                                                                        )),
                                                                        column(7, fluidRow(h6("Heatmaps of the top features correlated with each component"), align="center"))
                                                                      )
                                                                ),
                                                                panel(heading = "Plot",
                                                                      shinyjqui::jqui_resizable(plotOutput(outputId = "plot_heatmap_pca"), options = list(maxWidth = 700))
                                                                )
                                                          )
      ))

      withProgress(message = "Generating Heatmaps", max = 1, value = 1, {
        vals$counts@metadata$seurat$heatmap_pca <- seuratComputeHeatmap(inSCE = vals$counts,
                                                                        useAssay = "seuratScaledData",
                                                                        useReduction = "pca",
                                                                        dims = input$pca_no_components,
                                                                        nfeatures = input$pca_compute_heatmap_nfeatures,
                                                                        combine = FALSE,
                                                                        fast = FALSE)
        output$plot_heatmap_pca <- renderPlot({
          isolate({
            seuratHeatmapPlot(plotObject = vals$counts@metadata$seurat$heatmap_pca,
                              dims = input$pca_no_components,
                              ncol = 2,
                              labels = c("PC1", "PC2", "PC3", "PC4"))
          })
        })
        updatePickerInput(session = session, inputId = "picker_dimheatmap_components_pca", choices = singleCellTK:::.getComponentNames(vals$counts@metadata$seurat$count_pc, "PC"))
      })
    }
    updateCollapse(session = session, "SeuratUI", style = list("Dimensionality Reduction" = "success"))

    #Enable/Disable PCA plot panels not selected for computation (ElbowPlot, JackStraw or Heatmap)
    shinyjs::enable(
      selector = ".seurat_pca_plots a[data-value='PCA Plot']")

    shinyjs::toggleState(
      selector = ".seurat_pca_plots a[data-value='Elbow Plot']",
      condition = input$pca_compute_elbow)

    shinyjs::toggleState(
      selector = ".seurat_pca_plots a[data-value='JackStraw Plot']",
      condition = input$pca_compute_jackstraw)

    shinyjs::toggleState(
      selector = ".seurat_pca_plots a[data-value='Heatmap Plot']",
      condition = input$pca_compute_heatmap)

    shinyjs::enable(
      selector = "div[value='tSNE/UMAP']")

    shinyjs::show(selector = ".seurat_pca_plots")

    S4Vectors::metadata(vals$counts)$seuratMarkers <- NULL
    shinyjs::hide(
      selector = "div[value='Downstream Analysis']")

    showNotification("PCA Complete")
  })
  
  observeEvent(input$closeDropDownSeuratHM,{
    session$sendCustomMessage("close_dropDownSeuratHM", "")
  })

  #Run ICA
  observeEvent(input$run_ica_button, {
    req(vals$counts)

    #remove tabs if not generated
    removeTab(inputId = "seuratICAPlotTabset", target = "ICA Plot")
    removeTab(inputId = "seuratICAPlotTabset", target = "Heatmap Plot")

    withProgress(message = "Running ICA", max = 1, value = 1, {
      vals$counts <- seuratICA(inSCE = vals$counts,
                               useAssay = "seuratScaledData",
                               nics = input$ica_no_components)

      vals$counts@metadata$seurat$count_ic <- dim(convertSCEToSeurat(vals$counts)[["ica"]])[2]
      vals$counts <- singleCellTK:::.seuratInvalidate(inSCE = vals$counts, scaleData = FALSE, varFeatures = FALSE, PCA = FALSE, ICA = FALSE)
    })

    appendTab(inputId = "seuratICAPlotTabset", tabPanel(title = "ICA Plot",
                                                        panel(heading = "ICA Plot",
                                                              plotlyOutput(outputId = "plot_ica")
                                                        )
    ), select = TRUE)

    withProgress(message = "Plotting ICA", max = 1, value = 1, {
      output$plot_ica <- renderPlotly({
        isolate({
          plotly::ggplotly(seuratReductionPlot(inSCE = vals$counts,
                                               useReduction = "ica",
                                               showLegend = FALSE))
        })
      })
    })
    if (input$ica_compute_heatmap) {
      appendTab(inputId = "seuratICAPlotTabset", tabPanel(title = "Heatmap Plot",
                                                          panel(heading = "Heatmap Plot",
                                                                panel(heading = "Plot Options",
                                                                      fluidRow(
                                                                        column(6,
                                                                               pickerInput(inputId = "picker_dimheatmap_components_ica", label = "Select principal components to plot:", choices = c(), options = list(`actions-box` = TRUE, size = 10, `selected-text-format` = "count > 3"), multiple = TRUE)
                                                                        ),
                                                                        column(6,
                                                                               sliderInput(inputId = "slider_dimheatmap_ica", label = "Number of columns for the plot: ", min = 1, max = 4, value = 2)
                                                                        )
                                                                      ),
                                                                      actionButton(inputId = "plot_heatmap_ica_button", "Plot")
                                                                ),
                                                                panel(heading = "Plot",
                                                                      shinyjqui::jqui_resizable(plotOutput(outputId = "plot_heatmap_ica"), options = list(maxWidth = 700))
                                                                )
                                                          )
      ))

      withProgress(message = "Generating Heatmaps", max = 1, value = 1, {
        vals$counts@metadata$seurat$heatmap_ica <- seuratComputeHeatmap(inSCE = vals$counts,
                                                                        useAssay = "seuratScaledData",
                                                                        useReduction = "ica",
                                                                        dims = input$ica_no_components,
                                                                        nfeatures = input$ica_compute_heatmap_nfeatures,
                                                                        combine = FALSE,
                                                                        fast = FALSE)
        output$plot_heatmap_ica <- renderPlot({
          isolate({
            seuratHeatmapPlot(plotObject = vals$counts@metadata$seurat$heatmap_ica,
                              dims = input$ica_no_components,
                              ncol = 2,
                              labels = c("IC1", "IC2", "IC3", "IC4"))
          })
        })
        updatePickerInput(session = session, inputId = "picker_dimheatmap_components_ica", choices = singleCellTK:::.getComponentNames(vals$counts@metadata$seurat$count_ic, "IC"))
      })
    }
    updateCollapse(session = session, "SeuratUI", style = list("Dimensionality Reduction" = "success"))

    #Enable/Disable ICA plot panels not selected for computation (Heatmap)
    shinyjs::enable(
      selector = ".seurat_ica_plots a[data-value='ICA Plot']")

    shinyjs::toggleState(
      selector = ".seurat_ica_plots a[data-value='Heatmap Plot']",
      condition = input$ica_compute_heatmap)

    shinyjs::enable(
      selector = "div[value='tSNE/UMAP']")

    shinyjs::show(selector = ".seurat_ica_plots")

    S4Vectors::metadata(vals$counts)$seuratMarkers <- NULL
    shinyjs::hide(
      selector = "div[value='Downstream Analysis']")

    showNotification("ICA Complete")
  })

  #Find clusters
  observeEvent(input$find_clusters_button, {
    req(vals$counts)
    if(!is.null(slot(vals$counts@metadata$seurat$obj, "reductions")[[input$reduction_clustering_method]])){

      #Remove plot tabs if generated before
      removeTab(inputId = "seuratClusteringPlotTabset", target = "PCA Plot")
      removeTab(inputId = "seuratClusteringPlotTabset", target = "ICA Plot")
      removeTab(inputId = "seuratClusteringPlotTabset", target = "tSNE Plot")
      removeTab(inputId = "seuratClusteringPlotTabset", target = "UMAP Plot")


      withProgress(message = "Finding clusters", max = 1, value = 1, {
        vals$counts <- seuratFindClusters(inSCE = vals$counts,
                                          useAssay = "seuratScaledData",
                                          useReduction = input$reduction_clustering_method,
                                          dims = input$pca_significant_pc_counter,
                                          algorithm = input$algorithm.use,
                                          groupSingletons = input$group.singletons,
                                          resolution = input$resolution_clustering)
      })
      updateCollapse(session = session, "SeuratUI", style = list("Clustering" = "success"))
      showNotification("Find Clusters Complete")

      if(!is.null(slot(vals$counts@metadata$seurat$obj, "reductions")[["pca"]])){
        appendTab(inputId = "seuratClusteringPlotTabset", tabPanel(title = "PCA Plot",
                                                                   panel(heading = "PCA Plot",
                                                                         plotlyOutput(outputId = "plot_pca_clustering")
                                                                   )
        ), select = TRUE

        )
        withProgress(message = "Re-generating PCA plot with cluster labels", max = 1, value = 1,{
          output$plot_pca_clustering <- renderPlotly({
            isolate({
              plotly::ggplotly(seuratReductionPlot(inSCE = vals$counts,
                                                   useReduction = "pca",
                                                   showLegend = TRUE))
            })
          })
        })
        shinyjs::toggleState(
          selector = ".seurat_clustering_plots a[data-value='PCA Plot']",
          condition = !is.null(slot(vals$counts@metadata$seurat$obj, "reductions")[["pca"]]))
      }
      if(!is.null(slot(vals$counts@metadata$seurat$obj, "reductions")[["ica"]])){
        appendTab(inputId = "seuratClusteringPlotTabset", tabPanel(title = "ICA Plot",
                                                                   panel(heading = "ICA Plot",
                                                                         plotlyOutput(outputId = "plot_ica_clustering")
                                                                   )
        ), select = TRUE)
        withProgress(message = "Re-generating ICA plot with cluster labels", max = 1, value = 1,{
          output$plot_ica_clustering <- renderPlotly({
            isolate({
              plotly::ggplotly(seuratReductionPlot(inSCE = vals$counts,
                                                   useReduction = "ica",
                                                   showLegend = TRUE))
            })
          })
        })
        shinyjs::toggleState(
          selector = ".seurat_clustering_plots a[data-value='ICA Plot']",
          condition = !is.null(slot(vals$counts@metadata$seurat$obj, "reductions")[["ica"]]))
      }

      if(!is.null(slot(vals$counts@metadata$seurat$obj, "reductions")[["tsne"]])){
        appendTab(inputId = "seuratClusteringPlotTabset", tabPanel(title = "tSNE Plot",
                                                                   panel(heading = "tSNE Plot",
                                                                         plotlyOutput(outputId = "plot_tsne_clustering")
                                                                   )
        )
        )

        withProgress(message = "Re-generating tSNE plot with cluster labels", max = 1, value = 1,{
          output$plot_tsne_clustering <- renderPlotly({
            isolate({
              plotly::ggplotly(seuratReductionPlot(inSCE = vals$counts,
                                                   useReduction = "tsne",
                                                   showLegend = TRUE))
            })
          })
        })
        shinyjs::toggleState(
          selector = ".seurat_clustering_plots a[data-value='tSNE Plot']",
          condition = !is.null(slot(vals$counts@metadata$seurat$obj, "reductions")[["tsne"]]))
      }

      if(!is.null(slot(vals$counts@metadata$seurat$obj, "reductions")[["umap"]])){
        appendTab(inputId = "seuratClusteringPlotTabset", tabPanel(title = "UMAP Plot",
                                                                   panel(heading = "UMAP Plot",
                                                                         plotlyOutput(outputId = "plot_umap_clustering")
                                                                   )
        )
        )
        withProgress(message = "Re-generating UMAP plot with cluster labels", max = 1, value = 1,{
          output$plot_umap_clustering <- renderPlotly({
            isolate({
              plotly::ggplotly(seuratReductionPlot(inSCE = vals$counts,
                                                   useReduction = "umap",
                                                   showLegend = TRUE))
            })
          })
        })
        shinyjs::toggleState(
          selector = ".seurat_clustering_plots a[data-value='UMAP Plot']",
          condition = !is.null(slot(vals$counts@metadata$seurat$obj, "reductions")[["umap"]]))
      }

      shinyjs::show(selector = ".seurat_clustering_plots")

      #enable find marker selection
      shinyjs::enable(
        selector = "div[value='Find Markers']")

      #update colData names
      updateColDataNames()

      S4Vectors::metadata(vals$counts)$seuratMarkers <- NULL
      shinyjs::hide(
        selector = "div[value='Downstream Analysis']")

      #populate updated colData items for findMarkers tab
      updateSelectInput(session = session,
                        inputId = "seuratFindMarkerSelectPhenotype",
                        choices = colnames(colData(vals$counts)))

      #populate reducDim objects from seuratObject for findMarkers tab
      updateSelectInput(session = session,
                        inputId = "seuratFindMarkerReductionMethod",
                        choices = Seurat::Reductions(convertSCEToSeurat(vals$counts)))

    }
    else{
      showNotification(paste0("'", input$reduction_clustering_method, "' reduction not found in input object"))
    }
  })

  observeEvent(input$seuratFindMarkerSelectPhenotype,{
    if(!is.null(vals$counts)){
      updateSelectInput(
        session = session,
        inputId = "seuratFindMarkerGroup1",
        choices = unique(colData(vals$counts)[[input$seuratFindMarkerSelectPhenotype]])
      )
      updateSelectInput(
        session = session,
        inputId = "seuratFindMarkerGroup2",
        choices = unique(colData(vals$counts)[[input$seuratFindMarkerSelectPhenotype]])
      )
    }
  })

  observeEvent(input$seuratFindMarkerGroup1,{
    if(!is.null(vals$counts)){
      matchedIndex <- match(input$seuratFindMarkerGroup1,  unique(colData(vals$counts)[[input$seuratFindMarkerSelectPhenotype]]))
      if(!is.na(matchedIndex)){
        updateSelectInput(
          session = session,
          inputId = "seuratFindMarkerGroup2",
          choices = unique(colData(vals$counts)[[input$seuratFindMarkerSelectPhenotype]])[-matchedIndex]
        )
      }
    }
  })



  observeEvent(input$seuratFindMarkerRun,{
    withProgress(message = "Finding markers", max = 1, value = 1,{
      if(input$seuratFindMarkerType == "markerAll"){
        vals$counts <- seuratFindMarkers(inSCE = vals$counts,
                                         allGroup = input$seuratFindMarkerSelectPhenotype,
                                         test = input$seuratFindMarkerTest,
                                         onlyPos = input$seuratFindMarkerPosOnly)
      }
      else{
        indices1 <- which(colData(vals$counts)[[input$seuratFindMarkerSelectPhenotype]] == input$seuratFindMarkerGroup1, arr.ind = TRUE)
        indices2 <- which(colData(vals$counts)[[input$seuratFindMarkerSelectPhenotype]] == input$seuratFindMarkerGroup2, arr.ind = TRUE)
        cells1 <- colnames(vals$counts)[indices1]
        cells2 <- colnames(vals$counts)[indices2]
        if(input$seuratFindMarkerType == "markerConserved"){
          vals$counts <- seuratFindMarkers(inSCE = vals$counts,
                                           cells1 = cells1,
                                           cells2 = cells2,
                                           group1 = input$seuratFindMarkerGroup1,
                                           group2 = input$seuratFindMarkerGroup2,
                                           conserved = TRUE,
                                           test = input$seuratFindMarkerTest,
                                           onlyPos = input$seuratFindMarkerPosOnly)
        }
        else{
          vals$counts <- seuratFindMarkers(inSCE = vals$counts,
                                           cells1 = cells1,
                                           cells2 = cells2,
                                           group1 = input$seuratFindMarkerGroup1,
                                           group2 = input$seuratFindMarkerGroup2,
                                           test = input$seuratFindMarkerTest,
                                           onlyPos = input$seuratFindMarkerPosOnly)
        }
      }
    })


    shinyjs::show(selector = ".seurat_findmarker_table")
    shinyjs::show(selector = ".seurat_findmarker_jointHeatmap")
    shinyjs::show(selector = ".seurat_findmarker_plots")

    removeTab(inputId = "seuratFindMarkerPlotTabset", target = "Ridge Plot")
    removeTab(inputId = "seuratFindMarkerPlotTabset", target = "Violin Plot")
    removeTab(inputId = "seuratFindMarkerPlotTabset", target = "Feature Plot")
    removeTab(inputId = "seuratFindMarkerPlotTabset", target = "Dot Plot")
    removeTab(inputId = "seuratFindMarkerPlotTabset", target = "Heatmap Plot")

    appendTab(inputId = "seuratFindMarkerPlotTabset", tabPanel(title = "Ridge Plot",
                                                               panel(heading = "Ridge Plot",
                                                                     fluidRow(
                                                                       column(12, align = "center",
                                                                              panel(
                                                                                HTML(paste("<span style='color:red'>Select genes from the above table to plot!</span>"))
                                                                              )
                                                                              )
                                                                     )
                                                               )
    )
    )
    appendTab(inputId = "seuratFindMarkerPlotTabset", tabPanel(title = "Violin Plot",
                                                               panel(heading = "Violin Plot",
                                                                     fluidRow(
                                                                       column(12, align = "center",
                                                                              panel(
                                                                                HTML(paste("<span style='color:red'>Select genes from the above table to plot!</span>"))
                                                                              )
                                                                       )
                                                                     )
                                                               )
    )
    )
    appendTab(inputId = "seuratFindMarkerPlotTabset", tabPanel(title = "Feature Plot",
                                                               panel(heading = "Feature Plot",
                                                                     fluidRow(
                                                                       column(12, align = "center",
                                                                              panel(
                                                                                HTML(paste("<span style='color:red'>Select genes from the above table to plot!</span>"))
                                                                              )
                                                                       )
                                                                     )
                                                               )
    )
    )
    appendTab(inputId = "seuratFindMarkerPlotTabset", tabPanel(title = "Dot Plot",
                                                               panel(heading = "Dot Plot",
                                                                     fluidRow(
                                                                       column(12, align = "center",
                                                                              panel(
                                                                                HTML(paste("<span style='color:red'>Select genes from the above table to plot!</span>"))
                                                                              )
                                                                       )
                                                                     )
                                                               )
    )
    )
    appendTab(inputId = "seuratFindMarkerPlotTabset", tabPanel(title = "Heatmap Plot",
                                                               panel(heading = "Heatmap Plot",
                                                                     fluidRow(
                                                                       column(12, align = "center",
                                                                              panel(
                                                                                HTML(paste("<span style='color:red'>Select genes from the above table to plot!</span>"))
                                                                              )
                                                                       )
                                                                     )
                                                               )
    )
    )

    #df <- metadata(vals$counts)$seuratMarkers[which(metadata(vals$counts)$seuratMarkers$p_val_adj < 0.05, arr.ind = TRUE),]
    df <- metadata(vals$counts)$seuratMarkers
    seuratObject <- convertSCEToSeurat(vals$counts, scaledAssay = "seuratScaledData")
    indices <- list()
    cells <- list()
    groups <- unique(colData(vals$counts)[[input$seuratFindMarkerSelectPhenotype]])
    for(i in seq(length(groups))){
      indices[[i]] <- which(colData(vals$counts)[[input$seuratFindMarkerSelectPhenotype]] == groups[i], arr.ind = TRUE)
      cells[[i]] <- colnames(vals$counts)[indices[[i]]]
      cells[[i]] <- lapply(
        X = cells[[i]],
        FUN = function(t) gsub(
          pattern = "_",
          replacement = "-",
          x = t,
          fixed = TRUE)
      )
      Idents(seuratObject, cells = cells[[i]]) <- groups[i]
    }

    showTab(inputId = "seuratFindMarkerPlotTabset", target = "Joint Heatmap Plot")
    updateTabsetPanel(session = session, inputId = "seuratFindMarkerPlotTabset", selected = "Ridge Plot")
    shinyjs::show(selector = ".seurat_findmarker_plots")

     # Output the heatmap
     colnames(df)[which(startsWith(colnames(df), "avg") == TRUE)] <- "avg_log2FC"
     top10markers <- df %>% group_by(cluster1) %>% arrange(desc(avg_log2FC)) %>% slice_head(n=10)
     # Subset seuratObject to contain only cells available in selected clusters
     if(input$seuratFindMarkerType != "markerAll"){
       subsetIdents <- c(unique(top10markers$cluster1), unique(top10markers$cluster2))
       subsetIdents <- subsetIdents[subsetIdents!="all"]
       seuratObject <- subset(seuratObject, idents = subsetIdents) 
     }
     # Plot heatmap
     output$findMarkerHeatmapPlotFull <- renderPlot({
       isolate({
         DoHeatmap(seuratObject, features = top10markers$gene.id)
       })
     })

     # output$findMarkerHeatmapPlotFullTopText <- renderUI({
     #   h6(paste("Heatmap plotted across all groups against genes with adjusted p-values <", input$seuratFindMarkerPValAdjInput))
     # })

    showNotification("Find Markers Complete")

    # Show downstream analysis options
    callModule(module = nonLinearWorkflow, id = "nlw-seurat", parent = session, de = TRUE, pa = TRUE)

    updateCollapse(session = session, "SeuratUI", style = list("Find Markers" = "success"))

    updateCollapse(session = session, "SeuratUI", style = list("Downstream Analysis" = "info"))

    removeTab(inputId = "seuratFindMarkerPlotTabset", target = "Ridge Plot")
    removeTab(inputId = "seuratFindMarkerPlotTabset", target = "Violin Plot")
    removeTab(inputId = "seuratFindMarkerPlotTabset", target = "Feature Plot")
    removeTab(inputId = "seuratFindMarkerPlotTabset", target = "Dot Plot")
    removeTab(inputId = "seuratFindMarkerPlotTabset", target = "Heatmap Plot")

    appendTab(inputId = "seuratFindMarkerPlotTabset", tabPanel(title = "Ridge Plot",
                                                                                             panel(heading = "Ridge Plot",
                                                                                                   shinyjqui::jqui_resizable(
                                                                                                     plotOutput(outputId = "findMarkerRidgePlot")
                                                                                                   )
                                                                                             )
    )
    )
    appendTab(inputId = "seuratFindMarkerPlotTabset", tabPanel(title = "Violin Plot",
                                                                                             panel(heading = "Violin Plot",
                                                                                                   shinyjqui::jqui_resizable(
                                                                                                     plotOutput(outputId = "findMarkerViolinPlot")
                                                                                                   )
                                                                                             )
    )
    )
    appendTab(inputId = "seuratFindMarkerPlotTabset", tabPanel(title = "Feature Plot",
                                                                                             panel(heading = "Feature Plot",
                                                                                                   shinyjqui::jqui_resizable(
                                                                                                     plotOutput(outputId = "findMarkerFeaturePlot")
                                                                                                   )
                                                                                             )
    )
    )
    appendTab(inputId = "seuratFindMarkerPlotTabset", tabPanel(title = "Dot Plot",
                                                                                             panel(heading = "Dot Plot",
                                                                                                   shinyjqui::jqui_resizable(
                                                                                                     plotOutput(outputId = "findMarkerDotPlot")
                                                                                                   )
                                                                                             )
    )
    )
    appendTab(inputId = "seuratFindMarkerPlotTabset", tabPanel(title = "Heatmap Plot",
                                                                                             panel(heading = "Heatmap Plot",
                                                                                                   fluidRow(
                                                                                                     column(12, align = "center",
                                                                                                            panel(
                                                                                                              plotOutput(outputId = "findMarkerHeatmapPlot")
                                                                                                            )
                                                                                                     )
                                                                                                   )
                                                                                             )
    )

    )

    #singleCellTK:::.exportMetaSlot(vals$counts, "seuratMarkers")

    orderByLFCMarkers <- metadata(vals$counts)$seuratMarkers
    orderByLFCMarkers <- orderByLFCMarkers[order(-orderByLFCMarkers$avg_log2FC), ]
    # vals$fts <- callModule(
    #   module = filterTableServer,
    #   id = "filterSeuratFindMarker",
    #   dataframe = orderByLFCMarkers,
    #   defaultFilterColumns = c("p_val_adj"),
    #   defaultFilterOperators = c("<="),
    #   defaultFilterValues = c("0.05")
    #   )
    vals$fts <- callModule(
      module = filterTableServer,
      id = "filterSeuratFindMarker",
      dataframe = orderByLFCMarkers
    )

  })

  observeEvent(input$findMarkerHeatmapPlotFullNumericRun,{
    ##df <- metadata(vals$counts)$seuratMarkers[which(metadata(vals$counts)$seuratMarkers$p_val_adj < 0.05, arr.ind = TRUE),]
    df <- metadata(vals$counts)$seuratMarkers
    seuratObject <- convertSCEToSeurat(vals$counts, scaledAssay = "seuratScaledData")
    indices <- list()
    cells <- list()
    groups <- unique(colData(vals$counts)[[input$seuratFindMarkerSelectPhenotype]])
    for(i in seq(length(groups))){
      indices[[i]] <- which(colData(vals$counts)[[input$seuratFindMarkerSelectPhenotype]] == groups[i], arr.ind = TRUE)
      cells[[i]] <- colnames(vals$counts)[indices[[i]]]
      cells[[i]] <- lapply(
        X = cells[[i]],
        FUN = function(t) gsub(
          pattern = "_",
          replacement = "-",
          x = t,
          fixed = TRUE)
      )
      Idents(seuratObject, cells = cells[[i]]) <- groups[i]
    }
    colnames(df)[which(startsWith(colnames(df), "avg") == TRUE)] <- "avg_log2FC"
    topMarkers <- df %>% group_by(cluster1) %>% arrange(desc(avg_log2FC)) %>% slice_head(n=input$findMarkerHeatmapPlotFullNumeric)
    #topMarkers <- data.frame(df %>% group_by(cluster1) %>% top_n(input$findMarkerHeatmapPlotFullNumeric, avg_log2FC))
    # if(nrow(topMarkers) > (input$findMarkerHeatmapPlotFullNumeric * length(groups))){
    #   topMarkers <- data.frame(topMarkers %>% group_by(cluster1) %>% top_n(input$findMarkerHeatmapPlotFullNumeric, -p_val_adj))
    # }
    # Subset seuratObject to contain only cells available in selected clusters
    if(input$seuratFindMarkerType != "markerAll"){
      subsetIdents <- c(unique(topMarkers$cluster1), unique(topMarkers$cluster2))
      subsetIdents <- subsetIdents[subsetIdents!="all"]
      seuratObject <- subset(seuratObject, idents = subsetIdents) 
    }
    # Plot heatmap
    output$findMarkerHeatmapPlotFull <- renderPlot({
      isolate({
        DoHeatmap(seuratObject, features = topMarkers$gene.id)
      })
    })
  })

  observe({
    req(vals$fts$data)
    req(vals$fts$selectedRows)
    df <- vals$fts$data[vals$fts$selectedRows, ]
    output$findMarkerRidgePlot <- renderPlot({
      seuratGenePlot(
        inSCE = vals$counts,
        scaledAssayName = "seuratScaledData",
        plotType = "ridge",
        features = df$gene_id,
        groupVariable = input$seuratFindMarkerSelectPhenotype,
        ncol = 2
      )
    })
    output$findMarkerViolinPlot <- renderPlot({
      seuratGenePlot(
        inSCE = vals$counts,
        scaledAssayName = "seuratScaledData",
        plotType = "violin",
        features = df$gene_id,
        groupVariable = input$seuratFindMarkerSelectPhenotype,
        ncol = 2
      )
    })
    output$findMarkerFeaturePlot <- renderPlot({
      seuratGenePlot(
        inSCE = vals$counts,
        scaledAssayName = "seuratScaledData",
        plotType = "feature",
        features = df$gene_id,
        groupVariable = input$seuratFindMarkerSelectPhenotype,
        ncol = 2
      )
    })
    output$findMarkerDotPlot <- renderPlot({
      seuratGenePlot(
        inSCE = vals$counts,
        scaledAssayName = "seuratScaledData",
        plotType = "dot",
        features = df$gene_id,
        groupVariable = input$seuratFindMarkerSelectPhenotype
      )
    })
    output$findMarkerHeatmapPlot <- renderPlot({
      seuratGenePlot(
        inSCE = vals$counts,
        scaledAssayName = "seuratScaledData",
        plotType = "heatmap",
        features = df$gene_id,
        groupVariable = input$seuratFindMarkerSelectPhenotype
      )
    })
  })

  # observe({
  #   req(vals$fts$data)
  #   df <- vals$fts$data
  #   output$findMarkerHeatmapPlotFull <- renderPlot({
  #     seuratGenePlot(
  #       inSCE = vals$counts,
  #       scaledAssayName = "seuratScaledData",
  #       plotType = "heatmap",
  #       features = df$gene_id,
  #       groupVariable = input$seuratFindMarkerSelectPhenotype
  #     )
  #   })
  # })


  #Update PCA/ICA message in clustering tab
  output$display_message_clustering <- renderText({
    if(input$reduction_clustering_method == "pca"){
      if(input$pca_significant_pc_counter){
        paste("<p>Analysis will be performed with <span style='color:red'>", input$pca_significant_pc_counter," components</span> from PCA. This number can be changed in the 'Dimensionality Reduction' section. </p>")
      }
    }
    else{
      if(input$ica_significant_ic_counter){
        paste("<p>Analysis will be performed with <span style='color:red'>", input$ica_significant_ic_counter," components</span> from ICA. This number can be changed in the 'Dimensionality Reduction' section. </p>")
      }
    }
  })

  #Run tSNE
  observeEvent(input$run_tsne_button, {
    req(vals$counts)
    if(!is.null(slot(vals$counts@metadata$seurat$obj, "reductions")[[input$reduction_tsne_method]])){
      withProgress(message = "Running tSNE", max = 1, value = 1, {
        vals$counts <- seuratRunTSNE(inSCE = vals$counts,
                                     useReduction = input$reduction_tsne_method,
                                     reducedDimName = "seuratTSNE",
                                     dims = input$pca_significant_pc_counter,
                                     perplexity = input$perplexity_tsne)
        vals$counts <- singleCellTK:::.seuratInvalidate(inSCE = vals$counts, scaleData = FALSE, varFeatures = FALSE, PCA = FALSE, ICA = FALSE, tSNE = FALSE, UMAP = FALSE)
      })
      withProgress(message = "Plotting tSNE", max = 1, value = 1, {
        output$plot_tsne <- renderPlotly({
          isolate({
            plotly::ggplotly(seuratReductionPlot(inSCE = vals$counts,
                                                 useReduction = "tsne",
                                                 showLegend = FALSE))
          })
        })
      })
      updateCollapse(session = session, "SeuratUI", style = list("tSNE/UMAP" = "success"))
      shinyjs::enable(selector = "div[value='Clustering']")
      S4Vectors::metadata(vals$counts)$seuratMarkers <- NULL
      shinyjs::hide(
        selector = "div[value='Downstream Analysis']")
      showNotification("tSNE Complete")
    }
    else{
      showNotification(paste0("'", input$reduction_tsne_method, "' reduction not found in input object"))
    }
  })


  #Update PCA/ICA message in tSNE tab
  output$display_message_tsne <- renderText({
    if(input$reduction_tsne_method == "pca"){
      if(input$pca_significant_pc_counter){
        paste("<p>Analysis will be performed with <span style='color:red'>", input$pca_significant_pc_counter," components</span> from PCA. This number can be changed in the 'Dimensionality Reduction' section. </p>")
      }
    }
    else{
      if(input$ica_significant_ic_counter){
        paste("<p>Analysis will be performed with <span style='color:red'>", input$ica_significant_ic_counter," components</span> from ICA. This number can be changed in the 'Dimensionality Reduction' section. </p>")
      }
    }
  })

  #Run UMAP
  observeEvent(input$run_umap_button, {
    req(vals$counts)
    if(!is.null(slot(vals$counts@metadata$seurat$obj, "reductions")[[input$reduction_umap_method]])){
      withProgress(message = "Running UMAP", max = 1, value = 1, {
        vals$counts <- seuratRunUMAP(inSCE = vals$counts,
                                     useReduction = input$reduction_umap_method,
                                     reducedDimName = "seuratUMAP",
                                     dims = input$pca_significant_pc_counter,
                                     minDist = input$min_dist_umap,
                                     nNeighbors = input$n_neighbors_umap,
                                     spread = input$spread_umap)
        vals$counts <- singleCellTK:::.seuratInvalidate(inSCE = vals$counts, scaleData = FALSE, varFeatures = FALSE, PCA = FALSE, ICA = FALSE, tSNE = FALSE, UMAP = FALSE)
      })
      withProgress(message = "Plotting UMAP", max = 1, value = 1, {
        output$plot_umap <- renderPlotly({
          isolate({
            plotly::ggplotly(seuratReductionPlot(inSCE = vals$counts,
                                                 useReduction = "umap",
                                                 showLegend = FALSE))
          })
        })
      })
      updateCollapse(session = session, "SeuratUI", style = list("tSNE/UMAP" = "success"))
      shinyjs::enable(selector = "div[value='Clustering']")
      S4Vectors::metadata(vals$counts)$seuratMarkers <- NULL
      shinyjs::hide(
        selector = "div[value='Downstream Analysis']")
      showNotification("UMAP Complete")
    }
    else{
      showNotification(paste0("'", input$reduction_umap_method, "' reduction not found in input object"))
    }
  })

  #Update PCA/ICA message in UMAP tab
  output$display_message_umap <- renderText({
    if(input$reduction_umap_method == "pca"){
      if(input$pca_significant_pc_counter){
        paste("<p>Analysis will be performed with <span style='color:red'>", input$pca_significant_pc_counter," components</span> from PCA. This number can be changed in the 'Dimensionality Reduction' section. </p>")
      }
    }
    else{ #ICA to do
      if(input$ica_significant_ic_counter){
        paste("<p>Analysis will be performed with <span style='color:red'>", input$ica_significant_ic_counter," components</span> from ICA. This number can be changed in the 'Dimensionality Reduction' section. </p>")
      }
    }
  })

  #Update pca significant slider maximum value with total number of computed principal components
  observe({
    req(vals$counts)
    if (!is.null(vals$counts@metadata$seurat$count_pc)) {
      updateSliderInput(session = session, inputId = "pca_significant_pc_counter", max = vals$counts@metadata$seurat$count_pc)
    }
  })

  #Update ica significant slider maximum value with total number of computed independent components
  observe({
    req(vals$counts)
    if (!is.null(vals$counts@metadata$seurat$count_ic)) {
      updateNumericInput(session = session, inputId = "ica_significant_ic_counter", max = vals$counts@metadata$seurat$count_ic)
    }
  })

  #Update tsne, umap and clustering selected number of principal components input
  observe({
    if (input$reduction_umap_method == "pca") {
      updateTextInput(session = session, inputId = "reduction_umap_count", value = input$pca_significant_pc_counter)
    }
    else if (input$reduction_umap_method == "ica") {
      updateTextInput(session = session, inputId = "reduction_umap_count", value = vals$counts@metadata$seurat$count_ic)
    }
    if (input$reduction_clustering_method == "pca") {
      updateTextInput(session = session, inputId = "reduction_clustering_count", value = input$pca_significant_pc_counter)
    }
    else if (input$reduction_clustering_method == "ica") {
      updateTextInput(session = session, inputId = "reduction_clustering_count", value = vals$counts@metadata$seurat$count_ic)
    }
    if (input$reduction_tsne_method == "pca") {
      updateTextInput(session = session, inputId = "reduction_tsne_count", value = input$pca_significant_pc_counter)
    }
    else if (input$reduction_tsne_method == "ica") {
      updateTextInput(session = session, inputId = "reduction_tsne_count", value = vals$counts@metadata$seurat$count_ic)
    }
  })

  #Customize heatmap (pca) with selected options
  observeEvent(input$plot_heatmap_pca_button, {
    if (!is.null(input$picker_dimheatmap_components_pca)) {
      output$plot_heatmap_pca <- renderPlot({
        isolate({
          seuratHeatmapPlot(plotObject = vals$counts@metadata$seurat$heatmap_pca,
                            dims = length(input$picker_dimheatmap_components_pca),
                            ncol = input$slider_dimheatmap_pca,
                            labels = input$picker_dimheatmap_components_pca)
        })
      })
    }
    session$sendCustomMessage("close_dropDownSeuratHM", "")
  })

  #Customize heatmap (ica) with selected options
  observeEvent(input$plot_heatmap_ica_button, {
    if (!is.null(input$picker_dimheatmap_components_ica)) {
      output$plot_heatmap_ica <- renderPlot({
        isolate({
          seuratHeatmapPlot(plotObject = vals$counts@metadata$seurat$heatmap_ica,
                            dims = length(input$picker_dimheatmap_components_ica),
                            ncol = input$slider_dimheatmap_ica,
                            labels = input$picker_dimheatmap_components_ica)
        })
      })
    }
  })


  #Disable tabs & reset collapse panel tabs
  observe({
    if(!is.null(vals$counts)){
      #If data is uploaded in data tab, enable first tab i.e. Normalization tab in Seurat workflow
      shinyjs::enable(
        selector = "div[value='Normalize Data']")

      #Proceed only if sce object has metadata slot
      if(!is.null(vals$counts@metadata)){

        #Proceed only if sce object has seurat object stored in metadata slot
        if(!is.null(vals$counts@metadata$seurat$obj)){
          #If seuratScaledData has been removed from sce object, reset Scale Data tab and reset/lock its next tab
          if(!"seuratScaledData" %in% assayNames(vals$counts)){
            updateCollapse(session = session, "SeuratUI", style = list("Scale Data" = "primary"))
            updateCollapse(session = session, "SeuratUI", style = list("Highly Variable Genes" = "primary"))
            shinyjs::disable(selector = "div[value='Highly Variable Genes']")
          }

          #If variableFeatures have been removed from sce object, reset HVG tab and reset/lock next tab
          if(length(slot(vals$counts@metadata$seurat$obj, "assays")[["RNA"]]@var.features) <= 0){
            updateCollapse(session = session, "SeuratUI", style = list("Highly Variable Genes" = "primary"))
            updateCollapse(session = session, "SeuratUI", style = list("Dimensionality Reduction" = "primary"))
            shinyjs::disable(selector = "div[value='Dimensionality Reduction']")
          }

          #Proceed if reduction slot is present in seurat object in metadata slot
          if("reductions" %in% slotNames(vals$counts@metadata$seurat$obj)){

            #If PCA and ICA both removed from sce object, reset DR tab and reset/lock next tab
            if(is.null(slot(vals$counts@metadata$seurat$obj, "reductions")[["pca"]])
               && is.null(slot(vals$counts@metadata$seurat$obj, "reductions")[["ica"]])){
              updateCollapse(session = session, "SeuratUI", style = list("Dimensionality Reduction" = "primary"))
              updateCollapse(session = session, "SeuratUI", style = list("tSNE/UMAP" = "primary"))
              shinyjs::disable(selector = "div[value='tSNE/UMAP']")
            }

            #If TSNE and UMAP both removed from sce object, reset tSNE/UMAP tab and reset/lock next tab
            if(is.null(slot(vals$counts@metadata$seurat$obj, "reductions")[["tsne"]])
               && is.null(slot(vals$counts@metadata$seurat$obj, "reductions")[["umap"]])){
              updateCollapse(session = session, "SeuratUI", style = list("tSNE/UMAP" = "primary"))
              updateCollapse(session = session, "SeuratUI", style = list("Clustering" = "primary"))
              shinyjs::disable(selector = "div[value='Clustering']")
            }

            #If seurat cluster information removed from sce object, reset Clustering tab
            if(!"seurat_clusters" %in% names(vals$counts@metadata$seurat$obj@meta.data)){
              updateCollapse(session = session, "SeuratUI", style = list("Clustering" = "primary"))
              updateCollapse(session = session, "SeuratUI", style = list("Find Markers" = "primary"))
              shinyjs::disable(selector = "div[value='Find Markers']")
            }
          }

        }

      }
    }
    else{
      #If no data uploaded in data tab, disabled all tabs and plots.

      #Disable tabs
      shinyjs::disable(
        selector = "div[value='Normalize Data']")
      shinyjs::disable(
        selector = "div[value='Scale Data']")
      shinyjs::disable(
        selector = "div[value='Highly Variable Genes']")
      shinyjs::disable(
        selector = "div[value='Dimensionality Reduction']")
      shinyjs::disable(
        selector = "div[value='tSNE/UMAP']")
      shinyjs::disable(
        selector = "div[value='Clustering']")
      shinyjs::disable(
        selector = "div[value='Scale Data']")
      shinyjs::disable(
        selector = "div[value='Find Markers']")

      #Disable plots inside PCA subtab
      shinyjs::disable(
        selector = ".seurat_pca_plots a[data-value='PCA Plot']")
      shinyjs::disable(
        selector = ".seurat_pca_plots a[data-value='Elbow Plot']")
      shinyjs::disable(
        selector = ".seurat_pca_plots a[data-value='JackStraw Plot']")
      shinyjs::disable(
        selector = ".seurat_pca_plots a[data-value='Heatmap Plot']")

      #Disable plots inside ICA subtab
      shinyjs::disable(
        selector = ".seurat_ica_plots a[data-value='ICA Plot']")
      shinyjs::disable(
        selector = ".seurat_ica_plots a[data-value='Heatmap Plot']")

      #Disabled plots inside Clustering tab
      shinyjs::disable(
        selector = ".seurat_clustering_plots a[data-value='PCA Plot']")
      shinyjs::disable(
        selector = ".seurat_clustering_plots a[data-value='ICA Plot']")
      shinyjs::disable(
        selector = ".seurat_clustering_plots a[data-value='tSNE Plot']")
      shinyjs::disable(
        selector = ".seurat_clustering_plots a[data-value='UMAP Plot']")
    }
  })

  #-----------------------------------------------------------------------------
  # Page: Column Annotation (colData) ####
  #-----------------------------------------------------------------------------

  #populate colData from sce object when uploaded
  observe({
    if(!is.null(vals$counts)){
      if(!is.null(colData(vals$counts))){
        vals$columnAnnotation <- as.data.frame(colData(vals$counts))
      }
    }
  })

  #import colData from local storage
  observeEvent(input$importDataButton_colData, {
    withBusyIndicatorServer("importDataButton_colData",{
      if(!is.null(input$uploadFile_colData)){
        temp <- read.csv(input$uploadFile_colData$datapath, header = TRUE,sep = ",")
        if(nrow(colData(vals$counts)) == nrow(temp)){
          if(input$editorChoiceRadio_colData == "replace"){
            vals$columnAnnotation <- temp
          }
          else{
            x <- as.data.frame(colData(vals$counts))
            y <- as.data.frame(temp)
            commonCols <- intersect(colnames(x), colnames(y))
            x[, commonCols] <- y[,commonCols]
            y[, commonCols] <- NULL
            vals$columnAnnotation <- cbind(x, y)
          }
        }
        else{
          showNotification("Number of rows of the assay and the input colData must be equal", type = "error")
        }
      }
      else{
        showNotification("No file selected to upload", type = "error")
      }
    })

    #Render a warning message if there are unsaved changes to colData
    output$changesWarning_colData <- renderUI({
      HTML("<h5><span style='color:red'> You have made changes to the Cell Annotation data. Select 'Save' to finalize these changes or 'Reset' to discard the changes.</span></h5></br>")
    })
    showNotification("You have made changes to the Cell Annotation data. Select 'Save' to finalize these changes or 'Reset' to discard the changes.", type = "error")
  })

  #Render table with colData
  output$outputColumnAnnotationTable_colData <- renderUI({
    output$colOutTable <- DT::renderDataTable({
      DT::datatable(vals$columnAnnotation,
                    editable = 'cell',
                    options = list(pageLength = 5,
                                   scrollX = TRUE))
    })
    DT::dataTableOutput("colOutTable")
  })

  #create selectinput for selecting attribute with colnames from incoming dataset
  #create selectinput for selecting attribute value
  output$inputSelectAttribute_colData <- renderUI({
    if(!is.null(vals$columnAnnotation)){
      if(ncol(vals$columnAnnotation) > 0){
        selectInput("inputSelectAttribute_colData",
                    label = "select attribute",
                    choices = colnames(vals$columnAnnotation))
      }
    }
  })
  output$inputSelectAttributeDelete_colData <- renderUI({
    if(!is.null(vals$columnAnnotation)){
      if(ncol(vals$columnAnnotation) > 0){
        selectInput("inputSelectAttributeDelete_colData",
                    label = "select attribute to delete",
                    choices = colnames(vals$columnAnnotation))
      }
    }
  })

  #create selectinput for selecting column to delete
  output$inputSelectAttributeValue_colData <- renderUI({
    if(!is.null(vals$columnAnnotation) &&
       ncol(vals$columnAnnotation) > 0 &&
       !is.null(input$inputSelectAttribute_colData) &&
       input$inputSelectAttribute_colData %in% colnames(vals$columnAnnotation)){
      selectInput("inputSelectAttributeValue_colData",
                  label = "select attribute value",
                  choices = vals$columnAnnotation[, input$inputSelectAttribute_colData])
    }
  })

  #create selectinput for selecting merge_1 attribute
  #create selectinput for selecting merge_2 attribute
  output$inputSelectAttributeMerge1_colData <- renderUI({
    if(!is.null(vals$columnAnnotation)){
      if(ncol(vals$columnAnnotation) > 0){
        selectInput("inputSelectAttributeMerge1_colData",
                    label = "select first column",
                    choices = colnames(vals$columnAnnotation))
      }
    }
  })
  output$inputSelectAttributeMerge2_colData <- renderUI({
    if(!is.null(vals$columnAnnotation)){
      if(ncol(vals$columnAnnotation) > 0){
        selectInput("inputSelectAttributeMerge2_colData",
                    label = "select second column",
                    choices = colnames(vals$columnAnnotation))
      }
    }
  })

  #create selectinput for selecting fill_1 attribute
  #create selectinput for selecting fill_2 attribute
  output$inputSelectAttributeFill1_colData <- renderUI({
    if(!is.null(vals$columnAnnotation)){
      if(ncol(vals$columnAnnotation) > 0){
        selectInput("inputSelectAttributeFill1_colData",
                    label = "select attribute column",
                    choices = colnames(vals$columnAnnotation))
      }
    }
  })
  output$inputSelectAttributeFill2_colData <- renderUI({
    if(!is.null(vals$columnAnnotation)){
      if(ncol(vals$columnAnnotation) > 0){
        selectInput("inputSelectAttributeFill2_colData",
                    label = "select column to fill",
                    choices = colnames(vals$columnAnnotation))
      }
    }
  })

  #create selectinput for selecting attribute value for magic fill
  output$inputSelectAttributeFillvalue_colData <- renderUI({
    if(!is.null(vals$columnAnnotation)){
      if(ncol(vals$columnAnnotation) > 0){
        selectInput("inputSelectAttributeFillvalue_colData",
                    label = "select attribute value",
                    choices = vals$columnAnnotation[, match(input$inputSelectAttributeFill1_colData,
                                                            colnames(vals$columnAnnotation))])
      }
    }
  })

  #update criteria parameter text input when attribute value selectinput is changed
  observeEvent(input$inputSelectAttributeValue_colData, {
    updateTextInput(session = session,
                    "inputCriteria_colData",
                    value = input$inputSelectAttributeValue_colData)
  })

  #create selectinput for selecting attribute for clean operation
  output$inputSelectAttributeClean_colData <- renderUI({
    if(!is.null(vals$columnAnnotation)){
      if(ncol(vals$columnAnnotation) > 0){
        selectInput("inputSelectAttributeClean_colData",
                    label = "select attribute column",
                    choices = colnames(vals$columnAnnotation))
      }
    }
  })

  #confirm create bin button
  observeEvent(input$buttonConfirmBin_colData, {
    #getting variables
    selected_attribute <- input$inputSelectAttribute_colData
    bin_name <- input$inputBinName_colData
    selected_column_no <- match(selected_attribute, colnames(vals$columnAnnotation))
    criteria_term <- input$inputCriteria_colData
    operator_term <- input$inputOperator_colData

    #get df from reactive input, backup column datatypes and convert factor to character
    data <- singleCellTK:::.manageFactor(vals$columnAnnotation, operation = "backup")
    df <- data$df

    #perform operations
    if (operator_term == "=")
    {
      df[, selected_column_no][df[, selected_column_no] %in% criteria_term] <- bin_name
    }
    else if (operator_term == ">")
    {
      df[, selected_column_no][as.numeric(df[, selected_column_no]) > criteria_term] <- bin_name
    }
    else if (operator_term == "<")
    {
      df[, selected_column_no][as.numeric(df[, selected_column_no]) < criteria_term] <- bin_name
    }
    else if (operator_term == "<=")
    {
      df[, selected_column_no][as.numeric(df[, selected_column_no]) <= criteria_term] <- bin_name
    }
    else if (operator_term == ">=")
    {
      df[, selected_column_no][as.numeric(df[, selected_column_no]) >= criteria_term] <- bin_name
    }

    #restore datatypes
    data$df <- df
    data <- singleCellTK:::.manageFactor(data, operation = "restore")
    vals$columnAnnotation <- data$df

    output$changesWarning_colData <- renderUI({
      HTML("<h5><span style='color:red'> You have made changes to the Cell Annotation data. Select 'Save' to finalize these changes or 'Reset' to discard the changes.</span></h5></br>")
    })
    showNotification("You have made changes to the Cell Annotation data. Select 'Save' to finalize these changes or 'Reset' to discard the changes.", type = "error")
  })

  #confirm merge button
  observeEvent(input$buttonConfirmMerge_colData, {
    df <- vals$columnAnnotation
    colname1 <- input$inputSelectAttributeMerge1_colData
    colname2 <- input$inputSelectAttributeMerge2_colData
    df <- unite_(df, col = colname1, c(colname1, colname2),
                 sep = input$inputSelectSeparatorMerge_colData)

    vals$columnAnnotation <- df

    output$changesWarning_colData <- renderUI({
      HTML("<h5><span style='color:red'> You have made changes to the Cell Annotation data. Select 'Save' to finalize these changes or 'Reset' to discard the changes.</span></h5></br>")
    })
    showNotification("You have made changes to the Cell Annotation data. Select 'Save' to finalize these changes or 'Reset' to discard the changes.", type = "error")
  })

  #fill column button
  observeEvent(input$buttonConfirmFill_colData, {
    #get df from reactive input, backup column datatypes and convert factor to character
    data <- singleCellTK:::.manageFactor(vals$columnAnnotation, operation = "backup")
    df <- data$df

    #perform operation
    selected_attribute_1 <- input$inputSelectAttributeFill1_colData
    selected_attribute_2 <- input$inputSelectAttributeFill2_colData
    selected_column_no_1 <- match(selected_attribute_1, colnames(df))
    selected_column_no_2 <- match(selected_attribute_2, colnames(df))
    old_value <- input$inputSelectAttributeFillvalue_colData
    new_value <- input$inputReplaceText_colData

    for (i in 1:nrow(df))
    {
      if (df[i, selected_column_no_1] == old_value)
      {
        df[i, selected_column_no_2] <- new_value
      }
    }

    #restore datatypes
    data$df <- df
    data <- singleCellTK:::.manageFactor(data, operation = "restore")
    vals$columnAnnotation <- data$df

    output$changesWarning_colData <- renderUI({
      HTML("<h5><span style='color:red'> You have made changes to the Cell Annotation data. Select 'Save' to finalize these changes or 'Reset' to discard the changes.</span></h5></br>")
    })
    showNotification("You have made changes to the Cell Annotation data. Select 'Save' to finalize these changes or 'Reset' to discard the changes.", type = "error")
  })

  #confirm clean button
  observeEvent(input$buttonConfirmClean_colData, {
    #get df from reactive input, backup column datatypes and convert factor to character
    data <- singleCellTK:::.manageFactor(vals$columnAnnotation, operation = "backup")
    df <- data$df

    #perform operation
    selected_attribute <- input$inputSelectAttributeClean_colData
    selected_column_no <- match(selected_attribute, colnames(df))
    selected_choice <- input$inputRemovalOperation_colData
    selected_choice_no <- match(selected_choice, c("remove alphabets",
                                                   "remove digits",
                                                   "remove spaces",
                                                   "remove symbols"))

    if (selected_choice_no == 1)
    {
      for (i in 1:nrow(df))
      {
        df[i, selected_column_no] <- gsub("[A-z]", "", df[i, selected_column_no])
      }

    }
    else if (selected_choice_no == 2)
    {
      for (i in 1:nrow(df))
      {
        df[i, selected_column_no] <- gsub("[0-9]", "", df[i, selected_column_no])
      }
    }
    else if (selected_choice_no == 3)
    {
      for (i in 1:nrow(df))
      {
        df[i, selected_column_no] <- gsub(" ", "", df[i, selected_column_no])
      }
    }
    else if (selected_choice_no == 4)
    {
      for (i in 1:nrow(df))
      {
        df[i, selected_column_no] <- gsub("[[:punct:]]", "", df[i, selected_column_no])
      }
    }

    #restore datatypes
    data$df <- df
    data <- singleCellTK:::.manageFactor(data, operation = "restore")
    vals$columnAnnotation <- data$df

    output$changesWarning_colData <- renderUI({
      HTML("<h5><span style='color:red'> You have made changes to the Cell Annotation data. Select 'Save' to finalize these changes or 'Reset' to discard the changes.</span></h5></br>")
    })
    showNotification("You have made changes to the Cell Annotation data. Select 'Save' to finalize these changes or 'Reset' to discard the changes.", type = "error")
  })

  #add new empty column button
  observeEvent(input$buttonConfirmEmptyColumnName_colData, {
    df <- vals$columnAnnotation

    colname <- input$inputEmptyColumnName_colData
    df$newcolumn <- input$inputDefaultValueAddColumn_colData
    names(df)[ncol(df)] <- colname

    vals$columnAnnotation <- df

    output$changesWarning_colData <- renderUI({
      HTML("<h5><span style='color:red'> You have made changes to the Cell Annotation data. Select 'Save' to finalize these changes or 'Reset' to discard the changes.</span></h5></br>")
    })
    showNotification("You have made changes to the Cell Annotation data. Select 'Save' to finalize these changes or 'Reset' to discard the changes.", type = "error")
  })

  #delete column button
  observeEvent(input$buttonConfirmDeleteColumn_colData,{

    #getting variables
    selected_attribute <- input$inputSelectAttributeDelete_colData

    #get df from reactive input
    df <- vals$columnAnnotation

    #delete
    df[[selected_attribute]] <- NULL

    vals$columnAnnotation <- df

    output$changesWarning_colData <- renderUI({
      HTML("<h5><span style='color:red'> You have made changes to the Cell Annotation data. Select 'Save' to finalize these changes or 'Reset' to discard the changes.</span></h5></br>")
    })
    showNotification("You have made changes to the Cell Annotation data. Select 'Save' to finalize these changes or 'Reset' to discard the changes.", type = "error")
  })

  #restore saved/original colData
  observeEvent(input$buttonRestore_colData,{
    vals$columnAnnotation <- as.data.frame(colData(vals$counts))
    output$changesWarning_colData <- NULL
    showNotification("Changes reverted back to last checkpoint.")
  })

  #save changes to colData
  observeEvent(input$buttonSave_colData,{
    colData(vals$counts) <- DataFrame(vals$columnAnnotation)
    output$changesWarning_colData <- NULL
    updateColDataNames()
    showNotification("Changes saved successfully.")
  })

  #-----------------------------------------------------------------------------
  # Page: Row Annotation (rowData) ####
  #-----------------------------------------------------------------------------

  #populate colData from sce object when uploaded
  observe({
    if(!is.null(vals$counts)){
      if(!is.null(rowData(vals$counts))){
        vals$rowAnnotation <- as.data.frame(rowData(vals$counts))
      }
    }
  })

  #import rowData from local storage
  observeEvent(input$importDataButton_rowData, {
    withBusyIndicatorServer("importDataButton_rowData",{
      if(!is.null(input$uploadFile_rowData)){
        temp <- read.csv(input$uploadFile_rowData$datapath, header = TRUE,sep = ",")
        if(nrow(rowData(vals$counts)) == nrow(temp)){
          if(input$editorChoiceRadio_rowData == "replace"){
            vals$rowAnnotation <- temp
          }
          else{
            x <- as.data.frame(rowData(vals$counts))
            y <- as.data.frame(temp)
            commonCols <- intersect(colnames(x), colnames(y))
            x[, commonCols] <- y[,commonCols]
            y[, commonCols] <- NULL
            vals$rowAnnotation <- cbind(x, y)
          }
        }
        else{
          showNotification("Number of rows of the assay and the input rowData must be equal", type = "error")
        }
      }
      else{
        showNotification("No file selected to upload", type = "error")
      }
    })

    #Render a warning message if there are unsaved changes to rowData
    output$changesWarning_rowData <- renderUI({
      HTML("<h5><span style='color:red'> You have made changes to the Cell Annotation data. Select 'Save' to finalize these changes or 'Reset' to discard the changes.</span></h5></br>")
    })
    showNotification("You have made changes to the Cell Annotation data. Select 'Save' to finalize these changes or 'Reset' to discard the changes.", type = "error")
  })

  #Render table with rowData
  output$outputColumnAnnotationTable_rowData <- renderUI({
    output$rowOutTable <- DT::renderDataTable({
      DT::datatable(vals$rowAnnotation,
                    editable = 'cell',
                    options = list(pageLength = 5,
                                   scrollX = TRUE))
    })
    DT::dataTableOutput("rowOutTable")
  })

  #create selectinput for selecting attribute with colnames from incoming dataset
  #create selectinput for selecting attribute value
  output$inputSelectAttribute_rowData <- renderUI({
    if(!is.null(vals$rowAnnotation)){
      if(ncol(vals$rowAnnotation) > 0){
        selectInput("inputSelectAttribute_rowData",
                    label = "select attribute",
                    choices = colnames(vals$rowAnnotation))
      }
    }
  })
  output$inputSelectAttributeDelete_rowData <- renderUI({
    if(!is.null(vals$rowAnnotation)){
      if(ncol(vals$rowAnnotation) > 0){
        selectInput("inputSelectAttributeDelete_rowData",
                    label = "select attribute to delete",
                    choices = colnames(vals$rowAnnotation))
      }
    }
  })

  #create selectinput for selecting column to delete
  observeEvent(input$inputSelectAttribute_rowData, {
    if(!is.null(vals$rowAnnotation) &&
       ncol(vals$rowAnnotation) > 0 &&
       !is.null(input$inputSelectAttribute_rowData) &&
       input$inputSelectAttribute_rowData %in% colnames(vals$rowAnnotation)){
      updateSelectizeInput(session, "inputSelectAttributeValue_rowData",
                           choices = vals$rowAnnotation[, input$inputSelectAttribute_rowData],
                           server = TRUE)
    }
  })

  #create selectinput for selecting merge_1 attribute
  #create selectinput for selecting merge_2 attribute
  output$inputSelectAttributeMerge1_rowData <- renderUI({
    if(!is.null(vals$rowAnnotation)){
      if(ncol(vals$rowAnnotation) > 0){
        selectInput("inputSelectAttributeMerge1_rowData",
                    label = "select first column",
                    choices = colnames(vals$rowAnnotation))
      }
    }
  })
  output$inputSelectAttributeMerge2_rowData <- renderUI({
    if(!is.null(vals$rowAnnotation)){
      if(ncol(vals$rowAnnotation) > 0){
        selectInput("inputSelectAttributeMerge2_rowData",
                    label = "select second column",
                    choices = colnames(vals$rowAnnotation))
      }
    }
  })

  #create selectinput for selecting fill_1 attribute
  #create selectinput for selecting fill_2 attribute
  output$inputSelectAttributeFill1_rowData <- renderUI({
    if(!is.null(vals$rowAnnotation)){
      if(ncol(vals$rowAnnotation) > 0){
        selectInput("inputSelectAttributeFill1_rowData",
                    label = "select attribute column",
                    choices = colnames(vals$rowAnnotation))
      }
    }
  })
  output$inputSelectAttributeFill2_rowData <- renderUI({
    if(!is.null(vals$rowAnnotation)){
      if(ncol(vals$rowAnnotation) > 0){
        selectInput("inputSelectAttributeFill2_rowData",
                    label = "select column to fill",
                    choices = colnames(vals$rowAnnotation))
      }
    }
  })

  #create selectinput for selecting attribute value for magic fill
  observeEvent(input$inputSelectAttributeFill1_rowData, {
    if(!is.null(vals$rowAnnotation)){
      if(ncol(vals$rowAnnotation) > 0){
        updateSelectizeInput(session, "inputSelectAttributeFillvalue_rowData",
                             choices = vals$rowAnnotation[, match(input$inputSelectAttributeFill1_rowData,
                                                                  colnames(vals$rowAnnotation))],
                             server = TRUE)
      }
    }
  })

  #update criteria parameter text input when attribute value selectinput is changed
  observeEvent(input$inputSelectAttributeValue_rowData, {
    updateTextInput(session = session,
                    "inputCriteria_rowData",
                    value = input$inputSelectAttributeValue_rowData)
  })

  #create selectinput for selecting attribute for clean operation
  output$inputSelectAttributeClean_rowData <- renderUI({
    if(!is.null(vals$rowAnnotation)){
      if(ncol(vals$rowAnnotation) > 0){
        selectInput("inputSelectAttributeClean_rowData",
                    label = "select attribute column",
                    choices = colnames(vals$rowAnnotation))
      }
    }
  })

  #confirm create bin button
  observeEvent(input$buttonConfirmBin_rowData, {
    #getting variables
    selected_attribute <- input$inputSelectAttribute_rowData
    bin_name <- input$inputBinName_rowData
    selected_column_no <- match(selected_attribute, colnames(vals$rowAnnotation))
    criteria_term <- input$inputCriteria_rowData
    operator_term <- input$inputOperator_rowData

    #get df from reactive input, backup column datatypes and convert factor to character
    data <- singleCellTK:::.manageFactor(vals$rowAnnotation, operation = "backup")
    df <- data$df

    #operations
    if (operator_term == "=")
    {
      df[, selected_column_no][df[, selected_column_no] %in% criteria_term] <- bin_name
    }
    else if (operator_term == ">")
    {
      df[, selected_column_no][as.numeric(df[, selected_column_no]) > criteria_term] <- bin_name
    }
    else if (operator_term == "<")
    {
      df[, selected_column_no][as.numeric(df[, selected_column_no]) < criteria_term] <- bin_name
    }
    else if (operator_term == "<=")
    {
      df[, selected_column_no][as.numeric(df[, selected_column_no]) <= criteria_term] <- bin_name
    }
    else if (operator_term == ">=")
    {
      df[, selected_column_no][as.numeric(df[, selected_column_no]) >= criteria_term] <- bin_name
    }

    #restore datatypes
    data$df <- df
    data <- singleCellTK:::.manageFactor(data, operation = "restore")
    vals$rowAnnotation <- data$df

    output$changesWarning_rowData <- renderUI({
      HTML("<h5><span style='color:red'> You have made changes to the Cell Annotation data. Select 'Save' to finalize these changes or 'Reset' to discard the changes.</span></h5></br>")
    })
    showNotification("You have made changes to the Cell Annotation data. Select 'Save' to finalize these changes or 'Reset' to discard the changes.", type = "error")
  })

  #confirm merge button
  observeEvent(input$buttonConfirmMerge_rowData, {
    df <- vals$rowAnnotation
    colname1 <- input$inputSelectAttributeMerge1_rowData
    colname2 <- input$inputSelectAttributeMerge2_rowData
    df <- unite_(df, col = colname1, c(colname1, colname2),
                 sep = input$inputSelectSeparatorMerge_rowData)

    vals$rowAnnotation <- df

    output$changesWarning_rowData <- renderUI({
      HTML("<h5><span style='color:red'> You have made changes to the Cell Annotation data. Select 'Save' to finalize these changes or 'Reset' to discard the changes.</span></h5></br>")
    })
    showNotification("You have made changes to the Cell Annotation data. Select 'Save' to finalize these changes or 'Reset' to discard the changes.", type = "error")
  })

  #fill column button
  observeEvent(input$buttonConfirmFill_rowData, {
    #get df from reactive input, backup column datatypes and convert factor to character
    data <- singleCellTK:::.manageFactor(vals$rowAnnotation, operation = "backup")
    df <- data$df

    #operations
    selected_attribute_1 <- input$inputSelectAttributeFill1_rowData
    selected_attribute_2 <- input$inputSelectAttributeFill2_rowData
    selected_column_no_1 <- match(selected_attribute_1, colnames(df))
    selected_column_no_2 <- match(selected_attribute_2, colnames(df))
    old_value <- input$inputSelectAttributeFillvalue_rowData
    new_value <- input$inputReplaceText_rowData

    for (i in 1:nrow(df))
    {
      if (df[i, selected_column_no_1] == old_value)
      {
        df[i, selected_column_no_2] <- new_value
      }
    }

    #restore datatypes
    data$df <- df
    data <- singleCellTK:::.manageFactor(data, operation = "restore")
    vals$rowAnnotation <- data$df

    output$changesWarning_rowData <- renderUI({
      HTML("<h5><span style='color:red'> You have made changes to the Cell Annotation data. Select 'Save' to finalize these changes or 'Reset' to discard the changes.</span></h5></br>")
    })
    showNotification("You have made changes to the Cell Annotation data. Select 'Save' to finalize these changes or 'Reset' to discard the changes.", type = "error")
  })

  #confirm clean button
  observeEvent(input$buttonConfirmClean_rowData, {
    #get df from reactive input, backup column datatypes and convert factor to character
    data <- singleCellTK:::.manageFactor(vals$rowAnnotation, operation = "backup")
    df <- data$df

    #operations
    selected_attribute <- input$inputSelectAttributeClean_rowData
    selected_column_no <- match(selected_attribute, colnames(df))
    selected_choice <- input$inputRemovalOperation_rowData
    selected_choice_no <- match(selected_choice, c("remove alphabets",
                                                   "remove digits",
                                                   "remove spaces",
                                                   "remove symbols"))

    if (selected_choice_no == 1)
    {
      for (i in 1:nrow(df))
      {
        df[i, selected_column_no] <- gsub("[A-z]", "", df[i, selected_column_no])
      }

    }
    else if (selected_choice_no == 2)
    {
      for (i in 1:nrow(df))
      {
        df[i, selected_column_no] <- gsub("[0-9]", "", df[i, selected_column_no])
      }
    }
    else if (selected_choice_no == 3)
    {
      for (i in 1:nrow(df))
      {
        df[i, selected_column_no] <- gsub(" ", "", df[i, selected_column_no])
      }
    }
    else if (selected_choice_no == 4)
    {
      for (i in 1:nrow(df))
      {
        df[i, selected_column_no] <- gsub("[[:punct:]]", "", df[i, selected_column_no])
      }
    }

    #restore datatypes
    data$df <- df
    data <- singleCellTK:::.manageFactor(data, operation = "restore")
    vals$rowAnnotation <- data$df

    output$changesWarning_rowData <- renderUI({
      HTML("<h5><span style='color:red'> You have made changes to the Cell Annotation data. Select 'Save' to finalize these changes or 'Reset' to discard the changes.</span></h5></br>")
    })
    showNotification("You have made changes to the Cell Annotation data. Select 'Save' to finalize these changes or 'Reset' to discard the changes.", type = "error")
  })

  #add new empty column button
  observeEvent(input$buttonConfirmEmptyColumnName_rowData, {
    df <- vals$rowAnnotation
    colname <- input$inputEmptyColumnName_rowData
    df$newcolumn <- input$inputDefaultValueAddColumn_rowData
    names(df)[ncol(df)] <- colname

    vals$columnAnnotation <- df

    output$changesWarning_rowData <- renderUI({
      HTML("<h5><span style='color:red'> You have made changes to the Cell Annotation data. Select 'Save' to finalize these changes or 'Reset' to discard the changes.</span></h5></br>")
    })
    showNotification("You have made changes to the Cell Annotation data. Select 'Save' to finalize these changes or 'Reset' to discard the changes.", type = "error")
  })

  #delete column button
  observeEvent(input$buttonConfirmDeleteColumn_rowData,{

    #getting variables
    selected_attribute <- input$inputSelectAttributeDelete_rowData

    #get df from reactive input
    df <- vals$rowAnnotation

    #delete
    df[[selected_attribute]] <- NULL

    vals$rowAnnotation <- df

    output$changesWarning_rowData <- renderUI({
      HTML("<h5><span style='color:red'> You have made changes to the Cell Annotation data. Select 'Save' to finalize these changes or 'Reset' to discard the changes.</span></h5></br>")
    })
    showNotification("You have made changes to the Cell Annotation data. Select 'Save' to finalize these changes or 'Reset' to discard the changes.", type = "error")
  })

  #restore saved/original rowData
  observeEvent(input$buttonRestore_rowData,{
    vals$rowAnnotation <- as.data.frame(rowData(vals$counts))
    output$changesWarning_rowData <- NULL
    showNotification("Changes reverted back to last checkpoint.")
  })

  #save changes to rowData
  observeEvent(input$buttonSave_rowData,{
    rowData(vals$counts) <- DataFrame(vals$rowAnnotation)
    output$changesWarning_rowData <- NULL
    updateFeatureAnnots()
    showNotification("Changes saved successfully.")
  })


  #-----------------------------------------------------------------------------
  # Page Download ####
  #-----------------------------------------------------------------------------

  exportPath = '~'
  shinyDirChoose(input, 'outputDirectory', roots = roots)
  output$outputDirectoryPath <- renderText({
    dirPaths$outputDirectory
  })
  observeEvent(
    ignoreNULL = TRUE,
    eventExpr = {
      input$outputDirectory
    },
    handlerExpr = {
      if ("path" %in% names(input$outputDirectory)) {
        # condition prevents handler execution on initial app launch
        #path <<- choose.dir(default = readDirectoryInput(session, 'outputDirectory'))
        #updateDirectoryInput(session, 'outputDirectory', value = path)
        vol <- roots[[input$outputDirectory$root]]
        dirPaths$outputDirectory <- paste0(vol, paste(unlist(input$outputDirectory$path[-1]),
                                             collapse = .Platform$file.sep))
        exportPath <<- dirPaths$outputDirectory
      }
    }
  )

  output$exportFileName <- renderUI({
    defaultName <- paste0("SCE-", strftime(Sys.time(), format = "%y%m%d_%H%M"))
    if (input$exportChoice == "rds") {
      extName <- ".rds"
    } else if (input$exportChoice == "annData") {
      extName <- ".h5ad"
    } else if (input$exportChoice == "textfile") {
      extName <- ".txt"
    }
    if (input$exportChoice != "textfile") {
      tags$div(
        div(style = "display: inline-block;vertical-align:top; width: 160px;",
            textInput("exportPrefix", label = NULL,
                      value = defaultName, placeholder = "Required!",
                      width = '160px')),
        div(
          style = "display: inline-block;vertical-align:top; width: 50px;",
          p(extName, style = "margin-top: 8px; margin-left: 2px; font-size: 16px;")
        )
      )
    } else {
      tags$div(
        div(style = "display: inline-block;vertical-align:top; width: 160px;",
            textInput("exportPrefix", label = NULL,
                      value = defaultName, placeholder = "Required!",
                      width = '160px')),
      )
    }

  })

  addPopover(session, 'exportAssayLabel', '', "The name of assay of interests that will be set as the primary matrix of the output AnnData.", 'right')
  addPopover(session, 'compressionLabel', '', "If output file compression is required, this variable accepts 'gzip' or 'lzf' as inputs", 'right')
  addPopover(session, 'compressionOptsLabel', '', "Sets the compression level", 'right')
  addPopover(session, 'forceDenseLabel', '', "Default False. Write sparse data as a dense matrix. Refer anndata.write_h5ad documentation for details.", 'right')

  addPopover(session, 'gzipLabel', '', 'Set to true if output files are to be gzip compressed', 'right')
  addPopover(session, 'overwriteLabel', '', 'Overwrites the file if it already exists', 'right')

  observeEvent(input$exportData, {
    withBusyIndicatorServer("exportData", {
      if (is.null(vals$counts) && is.null(vals$original)) {
        shinyalert::shinyalert("Error!", "Upload data first.", type = "error")
      } else {
        if (input$exportChoice == "rds") {
          filename <- paste0(input$exportPrefix, ".rds")
          saveRDS(vals$counts, paste0(exportPath, "/", filename))
        } else if (input$exportChoice == "annData") {
          exportSCEtoAnnData(sce=vals$counts,
                             useAssay = input$exportAssay,
                             outputDir = exportPath,
                             prefix = input$exportPrefix,
                             overwrite = input$exportOverwrite,
                             compression = "gzip",
                             compressionOpts = input$compressionOpts,
                             forceDense = input$forceDense)
        } else if (input$exportChoice == "textfile") {
          exportSCEtoFlatFile(sce = vals$counts,
                              outputDir = exportPath,
                              overwrite = input$exportOverwrite,
                              gzipped = input$exportFlatGzip,
                              prefix = input$exportPrefix)
        }
      }
    })
  })

  ##############################################################################
  # Page: Cell Type Labeling ####
  ##############################################################################
  output$ctLabelLevelUI <- renderUI({
    if (input$ctLabelRef %in% c("hpca", "bpe", "dice", "immgen", "mouse")) {
      selectInput("ctLabelLevel", "Labeling level:",
                  c("main", "fine", "ont"), "main")
    } else {
      disabled(
        selectInput("ctLabelLevel", "Labeling level (not supported):",
                    choices = NULL, selected = NULL)
      )
    }
  })

  observeEvent(input$ctLabelRun, {
    if (!is.null(vals$counts)) {
      withBusyIndicatorServer("ctLabelRun", {
        if (input$ctLabelBy == "Clusters") {
          cluster <- input$ctLabelByCluster
          if (is.null(cluster)) {
            stop("Choose the clustering label for this condition!")
          }
        } else {
          cluster <- NULL
        }
        vals$counts <- runSingleR(vals$counts,
                                  useAssay = input$ctLabelAssay,
                                  useBltinRef = input$ctLabelRef,
                                  level = input$ctLabelLevel,
                                  featureType = input$ctLabelFeatureType,
                                  labelByCluster = cluster)
        updateColDataNames()
        shinyalert(
          title = "Labeled!",
          type = "success",
          text = "Cell type labeling stored as cell annotations. Users can visualize via CellViewer."
        )
      })
    }
  })

  ##############################################################################
  # Code for ShinyTest ####
  ##############################################################################
  observe({
    shinyBS::updateCollapse(session,
                            "SeuratUI",
                            open = input$activePanelSelectSeurat)
  })

  ##############################################################################
  # Code for PushBar ####
  ##############################################################################
  # observeEvent(input$interpretToggle, {
  #   pushbar_open(id = "myPushbar")
  # })
})

