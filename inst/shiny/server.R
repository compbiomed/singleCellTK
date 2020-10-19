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
# source("server_partials/server_01_data.R", local = TRUE) # functions for Data section

# Define server logic required to draw a histogram
shinyServer(function(input, output, session) {
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
    hvgCalculated = list(status = FALSE, method = NULL)
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
    updateSelectInput(session, "scMergeCT",
                      choices = c('None', pdataOptions))
    updateSelectInput(session, "combatCond",
                      choices = c('None', pdataOptions))
    updateSelectInput(session, "batchCorrVar",
                      choices = pdataOptions)
    updateSelectInput(session, "batchCheckVar",
                      choices = pdataOptions)
    updateSelectInput(session, "batchCheckCond",
                      choices = c('None', pdataOptions))
    updateSelectInput(session, "clustVisCol", choices = pdataOptions)
    updateSelectInput(session, "deC1Class",
                      choices = c('None', pdataOptions))
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
    updateSelectInput(session, "visCondn",
                      choices = c("none", pdataOptions))
    updateSelectInput(session, "hmCellCol",
                      choices = pdataOptions)
    updateSelectInput(session, "hmCellTextBy",
                      choices = c("Row Names", pdataOptions))
    updateSelectInput(session, 'hmAddCellLabel',
                      choices = c("Default cell IDs", pdataOptions))
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
    updateSelectInput(session, "Knumber",
                      choices = 1:numsamples)
    updateSelectInput(session, "Cnumber",
                      choices = 1:numsamples)
    # updateSelectInput(session, "pcX",
    #                   choices = paste("PC", 1:numsamples, sep = ""),
    #                   selected = "PC1")
    # updateSelectInput(session, "pcY",
    #                   choices = paste("PC", 1:numsamples, sep = ""),
    #                   selected = "PC2")
    updateNumericInput(session, "downsampleNum", value = numsamples,
                       max = numsamples)
  }

  updateAssayInputs <- function(){
    currassays <- names(assays(vals$counts))
    updateSelectInput(session, "dimRedAssaySelect", choices = currassays)
    updateSelectInput(session, "batchCorrAssay", choices = currassays)
    updateSelectInput(session, "batchCheckAssay", choices = currassays)
    updateSelectInput(session, "batchCheckOrigAssay", choices = currassays)
    updateSelectInput(session, "deAssay", choices = currassays)
    updateSelectInput(session, "fmAssay", choices = currassays)
    updateSelectInput(session, "fmHMAssay", choices = currassays)
    updateSelectInput(session, "pathwayAssay", choices = currassays)
    updateSelectInput(session, "modifyAssaySelect", choices = currassays)
    updateSelectInput(session, "normalizeAssaySelect", choices = currassays)
    updateSelectInput(session, "seuratSelectNormalizationAssay", choices = currassays)
    updateSelectInput(session, "assaySelectFS_Norm", choices = currassays)
    updateSelectInput(session, "filterAssaySelect", choices = currassays)
    updateSelectInput(session, "qcAssaySelect", choices = currassays)
    updateSelectInput(session, "visAssaySelect", choices = currassays)
    updateSelectInput(session, "enrichAssay", choices = currassays)
    updateSelectInput(session, "celdaAssay", choices = currassays)
    updateSelectInput(session, "celdaAssayGS", choices = currassays)
    updateSelectInput(session, "celdaAssaytSNE", choices = currassays)
    updateSelectInput(session, "celdaAssayProbabilityMap",
                      choices = currassays)
    updateSelectInput(session, "celdaAssayModuleHeatmap",
                      choices = currassays)
    updateSelectInput(session, "depthAssay", choices = currassays)
    updateSelectInput(session, "cellsAssay", choices = currassays)
    updateSelectInput(session, "snapshotAssay", choices = currassays)
    updateSelectInput(session, "exportAssay", choices = currassays)
    updateSelectInput(session, "hmAssay", choices = currassays)
  }

  observeEvent(vals$counts, {
    # vals$counts
    if (!is.null(vals$counts)) {
      updateAssayInputs()
    }
  })

  observeEvent(vals$original, {
    if (!is.null(vals$original)) {
      if (!is.null(metadata(vals$original)$sctk$genesets)) {
        newGSchoices <- sctkListGeneSetCollections(vals$original)
        updateSelectInput(session, "gsExisting", choices = c("None", newGSchoices))
        updateSelectInput(session, "QCMgeneSets", choices =c("None", newGSchoices))
        shinyjs::show(id = "gsAddToExisting", anim = FALSE)
      } else {
        shinyjs::hide(id = "gsAddToExisting", anim = FALSE)
        updateSelectInput(session, "gsExisting", choices = c("None"), selected = "None")
        updateSelectInput(session, "QCMgeneSets", choices =c("None"), selected = "None")
      }
    }
  })

  updateReddimInputs <- function(){
    currreddim <- names(reducedDims(vals$counts))
    updateSelectInput(session, "delRedDimType", choices = currreddim)
    updateSelectInput(session, "FastMNNReddim", choices = currreddim)
    updateSelectInput(session, "HarmonyReddim", choices = currreddim)
    updateSelectInput(session, "clustVisReddim", choices = currreddim)
    updateSelectInput(session, "clustKMeansReddim", choices = currreddim)
    updateSelectInput(session, "clustSeuratReddim", choices = currreddim)
  }

  updateAltExpInputs <- function(){
    options <- altExpNames(vals$counts)
    updateSelectInput(session, "clustScranSNNAltExp", choices = options)
    updateSelectInput(session, "dimRedAltExpSelect", choices = options)
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
  # Page 1: Upload
  #-----------------------------------------------------------------------------

  # Upload data through shiny app

  # Components for uploading directories if user is importing from a preprocessing step
  volumes <- c(Home = fs::path_home(), "R Installation" = R.home(), shinyFiles::getVolumes()())
  shinyFiles::shinyDirChoose(input, "base", roots = volumes, session = session, restrictions = system.file(package = "base"))
  shinyFiles::shinyDirChoose(input, "sample", roots = volumes, session = session, restrictions = system.file(package = "base"))

  sample <- reactive(input$sample)
  output$sample <- renderText({
    shinyFiles::parseDirPath(volumes, sample())
  })
  sampleFile <- reactive(input$sampleFile)
  output$sampleFile <- renderText({
    shinyFiles::parseFilePaths(volumes, sampleFile())$datapath
  })
  output$base = renderText({
    shinyDirectoryInput::readDirectoryInput(session, 'directory')
  })
  importCR2Files <- reactiveValues(files = list(), id_count = 0)
  importCR3Files <- reactiveValues(files = list(), id_count = 0)
  importSSFiles <- reactiveValues(files = list(), id_count = 0)
  importBUSFiles <- reactiveValues(files = list(), id_count = 0)
  importSEQFiles <- reactiveValues(files = list(), id_count = 0)
  importOptFiles <- reactiveValues(files = list(), id_count = 0)

  allImportEntries <- reactiveValues(samples=list(), id_count=0)


  # modal to import all preprocessed data except for CellRanger data
  importModal <- function(failed=FALSE, needsDir=FALSE) {

    modalDialog(
      h3("Sample Name"),
      textInput("sampleName", "*This is the name you would like to give your sample."),
      # only some functions need this input
      if (needsDir)
        h3("Sample ID"),
      if (needsDir)
        textInput("sampleID", "*This name must match your sample's directory name."),


      h3("Base Directory"),
      shinyDirectoryInput::directoryInput('directory', label = 'Choose Directory', value = '~'),

      if (failed)
        div(tags$b("Please fill out all the required fields", style = "color: red;")),

      footer = tagList(
        modalButton("Cancel"),
        actionButton("modalOk", "OK")
      )
    )
  }


  # modal to import CellRanger data
  importCRModal <- function() {
    modalDialog(
      h3("Add a Cell Ranger Sample"),
      tags$br(),
      h4("Option 1 - Select a directory containing multiple sample directories (and no other directories)."),
      actionButton("crOpt1", "Add"),
      tags$br(),
      h4("Option 2 - Select a single sample directory."),
      actionButton("crOpt2", "Add"),
      tags$br(),
      h4("Option 3 - Select a directory containing your data files (barcodes.tsv, features.tsv, matrix.mtx)."),
      actionButton("crOpt3", "Add"),

      footer = tagList(
        modalButton("Cancel"),
        actionButton("crOK", "OK")
      )
    )
  }
  # Upload a sample directory (CR) (parent of 'outs' directory)
  importCRSDir <- function(failed = FALSE) {
    modalDialog(
      h3("Sample Directory"),
      shinyDirectoryInput::directoryInput('sDirectory', label = 'Choose Directory', value = '~'),
      h3("Sample Name"),
      h5("If you do not provide an alternate sample name, the sample name will be set to the sample directory name."),
      textInput("sSampleID", ""),

      if (failed)
        div(tags$b("Please fill out all the required fields", style = "color: red;")),

      footer = tagList(
        modalButton("Cancel"),
        actionButton("SDirOK", "OK")
      )
    )
  }
  # Upload a data directory (CR) (parent of 'data files')
  importCRDDir <- function(failed = FALSE) {
    modalDialog(
      h3("Data Directory"),
      shinyDirectoryInput::directoryInput('directory', label = 'Choose Directory', value = '~'),
      h3("Sample Name"),

      textInput("dSampleID", "*This field is mandatory when uploading a data directory"),

      if (failed)
        div(tags$b("Please fill out all the required fields", style = "color: red;")),

      footer = tagList(
        modalButton("Cancel"),
        actionButton("DDirOK", "OK")
      )
    )
  }
  # Upload a base directory (CR) (parent of possibly multiple sample directories)
  importCRBDir <- function(failed = FALSE) {
    modalDialog(
      h3("Base Directory"),
      shinyDirectoryInput::directoryInput('bDirectory', label = 'Choose Directory', value = '~'),
      wellPanel(h5("*For any sample names that you do not provide, the sample name will be set to the sample directory name.")),

      tags$div(id = "bDirTable"),

      if (failed)
        div(tags$b("Please fill out all the required fields", style = "color: red;")),

      footer = tagList(
        modalButton("Cancel"),
        actionButton("BDirOK", "OK")
      )
    )
  }


  # see https://github.com/wleepang/shiny-directory-input
  observeEvent(
    ignoreNULL = TRUE,
    eventExpr = {
      input$directory
    },
    handlerExpr = {
      if (input$directory > 0) {
        # condition prevents handler execution on initial app launch
        path = shinyDirectoryInput::choose.dir(default = shinyDirectoryInput::readDirectoryInput(session, 'directory'),
                                               caption="Choose a directory")
        shinyDirectoryInput::updateDirectoryInput(session, 'directory', value = path)
      }
    }
  )

  # see https://github.com/wleepang/shiny-directory-input
  # for sample directory modal
  observeEvent(
    ignoreNULL = TRUE,
    eventExpr = {
      input$sDirectory
    },
    handlerExpr = {
      if (input$sDirectory > 0) {
        # condition prevents handler execution on initial app launch
        path = shinyDirectoryInput::choose.dir(default = shinyDirectoryInput::readDirectoryInput(session, 'sDirectory'),
                                               caption="Choose a directory")
        shinyDirectoryInput::updateDirectoryInput(session, 'sDirectory', value = path)
        if (!is.na(path)) {
          updateTextInput(session, "sSampleID", value = basename(path))
        }
      }
    }
  )

  # event listener for the base directory modal (need to populate table for sample names)
  # see https://github.com/wleepang/shiny-directory-input
  observeEvent(
    ignoreNULL = TRUE,
    eventExpr = {
      input$bDirectory
    },
    handlerExpr = {
      if (input$bDirectory > 0) {
        # condition prevents handler execution on initial app launch
        path = shinyDirectoryInput::choose.dir(default = shinyDirectoryInput::readDirectoryInput(session, 'bDirectory'),
                                               caption="Choose a directory")
        shinyDirectoryInput::updateDirectoryInput(session, 'bDirectory', value = path)
        # clear the previous table of sample names
        prevPath <- shinyDirectoryInput::readDirectoryInput(session, 'bDirectory')
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

  # event listeners for Cell Ranger import modals' OK buttons
  # sample directory
  observeEvent(input$SDirOK, {
    samplePath <- shinyDirectoryInput::readDirectoryInput(session, 'sDirectory')
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
    dataPath <- shinyDirectoryInput::readDirectoryInput(session, 'directory')
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
        entry <- list(type-"cellRanger3", id=id, params=list(dataDir = dataPath, sampleName = input$dSampleID))
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

  # base directory
  observeEvent(input$BDirOK, {
    basePath <- shinyDirectoryInput::readDirectoryInput(session, 'bDirectory')
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
          entry <- list(type="cellRanger2", id=id, params=list(cellRangerDirs = substr(basePath, 1, nchar(basePath)-1), sampleDirs = basename(sample), sampleNames = name))
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
          entry <- list(type="cellRanger3", id=id, params=list(cellRangerDirs = substr(basePath, 1, nchar(basePath)-1), sampleDirs = basename(sample), sampleNames = name))
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
            print(allImportEntries$samples)
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

  # event handler for pressing OK on the import modal
  observeEvent(input$modalOk, {
    samplePath <- shinyFiles::parseDirPath(volumes, input$sample)
    basePath <- shinyDirectoryInput::readDirectoryInput(session, 'directory')
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
        entry <- list(type="busTools", id = id, params=list(BUStoolsDirs = substr(basePath, 1, nchar(basePath)-1), samples = input$sampleName))
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
      if(is.null(colData(vals$original)$sample)) {
        colData(vals$original)$sample = "sample"
      }

      if (!is.null(vals$original)) {
        vals$counts <- vals$original

        # ToDo: Remove these automatic updates and replace with
        # observeEvents functions that activate upon the tab selection
        updateColDataNames()
        updateFeatureAnnots()
        updateNumSamples()
        # updateAssayInputs()
        updateGeneNames()
        updateReddimInputs()
        updateSelectInput(session, "qcSampleSelect", choices = c("None", names(colData(vals$original))), selected = "sample")
        shinyjs::show(id="annotationData")
        js$enableTabs();
      } else {
        shinyalert::shinyalert("Error!", "The data upload failed!",
                               type = "error")
      }
      vals$gsvaRes <- NULL
      vals$gsvaLimma <- NULL
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
    })
  })


  #-----------#
  # Gene Sets #
  #-----------#
  handleGSPasteOption <- function() {
    if (!nzchar(input$geneSetText)) {
      shinyjs::show(id = "gsUploadError", anim = FALSE)
    } else if ((!nzchar(input$gsCollectionNameText)) && (input$gsExisting == "None")) {
      shinyjs::show(id = "gsUploadError", anim = FALSE)
    } else {
      shinyjs::hide(id = "gsUploadError", anim = FALSE)
      setList <- formatGeneSetList(input$geneSetText)
      if (!nzchar(input$gsCollectionNameText)) {
        vals$original <- importGeneSetsFromList(vals$original, setList, collectionName = input$gsCollectionNameText)
      } else if (input$gsExisting != "None") {
        vals$original <- importGeneSetsFromList(vals$original, setList, collectionName = input$gsExisting)
      }
    }
  }

  observeEvent(input$uploadGS, {
    withBusyIndicatorServer("uploadGS", {
      if (input$geneSetSourceChoice == "gsGMTUpload") {
        if (is.null(input$geneSetGMT)) {
          shinyjs::show(id = "gsUploadError", anim = FALSE)
        } else if (!nzchar(input$gsCollectionNameGMT)){
          shinyjs::show(id = "gsUploadError", anim = FALSE)
        } else {
          shinyjs::hide(id = "gsUploadError", anim = FALSE)
          vals$original <- importGeneSetsFromGMT(vals$original, input$geneSetGMT$datapath, collectionName = input$gsCollectionNameGMT)
        }

      } else if (input$geneSetSourceChoice == "gsDBUpload") {
        if (is.null(input$geneSetDB)) {
          shinyjs::show(id = "gsUploadError", anim = FALSE)
        } else {
          shinyjs::hide(id = "gsUploadError", anim = FALSE)
          vals$original <- importGeneSetsFromMSigDB(vals$original, input$geneSetDB)
        }

      } else if (input$geneSetSourceChoice == "gsPasteUpload") {
        handleGSPasteOption()
      }
      newGSchoices <- sctkListGeneSetCollections(vals$original)
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

  #----#
  # QC #
  #----#
  # Hide and show parameters for QC functions
  shinyjs::onclick("QCMetrics", shinyjs::toggle(id = "QCMetricsParams",
                                                anim = FALSE), add = TRUE)
  shinyjs::onclick("decontX", shinyjs::toggle(id = "decontXParams",
                                              anim = FALSE), add = TRUE)
  shinyjs::onclick("doubletCells", shinyjs::toggle(id = "doubletCellsParams",
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

  qc_choice_list <- list("doubletCells", "cxds", "bcds",
                         "cxds_bcds_hybrid", "decontX", "QCMetrics", "scrublet", "doubletFinder")
  # holds all the input ids for the QC algorithm parameters by algorithm name
  qc_input_ids <- list(doubletCells = list(nNeighbors="DCnNeighbors", simDoublets="DCsimDoublets"),

                       cxds = list(ntop="CXntop", binThresh="CXbinThresh", verb="CXverb", retRes="CXretRes", estNdbl="CXestNdbl"),

                       bcds = list(ntop="BCntop", srat="BCsrat", verb="BCverb", retRes="BCretRes", nmax="BCnmax", varImp="BCvarImp", estNdbl="BCestNdbl"),

                       cxds_bcds_hybrid = list(cxdsArgs=list(ntop="CX2ntop", binThresh="CX2binThresh", retRes="CX2retRes"),
                                               bcdsArgs=list(ntop="BC2ntop", srat="BC2srat", retRes="BC2retRes", namx="BC2nmax", varImp="BC2varImp"),
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
  qc_algo_status = reactiveValues(doubletCells=NULL, cxds=NULL, bcds=NULL, cxds_bcds_hybrid=NULL, decontX=NULL,
                                  QCMetrics=NULL, scrublet=NULL, doubletFinder=NULL)

  qc_plot_ids = reactiveValues(doubletCells="DCplots", cxds="CXplots", bcds="BCplots", cxds_bcds_hybrid="CXBCplots", decontX="DXplots",
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
    showModal(doubletCellsHelpModal())
  })
  observeEvent(input$QCMhelp, {
    showModal(QCMHelpModal())
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
      if (input[[algo]]) {
        return(TRUE)
      }
    }
    return(FALSE)
  }

  observeEvent(input$runQC, withConsoleMsgRedirect({
    withBusyIndicatorServer("runQC", {
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
        if (qcSample == "None") {
          qcSample <- NULL
        }
        algoList = list()
        paramsList <- list()
        for (algo in qc_choice_list) {
          if (input[[algo]]) {
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
        vals$counts <- runCellQC(inSCE = vals$original,
                                 algorithms = algoList,
                                 sample = qcSample,
                                 useAssay = input$qcAssaySelect,
                                 paramsList = paramsList)
        redDimList <- strsplit(reducedDimNames(vals$counts), " ")
        # run getUMAP if no UMAP
        if (!("UMAP" %in% redDimList)) {
          vals$counts <- getUMAP(inSCE = vals$counts,
                                 sample = qcSample,
                                 useAssay = input$qcAssaySelect,
          )
        }
        updateSelectInput(session, "qcPlotRedDim", choices = c(redDimList, "UMAP"))
        shinyjs::show(id = "qcPlotSection", anim = FALSE)
      }
    })
  }))

  observeEvent(input$qcPlotRedDim, {
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
        if (input[[algo]]) {
          algoList <- c(algoList, algo)
        }
      }
      # only run getUMAP if there are no reducedDimNames
      redDimName <- input$qcPlotRedDim
      # show the tabs for the result plots  output[[qc_plot_ids[[a]]]]
      showQCResTabs(vals, algoList, qc_algo_status, qc_plot_ids)
      arrangeQCPlots(vals$counts, output, algoList, colData(vals$counts)[,"sample"], qc_plot_ids, qc_algo_status, redDimName)
      uniqueSampleNames = unique(colData(vals$counts)[,"sample"])
      for (algo in algoList) {
        qc_algo_status[[algo]] <- list(self="done")
        if (length(uniqueSampleNames) > 1) {
          for (s in uniqueSampleNames) {
            qc_algo_status[[algo]][[s]] = TRUE
          }
        }
      }
    }
  })

  #-----------#
  # FILTERING #
  #-----------#

  filteringParams <- reactiveValues(params = list(), id_count = 0)
  rowFilteringParams <- reactiveValues(params = list(), id_count = 0)

  observeEvent(input$addFilteringParam, {
    showModal(filteringModal(colNames = names(colData(vals$counts))))
  })

  observeEvent(input$addRowFilteringParam, {
    if (!is.null(names(assays(vals$counts)))) {
      showModal(rowFilteringModal(assayInput = names(assays(vals$counts))))
    }
  })

  observeEvent(input$filterColSelect, {
    removeUI(selector = "#newThresh")
    isNum <- is.numeric(vals$counts[[input$filterColSelect]][0])
    if (length(vals$counts[[input$filterColSelect]]) > 0) {
      if (isNum) {
        minCol <- min(vals$counts[[input$filterColSelect]])
        maxCol <- max(vals$counts[[input$filterColSelect]])
        label_str <- sprintf("Please pick a number between %.5f and %.5f as a filtering threshold", minCol, maxCol)
        insertUI(
          selector = "#filterCriteria",
          ui = tags$div(id="newThresh", numericInput("filterThresh", label_str, minCol, min = minCol, max = maxCol))
        )
      } else {
        insertUI(
          selector = "#filterCriteria",
          ui = tags$div(id="newThresh",
                        checkboxGroupInput("filterThresh", "Please select which columns to keep:",
                                           choiceNames = as.vector(unique(vals$counts[[input$filterColSelect]])),
                                           choiceValues = as.vector(unique(vals$counts[[input$filterColSelect]]))
                        ),
          )
        )
      }
    } else {
      insertUI(
        selector = "#filterCriteria",
        ui = tags$div(id="newThresh", tags$b("This column does not have any filtering criteria", style = "color: red;"))
      )
    }
  })


  observeEvent(input$filterAssaySelect, {
    removeUI(selector = "#newThresh")
    insertUI(
      selector = "#rowFilterCriteria",
      ui = tags$div(id="newThresh",
                    numericInput("filterThreshX", "Number of counts per cell", 0),
                    numericInput("filterThreshY", "Number of cells", 0),
      )
    )

  })

  observeEvent(input$filtModalOK, {
    if ((!nzchar(input$filterThresh)) || (is.null(input$filterColSelect))) {
      showModal(filteringModal(failed=TRUE, colNames = names(colData(vals$counts))))
    } else {
      id <- paste0("filteringParam", filteringParams$id_count)
      # new row in parameters table
      addToColFilterParams(input$filterColSelect, input$filterThresh, id, filteringParams)
      threshStr <- ""
      if (is.numeric(input$filterThresh)) {
        threshStr <- sprintf("> %.5f", input$filterThresh)
      } else {
        threshStr <- paste(input$filterThresh, collapse = ', ')
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

  formatFilteringCriteria <- function(paramsReactive) {
    criteria = list()
    for (entry in paramsReactive) {
      criteria <- c(criteria, entry$param)
    }
    return(criteria)
  }

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
        vals$counts <- subsetSCERows(vals$counts, rowData = rowInput, returnAsAltExp = FALSE)
      }
    })
  })

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
  #  observeEvent(input$downsampleGo, {
  #    req(vals$counts)
  #    withBusyIndicatorServer("downsampleGo", {
  #      vals$counts <- vals$counts[, sample(ncol(vals$counts), input$downsampleNum)]
  #      updateNumSamples()
  #      vals$gsvaRes <- NULL
  #      vals$gsvaLimma <- NULL
  #      vals$visplotobject <- NULL
  #      vals$enrichRes <- NULL
  #      vals$dimRedPlot <- NULL
  #      vals$dimRedPlot_geneExp <- NULL
  #      vals$dendrogram <- NULL
  #      vals$pcX <- NULL
  #      vals$pcY <- NULL
  #    })
  #  })

  #Render summary table
  output$summarycontents <- renderTable({
    req(vals$counts)

    # Setting 'useAssay=NULL' assumes that the first assay is the one to count
    singleCellTK::summarizeSCE(inSCE = vals$counts,
                               useAssay = NULL,
                               sampleVariableName = "sample")
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
      updateNumSamples()
      # updateAssayInputs()
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
        if (input$assayModifyAction == "log") {
          if (input$trimAssayCheckbox) {
            assay(vals$counts, input$modifyAssayOutname) <- log2(assay(vals$counts, input$modifyAssaySelect) + 1)
            assay(vals$counts, input$modifyAssayOutname) <- trimCounts(assay(vals$counts, input$modifyAssayOutname), c(input$trimUpperValueAssay, input$trimLowerValueAssay))
          }
          else {
            assay(vals$counts, input$modifyAssayOutname) <- log2(assay(vals$counts, input$modifyAssaySelect) + 1)
          }
        }
        else if (input$assayModifyAction == "log1p") {
          if (input$trimAssayCheckbox) {
            assay(vals$counts, input$modifyAssayOutname) <- log1p(assay(vals$counts, input$modifyAssaySelect))
            assay(vals$counts, input$modifyAssayOutname) <- trimCounts(assay(vals$counts, input$modifyAssayOutname), c(input$trimUpperValueAssay, input$trimLowerValueAssay))
          }
          else {
            assay(vals$counts, input$modifyAssayOutname) <- log1p(assay(vals$counts, input$modifyAssaySelect))
          }
        }
        else if (input$assayModifyAction == "z.score") {
          if (input$trimAssayCheckbox) {
            assay(vals$counts, input$modifyAssayOutname) <- computeZScore(assay(vals$counts, input$modifyAssaySelect))
            assay(vals$counts, input$modifyAssayOutname) <- trimCounts(assay(vals$counts, input$modifyAssayOutname), c(input$trimUpperValueAssay, input$trimLowerValueAssay))
          }
          else {
            assay(vals$counts, input$modifyAssayOutname) <- computeZScore(assay(vals$counts, input$modifyAssaySelect))
          }
        }
        else if(input$assayModifyAction == "trim"){
          assay(vals$counts, input$modifyAssayOutname) <- trimCounts(assay(vals$counts, input$modifyAssaySelect), c(input$trimUpperValueAssay, input$trimLowerValueAssay))
        }
        else {
          showNotification("Error during assay transformation!", type = "error")
        }
        # updateAssayInputs()
      }
    })
  })

  observeEvent(input$assayModifyAction,{
    if (input$assayModifyAction == "log"){
      updateTextInput(session = session, inputId = "modifyAssayOutname", value = paste0(input$modifyAssaySelect, "Log"))
    }
    else if (input$assayModifyAction == "log1p"){
      updateTextInput(session = session, inputId = "modifyAssayOutname", value = paste0(input$modifyAssaySelect, "Log1p"))
    }
    else if (input$assayModifyAction == "z.score") {
      updateTextInput(session = session, inputId = "modifyAssayOutname", value = paste0(input$modifyAssaySelect, "Scaled"))
    }
    else if(input$assayModifyAction == "trim"){
      updateTextInput(session = session, inputId = "modifyAssayOutname", value = paste0(input$modifyAssaySelect, "Trim"))
    }

  })

  observeEvent(input$modifyAssaySelect,{
    if (input$assayModifyAction == "log"){
      updateTextInput(session = session, inputId = "modifyAssayOutname", value = paste0(input$modifyAssaySelect, "Log"))
    }
    else if (input$assayModifyAction == "log1p"){
      updateTextInput(session = session, inputId = "modifyAssayOutname", value = paste0(input$modifyAssaySelect, "Log1p"))
    }
    else if (input$assayModifyAction == "z.score") {
      updateTextInput(session = session, inputId = "modifyAssayOutname", value = paste0(input$modifyAssaySelect, "Scaled"))
    }
    else if(input$assayModifyAction == "trim"){
      updateTextInput(session = session, inputId = "modifyAssayOutname", value = paste0(input$modifyAssaySelect, "Trim"))
    }
  })

  observeEvent(input$normalizeAssay, {
    req(vals$counts)
    withBusyIndicatorServer("normalizeAssay", {
      if(!(input$normalizeAssaySelect %in% names(assays(vals$counts)))){
        stop("Selected assay does not exist!")
      }
      else if(input$normalizeAssayOutname == ""){
        stop("Assay Name cannot be empty!")
      }
      else if(input$normalizeAssayOutname %in% names(assays(vals$counts))){
        stop("Your selected Assay Name already exists! Try another Assay Name!")
      }
      else if(input$normalizeAssaySelect == ""){
        stop("Please select an assay before proceeding with normalization!")
      }
      else if(is.na(as.numeric(input$normalizationScaleFactor))){
        stop("Scaling factor must be a numeric non-empty value!")
      }
      else{
        if (input$normalizeAssayMethodSelect == "LogNormalize"
            || input$normalizeAssayMethodSelect == "CLR"
            || input$normalizeAssayMethodSelect == "RC") {
          vals$counts <- seuratNormalizeData(inSCE = vals$counts,
                                             useAssay = input$normalizeAssaySelect,
                                             normAssayName = input$normalizeAssayOutname,
                                             normalizationMethod = input$normalizeAssayMethodSelect,
                                             scaleFactor = as.numeric(input$normalizationScaleFactor))
          # updateAssayInputs()
        }
        else if (input$normalizeAssayMethodSelect == "CPM") {
          assay(vals$counts, input$normalizeAssayOutname) <- scater::calculateCPM(
            x = assay(vals$counts, input$normalizeAssaySelect))
          # updateAssayInputs()
        }
        else if(input$normalizeAssayMethodSelect == "LNC"){
          vals$counts <- scater_logNormCounts(
            inSCE = vals$counts,
            logAssayName = input$normalizeAssayOutname,
            useAssay = input$normalizeAssaySelect
          )
          # updateAssayInputs()
        }
        else if(input$normalizeAssayMethodSelect == "SCT"){
          vals$counts <- seuratSCTransform(
            inSCE = vals$counts,
            normAssayName = input$normalizeAssayOutname,
            useAssay = input$normalizeAssaySelect
          )
          # updateAssayInputs()
        }
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
    } else if(input$normalizeAssayMethodSelect == "LNC"){
      updateTextInput(session = session, inputId = "normalizeAssayOutname", value = "ScaterLogNormCounts")
    } else if(input$normalizeAssayMethodSelect == "SCT"){
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
  # Page 3: DR & Clustering
  #-----------------------------------------------------------------------------

  #Sidebar buttons functionality - not an accordion
  shinyjs::onclick("c_button2", shinyjs::toggle(id = "c_collapse2",
                                                anim = TRUE), add = TRUE)
  shinyjs::onclick("c_button3", shinyjs::toggle(id = "c_collapse3",
                                                anim = TRUE), add = TRUE)
  shinyjs::addClass(id = "c_button1", class = "btn-block")
  shinyjs::addClass(id = "c_button2", class = "btn-block")
  shinyjs::addClass(id = "c_button3", class = "btn-block")
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

  observeEvent(input$dimRedAltExpSelect, {
    if (!is.null(vals$counts) &&
        !is.null(input$dimRedAltExpSelect)) {
      ae <- altExp(vals$counts, input$dimRedAltExpSelect)
      aeAssays <- assayNames(ae)
      output$dimRedAltExpAssayUI <- renderUI({
        selectInput("dimRedAltExpAssay", "Select the Assay in the subset",
                    aeAssays)
      })
    }
  })

  output$dimRedNameUI <- renderUI({
    if (input$dimRedAssayType == 1){
      defaultText <- paste(input$dimRedAssaySelect, input$dimRedPlotMethod,
                           sep = '_')
    } else if (input$dimRedAssayType == 2){
      defaultText <- paste(input$dimRedAltExpAssay, input$dimRedPlotMethod,
                           sep = '_')
    }
    textInput('dimRedNameInput', "reducedDim Name:", defaultText)
  })

  observeEvent(input$runDimred, {
    if (!is.null(vals$counts)){
      withBusyIndicatorServer("runDimred", {
        if (input$dimRedNameInput == ""){
          shinyalert::shinyalert("Error", "enter a reducedDim name", type = "error")
        } #check for named entered and if its a duplicate
        else if (!is.null(input$dimRedNameInput)){
          if (input$dimRedNameInput %in% names(reducedDims(vals$counts))){
            shinyalert(
              "Warning",
              "Name already exits. Overwrite?",
              "warning", showCancelButton = TRUE,
              confirmButtonText = "Overwrite",
              callbackR = function(x){if(isTRUE(x)){
                dimrednamesave <- gsub(" ", "_", input$dimRedNameInput)
                if (input$dimRedPlotMethod == "PCA"){
                  if (input$dimRedAssayType == 1) {
                    vals$counts <- getPCA(inSCE = vals$counts,
                                          useAssay = input$dimRedAssaySelect,
                                          reducedDimName = dimrednamesave)
                  } else if (input$dimRedAssayType == 2) {
                    vals$counts <- getPCA(inSCE = vals$counts,
                                          useAssay = input$dimRedAltExpAssay,
                                          useAltExp = input$dimRedAltExpSelect,
                                          reducedDimName = dimrednamesave)
                  }
                } else if (input$dimRedPlotMethod == "tSNE"){
                  if (input$dimRedAssayType == 1) {
                    vals$counts <- getTSNE(inSCE = vals$counts,
                                           useAssay = input$dimRedAssaySelect,
                                           reducedDimName = dimrednamesave,
                                           perplexity = input$perplexityTSNE,
                                           n_iterations = input$iterTSNE)
                  } else if (input$dimRedAssayType == 2) {
                    vals$counts <- getTSNE(inSCE = vals$counts,
                                           useAssay = input$dimRedAltExpAssay,
                                           useAltExp = input$dimRedAltExpSelect,
                                           reducedDimName = dimrednamesave,
                                           perplexity = input$perplexityTSNE,
                                           n_iterations = input$iterTSNE)
                  }
                } else {
                  if (is.na(input$alphaUMAP)) {
                    stop("Learning rate (alpha) must be a numeric non-empty value!")
                  }
                  if (input$dimRedAssayType == 1) {
                    vals$counts <- getUMAP(inSCE = vals$counts,
                                           useAssay = input$dimRedAssaySelect,
                                           reducedDimName = dimrednamesave,
                                           nNeighbors = input$neighborsUMAP,
                                           nIterations = input$iterUMAP,
                                           minDist = input$mindistUMAP,
                                           alpha = input$alphaUMAP)
                  } else if (input$dimRedAssayType == 2) {
                    vals$counts <- getUMAP(inSCE = vals$counts,
                                           useAssay = input$dimRedAltExpAssay,
                                           useAltExp = input$dimRedAltExpSelect,
                                           reducedDimName = dimrednamesave,
                                           nNeighbors = input$neighborsUMAP,
                                           nIterations = input$iterUMAP,
                                           minDist = input$mindistUMAP,
                                           alpha = input$alphaUMAP)
                  }
                }
                updateReddimInputs()
              }}
            )
          } else {
            dimrednamesave <- gsub(" ", "_", input$dimRedNameInput)
            if (input$dimRedPlotMethod == "PCA"){
              if (input$dimRedAssayType == 1) {
                vals$counts <- getPCA(inSCE = vals$counts,
                                      useAssay = input$dimRedAssaySelect,
                                      reducedDimName = dimrednamesave)
              } else if (input$dimRedAssayType == 2) {
                vals$counts <- getPCA(inSCE = vals$counts,
                                      useAssay = input$dimRedAltExpAssay,
                                      useAltExp = input$dimRedAltExpSelect,
                                      reducedDimName = dimrednamesave)
              }
            } else if (input$dimRedPlotMethod == "tSNE"){
              if (input$dimRedAssayType == 1) {
                vals$counts <- getTSNE(inSCE = vals$counts,
                                       useAssay = input$dimRedAssaySelect,
                                       reducedDimName = dimrednamesave,
                                       perplexity = input$perplexityTSNE,
                                       n_iterations = input$iterTSNE)
              } else if (input$dimRedAssayType == 2) {
                vals$counts <- getTSNE(inSCE = vals$counts,
                                       useAssay = input$dimRedAltExpAssay,
                                       useAltExp = input$dimRedAltExpSelect,
                                       reducedDimName = dimrednamesave,
                                       perplexity = input$perplexityTSNE,
                                       n_iterations = input$iterTSNE)
              }
            } else {
              if (is.na(input$alphaUMAP)) {
                stop("Learning rate (alpha) must be a numeric non-empty value!")
              }
              if (input$dimRedAssayType == 1) {
                vals$counts <- getUMAP(inSCE = vals$counts,
                                       useAssay = input$dimRedAssaySelect,
                                       reducedDimName = dimrednamesave,
                                       nNeighbors = input$neighborsUMAP,
                                       nIterations = input$iterUMAP,
                                       minDist = input$mindistUMAP,
                                       alpha = input$alphaUMAP)
              } else if (input$dimRedAssayType == 2) {
                vals$counts <- getUMAP(inSCE = vals$counts,
                                       useAssay = input$dimRedAltExpAssay,
                                       useAltExp = input$dimRedAltExpSelect,
                                       reducedDimName = dimrednamesave,
                                       nNeighbors = input$neighborsUMAP,
                                       nIterations = input$iterUMAP,
                                       minDist = input$mindistUMAP,
                                       alpha = input$alphaUMAP)
              }
            }
            updateReddimInputs()
          }
        }
      })
    }
  })

  output$usingReducedDims <- renderUI({
    req(vals$counts)
    selectInput("usingReducedDims", "Select Reduced Dimension Data:", names(reducedDims(vals$counts)))
  })

  output$dimRedAxisSettings <- renderUI({
    req(vals$counts)
    req(input$usingReducedDims)
    if (any(grepl("PC*", colnames(reducedDim(vals$counts, input$usingReducedDims))))){
      pcComponents <- colnames(reducedDim(vals$counts, input$usingReducedDims))
      pcComponentsSelectedX <- pcComponents[1]
      pcComponentsSelectedY <- pcComponents[2]
      tagList(
        checkboxInput("checkAxis", label = "Modify axes", value = FALSE),
        conditionalPanel(
          condition = "input.checkAxis == true",
          h4("Axis Settings"),
          selectInput("pcX", "X Axis:", pcComponents),
          selectInput("pcY", "Y Axis:", pcComponents, selected = pcComponentsSelectedY)
        )
      )
    }
  })

  observeEvent(input$cUpdatePlot, {
    if (is.null(vals$counts)){
      shinyalert::shinyalert("Error!", "Upload data first.", type = "error")
    }  else {
      withBusyIndicatorServer("cUpdatePlot", {
        if (input$colorBy != "Gene Expression") {
          if (input$axisNames == TRUE) {
            if (input$dimRedAxis1 == "" & input$dimRedAxis2 == "") {
              shinyalert::shinyalert("Error", text = "Enter axis names", type = "error")
            } else {
              comp1 <- input$dimRedAxis1
              comp2 <- input$dimRedAxis2
            }
          } else {
            comp1 <- NULL
            comp2 <- NULL
          }
          #shinyjs doesn't have any visibility functions so have used the following conditions
          if (any(grepl("PC*", colnames(reducedDim(vals$counts, input$usingReducedDims))))){
            vals$pcX <- input$pcX
            vals$pcY <- input$pcY
          } else {
            vals$pcX <- NULL
            vals$pcY <- NULL
          }
          vals$dimRedPlot <- singleCellTK::plotDimRed(inSCE = vals$counts,
                                                      colorBy = input$colorBy,
                                                      shape = input$shapeBy,
                                                      useAssay = input$dimRedAssaySelect,
                                                      reducedDimName = input$usingReducedDims,
                                                      comp1 = comp1,
                                                      comp2 = comp2,
                                                      pcX = vals$pcX,
                                                      pcY = vals$pcY
          )
        }
      })
    }
  })

  # This code is commented out b/c it was causing a major lag whenever anything else was being
  # updated. Maybe it needs to be changed to observeEvent? - Josh
  # observe({
  #   output$geneExpPlot <- renderPlot({
  #   if (input$colorGeneBy == "Manual Input") {
  #     if (is.null(input$colorGenes)){
  #       ggplot2::ggplot() + ggplot2::theme_bw() +
  #         ggplot2::theme(plot.background = ggplot2::element_rect(fill = "white")) +
  #         ggplot2::theme(panel.border = ggplot2::element_rect(colour = "white"))
  #     } else {
  #       if (input$axisNames == TRUE) {
  #         if (input$dimRedAxis1 == "" & input$dimRedAxis2 == "") {
  #           shinyalert::shinyalert("Error", text = "Enter axis names", type = "error")
  #         } else {
  #           comp1 <- input$dimRedAxis1
  #           comp2 <- input$dimRedAxis2
  #         }
  #       } else {
  #         comp1 <- NULL
  #         comp2 <- NULL
  #       }
  #       #shinyjs doesn't have any visibility functions so have used the following conditions
  #       if (any(grepl("PC*", colnames(reducedDim(vals$counts, input$usingReducedDims))))){
  #         vals$pcX <- input$pcX
  #         vals$pcY <- input$pcY
  #       } else {
  #         vals$pcX <- NULL
  #         vals$pcY <- NULL
  #       }
  #       vals$dimRedPlot_geneExp <- singleCellTK::plotBiomarker(inSCE = vals$counts,
  #                                                              gene = input$colorGenes,
  #                                                              binary = input$colorBinary,
  #                                                              shape = input$shapeBy,
  #                                                              useAssay = input$dimRedAssaySelect,
  #                                                              reducedDimName = input$usingReducedDims,
  #                                                              comp1 = comp1, comp2 = comp2,
  #                                                              x = vals$pcX, y = vals$pcY)
  #       vals$dimRedPlot_geneExp
  #     }
  #   }
  # })
  # })

  #-----------------------------------------------------------------------------
  # Page 3: Clustering ####
  #-----------------------------------------------------------------------------
  scranSNNMats <- reactiveValues()
  output$clustScranSNNMatUI <- renderUI({
    if(!is.null(vals$counts)){
      choices <- list()

      nAssay <- length(assayNames(vals$counts))
      assayId <- NULL
      if (nAssay >= 1) {
        assayValues <- list()
        for (i in seq(nAssay)) {
          assayValues[[assayNames(vals$counts)[i]]] <- i
          assayId <- c(assayId, i)
        }
        choices[["Assays (full expression matrix)"]] <- assayValues
      }

      nReddim <- length(reducedDimNames(vals$counts))
      reddimId <- NULL
      if (nReddim >= 1) {
        reddimValues <- list()
        for (i in seq(nReddim)) {
          reddimValues[[reducedDimNames(vals$counts)[i]]] <- nAssay + i
          reddimId <- c(reddimId, nAssay + i)
        }
        choices[["ReducedDim (dimension reduction)"]] <- reddimValues
      }

      nAltExp <- length(altExpNames(vals$counts))
      altExpId <- NULL
      if (nAltExp >= 1) {
        altExpValues <- list()
        for (i in seq(nAltExp)) {
          altExpValues[[altExpNames(vals$counts)[i]]] <- nAssay + nReddim + i
          altExpId <- c(altExpId, nAssay + nReddim + i)
        }
        choices[["AltExp (subset of expression matrix)"]] <- altExpValues
      }

      scranSNNMats$allChoices <- list(assay = assayId, reducedDim = reddimId,
                                      altExp = altExpId)
      scranSNNMats$allNames <- c(assayNames(vals$counts),
                                 reducedDimNames(vals$counts),
                                 altExpNames(vals$counts))

      selectInput("clustScranSNNMat", "Select Input Matrix:", choices)
    } else {
      selectInput("clustScranSNNMat", "Select Input Matrix:", NULL)
    }
  })

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

  observeEvent(input$clustScranSNNMat, {
    output$clustScranSNNAltExpAssayUI <- renderUI({
      if (input$clustScranSNNMat %in% scranSNNMats$allChoices$altExp) {
        altExpName <- scranSNNMats$allNames[as.integer(input$clustScranSNNMat)]
        choices <- assayNames(altExp(vals$counts, altExpName))
        selectInput("clustScranSNNAltExpAssay", "Select the Assay in AltExp:",
                    choices = choices)
      }
    })
  })

  clustResults <- reactiveValues(names = NULL)

  observeEvent(input$clustRun, {
    if (is.null(vals$counts)){
      shinyalert::shinyalert("Error!", "Upload data first.", type = "error")
      print(input$clustAlgo)
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
          matNum <- as.integer(input$clustScranSNNMat)
          for (i in seq_along(names(scranSNNMats$allChoices))) {
            range <- scranSNNMats$allChoices[[i]]
            if (matNum %in% range) {
              matType <- names(scranSNNMats$allChoices)[i]
              break
            }
          }

          if (matType == "assay") {
            params$useAssay = scranSNNMats$allNames[matNum]
            params$nComp = input$clustScranSNNd
          } else if (matType == "reducedDim") {
            params$useReducedDim = scranSNNMats$allNames[matNum]
            updateSelectInput(session, "clustVisReddim",
                              selected = scranSNNMats$allNames[matNum])
          } else if (matType == "altExp") {
            params$useAltExp = scranSNNMats$allNames[matNum]
            params$altExpAssay = input$clustScranSNNAltExpAssay
            params$nComp = input$clustScranSNNd
          }
          vals$counts <- do.call(runScranSNN, params)
        } else if (input$clustAlgo %in% seq(7, 9)) {
          # K-Means
          if(input$clustKMeansReddim == ""){
            stop("Must select a reducedDim! If none available, ",
                 "compute one in the Dimensionality Reduction tab.")
          }
          if(is.na(input$clustKMeansN)){
            stop("Number of clusters/centers must be ",
                 "a numeric non-empty value!")
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
          updateSelectInput(session, "clustVisReddim", input$clustKMeansReddim)
        } else if (input$clustAlgo %in% seq(10, 12)) {
          # Seurat
          if(input$clustSeuratReddim == ""){
            stop("Must select a reducedDim! If none available, ",
                 "compute one in the Dimensionality Reduction tab.")
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
          updateSelectInput(session, "clustVisReddim", input$clustSeuratReddim)
        }
        updateColDataNames()
        clustResults$names <- c(clustResults$names, saveClusterName)
        updateSelectInput(session, "clustVisRes", choices = clustResults$names)
      })
    }
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
        shinyalert::shinyalert(
          "Error!",
          "No reduction selected. Select one or run dimension reduction first",
          type = "error")
      }
      inSCE <- vals$counts
      reducedDimName <- input$clustVisReddim
      output$clustVisPlot <- renderPlot({
        if (!is.null(choice) && choice != "" &&
            !is.null(reducedDimName) && reducedDimName != "") {
          plotSCEDimReduceColData(inSCE = inSCE,
                                  colorBy = choice,
                                  conditionClass = "factor",
                                  reducedDimName = reducedDimName,
                                  labelClusters = TRUE,
                                  dim1 = 1, dim2 = 2,
                                  legendTitle = choice)
        }
      })
    }
  })

  #-----------------------------------------------------------------------------
  # Page 3.2: Celda
  #-----------------------------------------------------------------------------

  # onclick buttons
  shinyjs::onclick("celdaBasicSet",
                   shinyjs::toggle(id = "celdaCollapse1",
                                   anim = TRUE), add = TRUE)
  shinyjs::onclick("celdaAdvSet",
                   shinyjs::toggle(id = "celdaCollapse2",
                                   anim = TRUE), add = TRUE)
  shinyjs::addClass(id = "celdaBasicSet", class = "btn-block")
  shinyjs::addClass(id = "celdaAdvSet", class = "btn-block")

  shinyjs::onclick("celdaBasicSetGS",
                   shinyjs::toggle(id = "celdaCollapseGS1",
                                   anim = TRUE), add = TRUE)
  shinyjs::onclick("celdaAdvSetGS",
                   shinyjs::toggle(id = "celdaCollapseGS2",
                                   anim = TRUE), add = TRUE)
  shinyjs::addClass(id = "celdaBasicSetGS", class = "btn-block")
  shinyjs::addClass(id = "celdaAdvSetGS", class = "btn-block")

  # shinyjs::onclick("celdatSNESet",
  #   shinyjs::toggle(id = "celdaCollapsetSNE",
  #     anim = TRUE), add = TRUE)
  # shinyjs::addClass(id = "celdatSNESet", class = "btn-block")

  observeEvent(input$celdamodsplit, {
    removeTab(inputId = "celdaModsplitTabset", target = "Perplexity Plot")
    removeTab(inputId = "celdaModsplitTabset", target = "Perplexity Diff Plot")
    appendTab(inputId = "celdaModsplitTabset", tabPanel(title = "Perplexity Plot",
      panel(heading = "Perplexity Plot",
        plotlyOutput(outputId = "plot_modsplit_perp", height = 300)
      )
    ), select = TRUE)
    appendTab(inputId = "celdaModsplitTabset", tabPanel(title = "Perplexity Difference Plot",
      panel(heading = "Perplexity Diff Plot",
        plotlyOutput(outputId = "plot_modsplit_perpdiff", height = 300)
      )
    ))
    withProgress(message = "Clustering Features", max = 1, value = 1, {
        vals$counts <- selectFeatures(as.matrix(counts(vals$counts)))
        assay(altExp(vals$counts), "normalizedCounts") <- normalizeCounts(counts(altExp(vals$counts)))
        updateNumericInput(session, "celdaLselect", min = input$celdaLinit, max = input$celdaLmax, value = input$celdaLinit)
        vals$counts <- recursiveSplitModule(vals$counts, initialL = input$celdaLinit, maxL = input$celdaLmax)
        output$plot_modsplit_perp <- renderPlotly({plotGridSearchPerplexity(vals$counts)})
        output$plot_modsplit_perpdiff <- renderPlotly({plotGridSearchPerplexityDiff(vals$counts)})
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
    vals$counts <- subsetCeldaList(vals$counts, params = list(L = input$celdaLselect))
    showNotification("Number of Feature Modules Selected.")
    updateCollapse(session = session, "CeldaUI", style = list("Identify # of Feature Modules" = "danger"))
    shinyjs::enable(
      selector = "div[value='Identify # of Cell Clusters']")
  })

  observeEvent(input$celdacellsplit, {
    removeTab(inputId = "celdaCellsplitTabset", target = "Perplexity Plot")
    appendTab(inputId = "celdaCellsplitTabset", tabPanel(title = "Perplexity Plot",
      panel(heading = "Perplexity Plot",
        plotlyOutput(outputId = "plot_cellsplit_perp", height = 300)
      )

    ), select = TRUE)
    withProgress(message = "Clustering Cells", max = 1, value = 1, {
      temp_umap <- celdaUmap(vals$counts)
      vals$counts <- recursiveSplitCell(vals$counts, initialK = input$celdaKinit, maxK = input$celdaKmax,
                                        yInit = celdaModules(vals$counts))
      output$plot_cellsplit_perp <- renderPlotly({plotGridSearchPerplexity(vals$counts)})
    })
    for (i in runParams(vals$counts)$K) {
      removeTab(inputId = "celdaCellsplitTabset", target = sprintf("Cluster %s", i))
      appendTab(inputId = "celdaCellsplitTabset", tabPanel(title = sprintf("Cluster %s", i),
        panel(heading = sprintf("Cluster %s", i),
          plotlyOutput(outputId = sprintf("plot_K_umap_%s", i), height = 300)
          )
      ))
      withProgress(message = "Plotting Clusters", max = 1, value = 1, {
        temp_model <- subsetCeldaList(vals$counts, params = list(K = i))
        output[[sprintf("plot_K_umap_%s", i)]] <- renderPlotly({plotDimReduceCluster(temp_model, dim1= reducedDim(altExp(temp_umap), "celda_UMAP")[, 1],
          dim2 = reducedDim(altExp(temp_umap), "celda_UMAP")[, 2], labelClusters = TRUE)})
      })
    shinyjs::enable(
        selector = sprintf(".celda_cellsplit_plots a[data-value='Cluster %s']", i))
    }
    shinyjs::enable(
      selector = ".celda_cellsplit_plots a[data-value='Perplexity Plot']")
    shinyjs::show(selector = ".celda_cellsplit_plots")
    showNotification("Cell Clustering Complete.")
    updateNumericInput(session, "celdaKselect", min = input$celdaKinit, max = input$celdaKmax, value = input$celdaKinit)
    shinyjs::show(id = "celdaKselect")
    shinyjs::show(id = "celdaKbtn")
  })

  observeEvent(input$celdaKbtn, {
    vals$counts <- subsetCeldaList(vals$counts, params = list(K = input$celdaKselect))
    showNotification("Number of Cell Clusters Selected.")
    updateCollapse(session = session, "CeldaUI", style = list("Identify # of Cell Clusters" = "danger"))
    shinyjs::enable(
      selector = "div[value='Visualization']")
    updateNumericInput(session, "celdamodheatmapnum", max = input$celdaKselect, value = 1)
  })

  observeEvent(input$CeldaUmap, {
    withProgress(message = "Computing Umap", max = 1, value = 1, {
      vals$counts <- celdaUmap(vals$counts)
      output$celdaumapplot <- renderPlotly({plotDimReduceCluster(vals$counts, reducedDimName = "celda_UMAP", xlab = "UMAP_1",
        ylab = "UMAP_2", labelClusters = TRUE)})
    })
    showNotification("Umap complete.")
    shinyjs::enable("CeldaTsne")
  })

  observeEvent(input$CeldaTsne, {
    withProgress(message = "Computing Tsne", max = 1, value = 1, {
      vals$counts <- celdaTsne(vals$counts)
      output$celdatsneplot <- renderPlotly({plotDimReduceCluster(vals$counts, reducedDimName = "celda_tSNE", xlab = "tSNE_1",
        ylab = "tSNE_2", labelClusters = TRUE)})
    })
    showNotification("Tsne complete.")
  })

  observeEvent(input$celdaheatmapbtn, {
    removeTab(inputId = "celdaHeatmapTabset", target = "Heatmap")
    removeTab(inputId = "celdaHeatmapTabset", target = "Module Heatmap")
    appendTab(inputId = "celdaHeatmapTabset", tabPanel(title = "Heatmap",
      panel(heading = "Heatmap",
        plotOutput(outputId = "celdaheatmapplt", height = 300)
      )
    ), select = TRUE)
    withProgress(message = "Plotting Heatmap", max = 1, value = 1, {
      output$celdaheatmapplt <- renderPlot({plot(celdaHeatmap(vals$counts))})
    })
    if (input$heatmap_module){
      appendTab(inputId = "celdaHeatmapTabset", tabPanel(title = "Module Heatmap",
        panel(heading = "Module Heatmap",
          plotOutput(outputId = "celdamodheatmapplt", height = 300)
        )
      ))
      withProgress(message = "Plotting Module Heatmap", max = 1, value = 1, {
        output$celdamodheatmapplt <- renderPlot({moduleHeatmap(vals$counts, featureModule = input$celdamodheatmapnum)})
      })
    }
    shinyjs::enable(
      selector = ".celda_heatmap_plots a[data-value='Heatmap']")
    shinyjs::toggleState(
      selector = ".celda_heatmap_plots a[data-value='Module Heatmap']",
      condition = input$celdamodheatmap)
    shinyjs::show(selector = ".celda_heatmap_plots")
    showNotification("Module heatmap complete.")
  })

  observeEvent(input$celdaprobplotbtn, {
    withProgress(message = "Plotting Module Heatmap", max = 1, value = 1, {
      output$celdaprobmapplt <- renderPlot({celdaProbabilityMap(vals$counts)})
    })
    showNotification("Probability map complete.")
  })

  observe({
    if(!is.null(vals$counts)){
      #If data is uploaded in data tab, enable first tab i.e. Normalization tab in Seurat workflow
      shinyjs::enable(
        selector = "div[value='Identify # of Feature Modules']")
    }else{
      #If no data uploaded in data tab, disabled all tabs and plots.

      #Disable tabs
      shinyjs::disable(
        selector = "div[value='Identify # of Feature Modules']")
      shinyjs::disable(
        selector = "div[value='Identify # of Cell Clusters']")
      shinyjs::disable(
        selector = "div[value='Visualization']")

      #Disable plots inside Modsplit subtab
      shinyjs::disable(
        selector = ".celda_modsplit_plots a[data-value='Perplexity Plot']")
      shinyjs::disable(
        selector = ".celda_modsplit_plots a[data-value='Perplexity Diff Plot']")

      #Disable plots inside Cellsplit subtab
      shinyjs::disable(
        selector = ".celda_cellsplit_plots a[data-value='Perplexity Plot']")

      #Disable plots inside Heatmap subtab
      shinyjs::disable(
        selector = ".celda_heatmap_plots a[data-value='Heatmap']")
      shinyjs::disable(
        selector = ".celda_heatmap_plots a[data-value='Module Heatmap']")
    }
  })

  # celda clustering tab
  observeEvent(input$runCelda2, {
    # is there an error or not
    if (is.null(cellsplit)) {
      shinyalert::shinyalert("Error!", "Upload data first.", type = "error")
    } else {
      # selected count matrix
      cm <- assay(cellsplit(), input$celdaAssay)
    }

    # And each row/column of the count matrix must have at least one count
    if (sum(rowSums(cm) == 0) >= 1 | sum(colSums(cm) == 0) >= 1) {
      shinyalert::shinyalert("Error!",
                             "Each row and column of the count matrix must have at least one count.
        Filter the data first.",
                             type = "error")
    }

    # Ensure that number of genes / cells is never more than
    # the number of requested clusters for each
    if (!is.null(input$cellClusterC) && ncol(cm) < input$cellClusterC) {
      shinyalert::shinyalert("Error!",
                             "Number of cells (columns) in count matrix must be >= K",
                             type = "error")
    }

    if (!is.null(input$geneModuleG) && nrow(cm) < input$geneModuleG) {
      shinyalert::shinyalert("Error!",
                             "Number of genes (rows) in count matrix must be >= L",
                             type = "error")
    }

    if (!is.null(input$geneModuleCG) && nrow(cm) < input$geneModuleCG) {
      shinyalert::shinyalert("Error!",
                             "Number of genes (rows) in count matrix must be >= L",
                             type = "error")
    }

    if (!is.null(input$cellClusterCG) && ncol(cm) < input$cellClusterCG) {
      shinyalert::shinyalert("Error!",
                             "Number of cells (columns) in count matrix must be >= K",
                             type = "error")
    } else {

      withBusyIndicatorServer("runCelda", {
        if (input$celdaModel == "celda_C") {
          vals$celdaMod <- celda_C(counts = assay(vals$counts,
                                                  input$celdaAssay),
                                   K = input$cellClusterC,
                                   alpha = input$celdaAlpha,
                                   beta = input$celdaBeta,
                                   algorithm = input$celdaAlgorithm,
                                   stopIter = input$celdaStopIter,
                                   maxIter = input$celdaMaxIter,
                                   splitOnIter = input$celdaSplitIter,
                                   seed = input$celdaSeed,
                                   nchains = input$celdaNChains)
          #cores = input$celdaCores)
          colData(vals$counts)$celdaCellCluster <- vals$celdaMod@clusters$z
          updateColDataNames()

        } else if (input$celdaModel == "celda_G") {
          vals$celdaMod <- celda_G(counts = assay(vals$counts,
                                                  input$celdaAssay),
                                   L = input$geneModuleG,
                                   beta = input$celdaBeta,
                                   delta = input$celdaDelta,
                                   gamma = input$celdaGamma,
                                   stopIter = input$celdaStopIter,
                                   maxIter = input$celdaMaxIter,
                                   splitOnIter = input$celdaSplitIter,
                                   seed = input$celdaSeed,
                                   nchains = input$celdaNChains)
          #cores = input$celdaCores)
          rowData(vals$counts)$celdaGeneModule <- vals$celdaMod@clusters$y
          updateFeatureAnnots()
          # update feature modules in module heatmap panel
          updateSelectInput(session, "celdaFeatureModule",
                            choices = seq_len(vals$celdaMod@params$L))

        } else if (input$celdaModel == "celda_CG") {
          vals$celdaMod <- celda_CG(counts = assay(vals$counts,input$celdaAssay),
                                    K = input$cellClusterCG,
                                    L = input$geneModuleCG,
                                    alpha = input$celdaAlpha,
                                    beta = input$celdaBeta,
                                    delta = input$celdaDelta,
                                    gamma = input$celdaGamma,
                                    algorithm = input$celdaAlgorithm,
                                    stopIter = input$celdaStopIter,
                                    maxIter = input$celdaMaxIter,
                                    splitOnIter = input$celdaSplitIter,
                                    seed = input$celdaSeed,
                                    nchains = input$celdaNChains)
          #cores = input$celdaCores)
          colData(vals$counts)$celdaCellCluster <- vals$celdaMod@clusters$z
          rowData(vals$counts)$celdaGeneModule <- rep(NA, nrow(vals$counts))
          rowData(vals$counts)$celdaGeneModule[ix] <- vals$celdaMod@clusters$y
          updateColDataNames()
          updateFeatureAnnots()
          # update feature modules in module heatmap panel
          updateSelectInput(session, "celdaFeatureModule",
                            choices = seq_len(vals$celdaMod@params$L))
          updateSelectInput(session, "celdatSNEFeature",
                            choices = vals$celdaMod@names$row)
        }
      })
    }
  })

  # disable the renderCeldaHeatmap button if no celda model is present
  isCeldaModelPresent <- reactive(is.null(vals$celdaMod))
  observe({
    if (isCeldaModelPresent()) {
      shinyjs::disable("renderHeatmap")
    } else {
      shinyjs::enable("renderHeatmap")
    }
  })


  # show celda heatmap
  observeEvent(input$renderHeatmap, {
    #celdaHeatmap <- eventReactive(input$showHeatmap, {

    # is there an error or not
    if (is.null(vals$celdaMod)) {
      shinyalert::shinyalert("Error!",
                             "Celda Model Not Found.",
                             type = "error")
    } else {

      withBusyIndicatorServer("renderHeatmap",
                              output$celdaHeatmap <- renderPlot({
                                g <- celdaHeatmap(counts = assay(vals$counts,input$celdaAssay),
                                                  celdaMod = vals$celdaMod)
                                g
                              }, height = 600)
      )}
  })


  #disable the downloadSCECelda button if no object is loaded
  isAssayResultCelda <- reactive(is.null(vals$counts))
  observe({
    if (isAssayResultCelda()) {
      shinyjs::disable("downloadSCECelda")
    } else {
      shinyjs::enable("downloadSCECelda")
    }
  })

  output$downloadSCECelda <- downloadHandler(
    filename <- function() {
      paste("SCE_", Sys.Date(), ".rds", sep = "")
    },
    content <- function(file) {
      saveRDS(vals$counts, file)
    })

  observeEvent(input$celdaSelectGSList, {
    updateSelectInput(session, "celdaSelectGSMod",
                      choices = vals$celdaListAllNames[[input$celdaSelectGSList]])
  })


  # show celda perplexity plot
  observeEvent(input$renderPerplexityPlot, {
    # is there an error or not
    if (is.null(vals$celdaList)) {
      shinyalert::shinyalert("Error!",
                             "Celda List Not Found.",
                             type = "error")
    } else {
      withBusyIndicatorServer("renderPerplexityPlot", {
        vals$celdaList <- resamplePerplexity(counts = assay(vals$counts,
                                                            input$celdaAssayGS),
                                             celdaList = vals$celdaList)
        output$celdaPerplexityPlot <- renderPlot({
          g <- plotGridSearchPerplexity(celdaList = vals$celdaList)
          g
        }, height = 600)
      })}
  })

  # Confirm celda model
  # disable the Confirm Selection button if no celda list is present
  isCeldaModelListPresent <- reactive(is.null(vals$celdaListAll))
  observe({
    if (isCeldaModelListPresent()) {
      shinyjs::disable("confirmCeldaModel")
    } else {
      shinyjs::enable("confirmCeldaModel")
    }
  })


  observeEvent(input$confirmCeldaModel, {
    withBusyIndicatorServer("confirmCeldaModel", {
      vals$celdaMod <- vals$celdaListAll[[
        input$celdaSelectGSList]][[input$celdaSelectGSMod]]
      # update data annotations
      if (paste(strsplit(input$celdaSelectGSMod,
                         "_")[[1]][seq(2)], collapse  = "_") == "celda_C") {
        colData(vals$counts)$celdaCellCluster <- vals$celdaMod@clusters$z
        updateColDataNames()
      } else if (paste(strsplit(input$celdaSelectGSMod,
                                "_")[[1]][seq(2)], collapse  = "_") == "celda_G") {
        rowData(vals$counts)$celdaGeneModule <- vals$celdaMod@clusters$y
        updateFeatureAnnots()
      } else if (paste(strsplit(input$celdaSelectGSMod,
                                "_")[[1]][seq(2)], collapse  = "_") == "celda_CG") {
        colData(vals$counts)$celdaCellCluster <- vals$celdaMod@clusters$z
        rowData(vals$counts)$celdaGeneModule <- vals$celdaMod@clusters$y
        updateColDataNames()
        updateFeatureAnnots()
        updateSelectInput(session, "celdatSNEFeature",
                          choices = vals$celdaMod@names$row)
      }
      # update feature modules in module heatmap panel
      updateSelectInput(session, "celdaFeatureModule",
                        choices = seq_len(vals$celdaMod@params$L))
    })
  })

  # Visualize tab
  # tSNE sub-panel
  # disable the runCeldatSNE button if no celda model is present
  observe({
    if (isCeldaModelPresent()) {
      shinyjs::disable("runCeldatSNE")
    } else {
      shinyjs::enable("runCeldatSNE")
    }
  })

  # run celda tSNE
  observeEvent(input$runCeldatSNE, {
    withBusyIndicatorServer("runCeldatSNE", {
      vals$celdatSNE <- celdaTsne(counts = assay(vals$counts,
                                                 input$celdaAssaytSNE),
                                  celdaMod = vals$celdaMod,
                                  maxCells = input$celdatSNEmaxCells,
                                  minClusterSize = input$celdatSNEminClusterSize,
                                  perplexity = input$celdatSNEPerplexity,
                                  maxIter = input$celdatSNEmaxIter,
                                  seed = input$celdatSNESeed)
    })
  })

  # disable the renderCeldatSNE button if no celda tSNE result is present
  isCeldatSNEResult <- reactive(is.null(vals$celdatSNE))
  observe({
    if (isCeldatSNEResult()) {
      shinyjs::disable("renderCeldatSNEByCellCluster")
      shinyjs::disable("renderCeldatSNEModule")
      shinyjs::disable("renderCeldatSNEFeature")
    } else {
      shinyjs::enable("renderCeldatSNEByCellCluster")
      shinyjs::enable("renderCeldatSNEModule")
      shinyjs::enable("renderCeldatSNEFeature")
    }
  })

  # show celda tSNE colored by cell cluster
  observeEvent(input$renderCeldatSNEByCellCluster, {
    withBusyIndicatorServer("renderCeldatSNEByCellCluster",
                            output$celdatSNECellClusterPlot <- renderPlot({
                              g <- plotDimReduceCluster(
                                dim1 = vals$celdatSNE[, 1],
                                dim2 = vals$celdatSNE[, 2],
                                cluster = clusters(vals$celdaMod)$z)
                              g}, height = 600)
    )}
  )

  # show celda tSNE colored by module probabilities
  observeEvent(input$renderCeldatSNEModule, {
    withBusyIndicatorServer("renderCeldatSNEModule",
                            output$celdatSNEModulePlot <- renderPlot({
                              g <- plotDimReduceModule(
                                dim1 = vals$celdatSNE[, 1],
                                dim2 = vals$celdatSNE[, 2],
                                counts = assay(vals$counts, input$celdaAssaytSNE),
                                celdaMod = vals$celdaMod)
                              g}, height = 600)
    )}
  )

  # show celda tSNE colored by feature expression
  observeEvent(input$renderCeldatSNEFeature, {
    withBusyIndicatorServer("renderCeldatSNEFeature",
                            output$celdatSNEFeaturePlot <- renderPlot({
                              g <- plotDimReduceFeature(
                                dim1 = vals$celdatSNE[, 1],
                                dim2 = vals$celdatSNE[, 2],
                                counts = assay(vals$counts, input$celdaAssaytSNE),
                                features = input$celdatSNEFeature)
                              g}, height = 600)
    )}
  )

  # Probability Map sub-panel
  # disable the renderCeldaProbabilityMap button if no celda model is present
  observe({
    if (isCeldaModelPresent()) {
      shinyjs::disable("renderCeldaProbabilityMap")
    } else {
      shinyjs::enable("renderCeldaProbabilityMap")
    }
  })

  # show celda Probability Map
  observeEvent(input$renderCeldaProbabilityMap, {
    withBusyIndicatorServer("renderCeldaProbabilityMap",
                            output$celdaProbabilityMapPlot <- renderPlot({
                              g <- celdaProbabilityMap(counts = assay(vals$counts,
                                                                      input$celdaAssayProbabilityMap),
                                                       celdaMod = vals$celdaMod)
                              g}, height = 600)
    )}
  )


  # Module Heatmap sub-panel
  # disable the renderCeldaModuleHeatmap button if no celda model is present
  observe({
    if (isCeldaModelPresent()) {
      shinyjs::disable("renderCeldaModuleHeatmap")
    } else {
      shinyjs::enable("renderCeldaModuleHeatmap")
    }
  })

  # show celda Module Heatmap
  observeEvent(input$renderCeldaModuleHeatmap, {
    withBusyIndicatorServer("renderCeldaModuleHeatmap",
                            output$celdaModuleHeatmapPlot <- renderPlot({
                              g <- moduleHeatmap(counts = assay(vals$counts,
                                                                input$celdaAssayModuleHeatmap),
                                                 celdaMod = vals$celdaMod,
                                                 featureModule = as.integer(input$celdaFeatureModule),
                                                 topCells = input$celdaModuleTopCells,
                                                 #normalize = as.logical(input$celdaFeatureModuleNormalize),
                                                 showFeaturenames = as.logical(input$celdaModuleFeatureNames))
                              g}, height = 600)
    )}
  )


  # download all celda lists
  # disable the downloadCeldaList button if no celda list is present
  isAssayResultCeldaList <- reactive(is.null(vals$celdaListAll))
  observe({
    if (isAssayResultCeldaList()) {
      shinyjs::disable("downloadAllCeldaLists")
    } else {
      shinyjs::enable("downloadAllCeldaLists")
    }
  })

  output$downloadAllCeldaLists <- downloadHandler(
    filename <- function() {
      paste("All_Celda_Lists_", Sys.Date(), ".rds", sep = "")
    },
    content <- function(file) {
      saveRDS(vals$celdaListAll, file)
    })


  #-----------------------------------------------------------------------------
  # Page 3.3: Cell Viewer
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
        cell_list <- BiocGenerics::colnames(vals$counts)
        gene_list <- BiocGenerics::rownames(vals$counts)
        method_list <- names(assays(vals$counts))
        approach_list <- names(reducedDims(vals$counts))
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
        updateSelectInput(session, "QuickAccess",
                          choices = c("",approach_list, "Custom"))
        updateSelectInput(session, "ApproachSelect_Xaxis",
                          choices = c(approach_list))
        updateSelectInput(session, "AdvancedMethodSelect_Xaxis",
                          choices = c(method_list))
        updateSelectInput(session, "GeneSelect_Assays_Xaxis",
                          choices = c(gene_list))
        updateSelectInput(session, "AnnotationSelect_Xaxis",
                          choices = c(annotation_list))
        updateSelectInput(session, "ApproachSelect_Yaxis",
                          choices = c(approach_list))
        updateSelectInput(session, "AdvancedMethodSelect_Yaxis",
                          choices = c(method_list))
        updateSelectInput(session, "GeneSelect_Assays_Yaxis",
                          choices = c(gene_list))
        updateSelectInput(session, "AnnotationSelect_Yaxis",
                          choices = c(annotation_list))
        updateSelectInput(session, "ApproachSelect_Colorby",
                          choices = c(approach_list))
        updateSelectInput(session, "AdvancedMethodSelect_Colorby",
                          choices = c(method_list))
        updateSelectInput(session, "GeneSelect_Assays_Colorby",
                          choices = c(gene_list))
        updateSelectInput(session, "AnnotationSelect_Colorby",
                          choices = c(annotation_list))
        updateSelectizeInput(session, "adjustgroupby", label = NULL, choices = c("None", annotation_list))
        updateSelectizeInput(session,"adjustbrewer", label = "Color Palettes:",
                             choices = c("RdYlBu",color_seqdiv))
      }
    }
  })

  hide_TypeSelect <- reactiveVal("hide")
  hide_bins <- reactiveVal()

  observeEvent(input$viewertabs, {
    if(!is.null(vals$counts)) {
      if(!is.null(reducedDims(vals$counts))) {
        approach_list <- names(reducedDims(vals$counts))
        if(input$viewertabs != "Scatter Plot"){
          updateSelectInput(session, "QuickAccess",
                            choices = c("Custom"))
          shinyjs::delay(5,shinyjs::disable("QuickAccess"))
        }else{
          updateSelectInput(session, "QuickAccess",
                            choices = c("", approach_list, "Custom"))
          shinyjs::delay(5,shinyjs::enable("QuickAccess"))
        }
        if(input$viewertabs == "Violin/Box Plot" || input$viewertabs == "Bar Plot"){
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
        }else{
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
        }
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
  cellviewer <- eventReactive(input$runCellViewer,{

    colors <- c()
    if (!is.null(numColors) && input$SelectColorType == 'Categorical') {
      for (i in 1: numColors) {
        colors[i] <- input[[ paste0(i,"_color")]]
      }
      names(colors) = colorLabels
    }
    #-+-+-+-+-+-cellviewer prepare3 : prepare Axis Label Name#####################
    ###Xaxis label name
    if(input$QuickAccess != "Custom" & input$QuickAccess != "" & input$adjustxlab == ""){
      xname <- paste0(input$QuickAccess, 1)
    }else if(input$QuickAccess != "Custom" & input$QuickAccess != ""& input$adjustxlab != ""){
      xname <- input$adjustxlab
    }else if(input$TypeSelect_Xaxis == 'Reduced Dimensions'){
      xname <- paste0(input$ApproachSelect_Xaxis,"_",substr(input$ColumnSelect_Xaxis,
                                                            str_length(input$ColumnSelect_Xaxis),str_length(input$ColumnSelect_Xaxis)))
    }else if(input$TypeSelect_Xaxis == 'Expression Assays'){
      xname <- input$GeneSelect_Assays_Xaxis
    }else{
      xname <- input$AnnotationSelect_Xaxis
    }
    ###Yaxis label name
    if(input$QuickAccess != "Custom" & input$QuickAccess != "" & input$adjustylab == ""){
      yname <- paste0(input$QuickAccess, 2)
    }else if(input$QuickAccess != "Custom" & input$QuickAccess != "" & input$adjustylab != ""){
      yname <- input$adjustylab
    }else if(input$TypeSelect_Yaxis == 'Reduced Dimensions'){
      yname <- paste0(input$ApproachSelect_Yaxis,"_",substr(input$ColumnSelect_Yaxis,
                                                            str_length(input$ColumnSelect_Yaxis),str_length(input$ColumnSelect_Yaxis)))
    }else if(input$TypeSelect_Yaxis == 'Expression Assays'){
      yname <- input$GeneSelect_Assays_Yaxis
    }else{
      yname <- input$AnnotationSelect_Yaxis
    }

    ###Yaxis label name
    if(input$TypeSelect_Colorby != 'Pick a Color'){
      if(input$TypeSelect_Colorby == 'Reduced Dimensions' && input$adjustlegendtitle == ""){
        legendname <- paste0(input$ApproachSelect_Colorby,"_",substr(input$ColumnSelect_Colorby,
                                                                     str_length(input$ColumnSelect_Colorby),str_length(input$ColumnSelect_Colorby)))
      }else if(input$TypeSelect_Colorby == 'Expression Assays' && input$adjustlegendtitle == ""){
        legendname <- input$GeneSelect_Assays_Colorby
      }else if(input$adjustlegendtitle == ""){
        legendname <- input$AnnotationSelect_Colorby
      }else{
        legendname <- input$adjustlegendtitle
      }
    }

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
                                      xlab = xname, ylab = yname, legendTitle = legendname, title = input$adjustitle,
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
        a <- plotSCEBarAssayData(vals$counts, title = input$adjusttitle,
          useAssay = input$AdvancedMethodSelect_Yaxis, groupBy = pltVars$groupby,
          feature = input$GeneSelect_Assays_Yaxis, transparency = input$adjustalpha,
          dotSize = input$adjustsize, combinePlot = "none", axisSize = input$adjustaxissize,
          axisLabelSize = input$adjustaxislabelsize, defaultTheme = as.logical(pltVars$defTheme))
      }else if(input$TypeSelect_Yaxis == "Cell Annotation"){
        a <- plotSCEBarColData(vals$counts, title = input$adjusttitle,
          coldata = input$AnnotationSelect_Yaxis, groupBy = pltVars$groupby,
          transparency = input$adjustalpha, dotSize = input$adjustsize, combinePlot = "none",
          axisSize = input$adjustaxissize, axisLabelSize = input$adjustaxislabelsize,
          defaultTheme = as.logical(pltVars$defTheme))
      }
    }else if(input$viewertabs == "Violin/Box Plot"){
      if(input$vlnboxcheck == TRUE){
        vln <- TRUE
        bx <- FALSE
      }else if(input$vlnboxcheck == FALSE){
        vln <- FALSE
        bx <- TRUE
      }
      if(input$TypeSelect_Yaxis == "Expression Assays"){
        a <- plotSCEViolinAssayData(vals$counts, violin = vln, box = bx,
          useAssay = input$AdvancedMethodSelect_Yaxis, title = input$adjusttitle,
          feature = input$GeneSelect_Assays_Yaxis, groupBy = pltVars$groupby,
          transparency = input$adjustalpha, dotSize = input$adjustsize, combinePlot = "none",
          axisSize = input$adjustaxissize, axisLabelSize = input$adjustaxislabelsize,
          defaultTheme = as.logical(pltVars$defTheme))
      }else if(input$TypeSelect_Yaxis == "Cell Annotation"){
        a <- plotSCEViolinColData(vals$counts, title = input$adjusttitle,
          coldata = input$AnnotationSelect_Yaxis, violin = vln, box = bx,
          groupBy = pltVars$groupby, transparency = input$adjustalpha,
          dotSize = input$adjustsize, combinePlot = "none", axisSize = input$adjustaxissize,
          axisLabelSize = input$adjustaxislabelsize, defaultTheme = as.logical(pltVars$defTheme))
      }
    }
    if (input$TypeSelect_Colorby == "Single Color"){
      a$layers[[1]]$aes_params$colour <- input$Col
    }
    if (input$adjustgridlines == TRUE){
      a <- a + ggplot2::theme_bw()
    }
    a <- plotly::ggplotly(a)
    plotly::subplot(plotlist = a, titleX = TRUE, titleY = TRUE)
  })
  output$scatter <- renderPlotly({cellviewer()})
  #
  #
  #-+-+-+-+-+-cellviewer prepare done: plot#####################
  ###plotly_after_reactive

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
    if(!is.null(vals$counts)){
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
              tmpSCE <- hmTemp$sce
              if(!is.null(hmTemp$colDataName)){
                cellAnnColor <- list()
                for(i in hmTemp$colDataName){
                  uniqs <- as.vector(unique(colData(tmpSCE)[[i]]))
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
                    if(is.numeric(colData(tmpSCE)[[i]])){
                      if(input[[paste0('hmcol', i, 'type')]] == 'Continuous'){
                        cFun <- circlize::colorRamp2(
                          c(min(colData(tmpSCE)[[i]]),
                            max(colData(tmpSCE)[[i]])),
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
                  uniqs <- as.vector(unique(rowData(tmpSCE)[[i]]))
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
                    if(is.numeric(rowData(tmpSCE)[[i]])){
                      if(input[[paste0('hmrow', i, 'type')]] == 'Continuous'){
                        cFun <- circlize::colorRamp2(
                          c(min(rowData(tmpSCE)[[i]]),
                            max(rowData(tmpSCE)[[i]])),
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
              if(is.null(hmTemp$rowSplitBy)){
                hmRowSplit <- NULL
              } else {
                hmRowSplit <- hmTemp$rowSplitBy
              }
              if(is.null(hmTemp$colSplitBy)){
                hmColSplit <- NULL
              } else {
                hmColSplit <- hmTemp$colSplitBy
              }
              trim <- input$hmTrim
              cs <- circlize::colorRamp2(
                c(trim[1], mean(trim), trim[2]),
                c(input$hmCSLow, input$hmCSMedium, input$hmCSHigh)
              )
              useAssay <- input$hmAssay
              cellIndex <- hmTemp$cellIndex
              featureIndex <- hmTemp$geneIndex
              rowDataName <- hmTemp$rowDataName
              colDataName <- hmTemp$colDataName
              scale <- input$hmScale
              output$Heatmap <- renderPlot({
                plotSCEHeatmap(
                  inSCE = tmpSCE, useAssay = useAssay, colorScheme = cs,
                  featureIndex = featureIndex, cellIndex = cellIndex,
                  rowDataName = rowDataName, colDataName = colDataName,
                  rowSplitBy = hmRowSplit, colSplitBy = hmColSplit,
                  rowLabel = hmAddLabel$gene, colLabel = hmAddLabel$cell,
                  rowDend = hmShowDendro[2], colDend = hmShowDendro[1],
                  scale = scale, trim = trim,
                  width = unit(20, 'cm'), height = unit(20, 'cm'),
                  featureAnnotationColor = geneAnnColor,
                  cellAnnotationColor = cellAnnColor
                )
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

  observeEvent(input$toggleNormalization, {
    if (vals$showAssayDetails == FALSE) {
      vals$showAssayDetails <- TRUE
      shinyjs::show(id="normalization", anim=TRUE, animType="slide", time=0.2)
      updateActionButton(session, "toggleAssayDetails", icon=icon("caret-up", lib="font-awesome"))
    } else {
      vals$showAssayDetails <- FALSE
      shinyjs::hide(id="normalization", anim=TRUE, animType="slide", time=0.2)
      updateActionButton(session, "toggleAssayDetails", icon=icon("caret-down", lib="font-awesome"))
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
        # Get "input" outside "renderPlot" expression so the plots aren't update
        # automatically but after pressing "Plot" button.
        ## Generals
        print(1)
        useAssay <- input$batchCheckOrigAssay
        batch <- input$batchCheckVar
        if(input$batchCheckCond == "None"){
          shapeBy <- NULL
        } else {
          shapeBy <- input$batchCheckCond
        }
        ## Original assay PCA
        oriAssayPCAName <- paste0(input$batchCheckOrigAssay, "_PCA")
        if(!oriAssayPCAName %in% names(reducedDims(vals$counts))){
          # TODO: Think about whether to perform this only on temp SCE
          # instead of vals$counts
          vals$counts <- getPCA(vals$counts,
                                useAssay = input$batchCheckOrigAssay,
                                reducedDimName = oriAssayPCAName)
          updateReddimInputs()
        }
        resName <- input$batchCheckCorrName
        ## Corrected assay/altExp PCA
        if (vals$batchRes[[resName]] == 'assay'){
          corrAssayPCAName = paste0(resName, "_PCA")
          vals$counts <- getPCA(vals$counts, useAssay = resName,
                                reducedDimName = corrAssayPCAName)
          updateReddimInputs()
        } else if (vals$batchRes[[resName]] == 'altExp'){
          ae <- altExp(vals$counts, resName)
          corrAltExpPCAName <- paste0(resName, "_PCA")
          ae <- getPCA(ae, useAssay = resName,
                       reducedDimName = corrAltExpPCAName)
          reducedDim(vals$counts, corrAltExpPCAName) <-
            reducedDim(ae, corrAltExpPCAName)
          updateReddimInputs()
        }
        inSCE <- vals$counts
        ## Update plots
        output$batchOriVars <- renderPlot({
          plotBatchVariance(inSCE = inSCE, useAssay = useAssay, batch = batch,
                            condition = shapeBy)
        })
        output$batchOriPCA <- renderPlot({
          plotSCEDimReduceColData(inSCE, colorBy = batch, shape = shapeBy,
                                  reducedDimName = oriAssayPCAName,
                                  title = paste0("Original ", useAssay, " PCA"))
        })
        output$batchCorrVars <- renderPlot({
          if (vals$batchRes[[resName]] == 'reddim'){
            plotBatchVariance(inSCE = inSCE, useReddim = resName, batch = batch,
                              condition = shapeBy)
          } else if (vals$batchRes[[resName]] == 'assay'){
            plotBatchVariance(inSCE = inSCE, useAssay = resName, batch = batch,
                              condition = shapeBy)
          } else if (vals$batchRes[[resName]] == 'altExp'){
            plotBatchVariance(inSCE = inSCE, useAltExp = resName, batch = batch,
                              condition = shapeBy)
          }
        })
        output$batchCorrReddim <- renderPlot({
          if (vals$batchRes[[resName]] == 'reddim'){
            plotSCEDimReduceColData(inSCE, colorBy = batch, shape = shapeBy,
                                    reducedDimName = resName,
                                    conditionClass = "character",
                                    title = paste0(resName, " corrected"))
          } else if (vals$batchRes[[resName]] == 'assay'){
            plotSCEDimReduceColData(inSCE, colorBy = batch, shape = shapeBy,
                                    reducedDimName = corrAssayPCAName,
                                    conditionClass = "character",
                                    title = paste0(resName, " corrected"))
          } else if (vals$batchRes[[resName]] == 'altExp'){
            plotSCEDimReduceColData(inSCE, colorBy = batch, shape = shapeBy,
                                    reducedDimName = corrAltExpPCAName,
                                    title = paste0(resName, " corrected"))
          }
        })
      })
    }
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

  output$selectCombatRefBatchUI <- renderUI({
    if (!is.null(vals$counts)){
      if (input$combatRef){
        selectInput("combatRefBatch", "Choose Reference Batch:",
                    unique(sort(colData(vals$counts)[, input$batchCorrVar])))
      }
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
          if(input$combatCond == "None"){
            cov <- NULL
          } else {
            cov <- input$combatCond
          }
          print(input$combatParametric)
          par.prior <- ifelse(input$combatParametric == "Parametric",
                              TRUE, FALSE)
          print(par.prior)
          vals$counts <- runComBat(inSCE = vals$counts,
                                   batch = input$batchCorrVar,
                                   useAssay = input$batchCorrAssay,
                                   par.prior = par.prior, covariates = cov,
                                   mean.only = input$combatMeanOnly,
                                   ref.batch = input$combatRefBatch,
                                   assayName = saveassayname)
          vals$batchRes[[saveassayname]] <- 'assay'
          updateAssayInputs()
          shinyalert::shinyalert('Success!', 'ComBat completed.',
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

  observeEvent(input$HarmonyRun, {
    if (is.null(vals$counts)){
      shinyalert::shinyalert("Error!", "Upload data first.", type = "error")
    } else {
      withBusyIndicatorServer("HarmonyRun", {
        saveassayname <- gsub(" ", "_", input$HarmonySaveReddim)
        if(isTRUE(input$HarmonyPcInput)){
          useAssay <- input$HarmonyReddim
        } else {
          useAssay <- input$batchCorrAssay
        }
        if(is.na(as.numeric(input$HarmonyTheta))){
          stop("Theta value must be numeric.")
        } else {
          theta <- as.numeric(input$HarmonyTheta)
        }
        vals$counts <- runHarmony(vals$counts, useAssay = useAssay,
                                  pcInput = input$HarmonyPcInput,
                                  batch = input$batchCorrVar,
                                  reducedDimName = saveassayname,
                                  nComponents = input$HarmonyNComp,
                                  theta = theta, nIter = input$HarmonyNIter)
        shinyalert::shinyalert('Success!', 'Harmony completed.',
                               type = 'success')
        vals$batchRes[[saveassayname]] <- 'reddim'
        updateReddimInputs()
      })
    }
  })

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
                     value = 0, min = 0, step = 1)

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
        updateAltExpInputs()
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
  # Page 4.1: Feature Selection
  #-----------------------------------------------------------------------------

  observeEvent(input$findHvgButtonFS, {
    withBusyIndicatorServer("findHvgButtonFS", {
      if (!is.null(vals$counts)) {
        if (input$hvgMethodFS == "vst"
            || input$hvgMethodFS == "mean.var.plot"
            || input$hvgMethodFS == "dispersion") {
          withProgress(
            message = "Finding highly variable genes",
            max = 1, value = 1, {
              tryCatch(vals$counts <- seuratFindHVG(
                                        inSCE = vals$counts,
                                        useAssay = input$assaySelectFS_Norm,
                                        hvgMethod = input$hvgMethodFS,
                                        hvgNumber = 100),
                       error = function(e) {
                         stop("HVG computation failed. ",
                              "Try re-computing with a normalized assay!")
                       }
              )
              # vals$counts <- seuratFindHVG(vals$counts,
              # useAssay = input$assaySelectFS_Norm,
              # seuratWorkflow$geneNamesSeurat,
              # input$hvgMethodFS,
              # as.numeric(input$hvgNoFeaturesFS))
            })
        } else if (input$hvgMethodFS == "modelGeneVar") {
          vals$counts <- scran_modelGeneVar(inSCE = vals$counts,
                                           assayName = input$assaySelectFS_Norm)
        }
        vals$hvgCalculated$status <- TRUE
        vals$hvgCalculated$method <- input$hvgMethodFS
      }
    })
  })

  observeEvent(input$showHVG, {
    withBusyIndicatorServer("showHVG", {
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
        if (input$hvgMethodFS == "vst") {
          x <- rowData(vals$counts)$seurat_variableFeatures_vst_mean
          y <- rowData(vals$counts)$seurat_variableFeatures_vst_varianceStandardized
          labeling <- "Standardized Variance"
        } else if (input$hvgMethodFS == "mean.var.plot") {
          x <- rowData(vals$counts)$seurat_variableFeatures_mvp_mean
          y <- rowData(vals$counts)$seurat_variableFeatures_mvp_dispersionScaled
          labeling <- "Dispersion"
        } else if (input$hvgMethodFS == "dispersion") {
          x <- rowData(vals$counts)$seurat_variableFeatures_dispersion_mean
          y <- rowData(vals$counts)$seurat_variableFeatures_dispersion_dispersionScaled
          labeling <- "Dispersion"
        } else if (input$hvgMethodFS == "modelGeneVar") {
          x <- rowData(vals$counts)$scran_modelGeneVar_mean
          y <- rowData(vals$counts)$scran_modelGeneVar_totalVariance
          labeling <- "Variance"
        }
        vals$vfplot <- ggplot() +
          geom_point(aes(x = x, y = y)) +
          geom_point(aes(x = subset(x, rownames(vals$counts) %in% HVGs),
                         y = subset(y, rownames(vals$counts) %in% HVGs)),
                     colour = "red") +
          geom_label(aes(x = subset(x, rownames(vals$counts) %in% HVGs),
                         y = subset(y, rownames(vals$counts) %in% HVGs),
                         label = subset(rownames(vals$counts),
                                        rownames(vals$counts) %in% HVGs)),
                     colour = "red",
                     size = 2) +
          labs(x = "Mean", y = labeling)
        output$plotFS <- renderPlot({
          if (!is.null(vals$vfplot)) {
            vals$vfplot
          }
        }, width = 400, height = 400)
        output$hvgOutputFS <- renderText({HVGs})
      } else {
        shinyalert::shinyalert(
          "Error",
          text = "Please compute the variance before the visualization!",
          type = "error"
        )
      }
    })
  })

  addAltExp <- function(inSCE, useAssay, geneSet, altExpName,
                        overwrite = FALSE){
    if ((!altExpName %in% altExpNames(inSCE)) ||
        isTRUE(overwrite)) {
      mat <- assay(inSCE[geneSet,], useAssay)
      assayList <- list()
      assayList[[paste0(altExpName, useAssay)]] <- mat
      ae <- SingleCellExperiment(assays = assayList)
      altExp(inSCE, altExpName) <- ae
    }
    return(inSCE)
  }

  observeEvent(input$hvgSubsetRun, {
    withBusyIndicatorServer("hvgSubsetRun", {
    if (isTRUE(vals$hvgCalculated$status) &&
        !is.null(vals$hvgCalculated$method)) {
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
                vals$counts <- addAltExp(vals$counts, input$assaySelectFS_Norm,
                                         HVGs, input$hvgAltExpName, x)
                updateAltExpInputs()
            }
          })
      } else {
          HVGs <- getTopHVG(inSCE = vals$counts,
                            method = input$hvgMethodFS,
                            n = input$hvgNumberSelect)
          #make sure no NA's are introduced in HVGs
          HVGs <- stats::na.omit(HVGs)
          vals$counts <- addAltExp(vals$counts, input$assaySelectFS_Norm, HVGs,
                                   input$hvgAltExpName)
          updateAltExpInputs()
      }
    } else {
      shinyalert::shinyalert(
        "Error",
        text = "Please compute the variance before the subsetting!",
        type = "error"
      )
    }
  })
  })
  #-----------------------------------------------------------------------------
  # Page 5.1: Differential Expression ####
  #-----------------------------------------------------------------------------
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
                                     log2fcThreshold = input$mastFCThresh,
                                     fdrThreshold = input$mastFDRThresh,
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
    })
  }

  observeEvent(input$runDE, {
    if (is.null(vals$counts)){
      shinyalert("Error!", "Upload data first.", type = "error")
    } else if(input$deAnalysisName == ""){
      shinyalert("Error!",
                 "Please enter differential expression analysis name.",
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
    }
  })

  # DE: Result visualize ####
  output$deResSelUI <- renderUI({
    if(!is.null(vals$counts)){
      res <- names(metadata(vals$counts)$diffExp)
      selectInput("deResSel", "Select Differential Expression Analysis", res)
    }
  })
  # Threshold adapting plot
  observeEvent(input$deAssay, {
    if(!is.null(vals$counts)){
      # MAST style sanity check for whether logged or not
      x <- assay(vals$counts, input$deAssay)
      if (!all(floor(x) == x, na.rm = TRUE) & max(x, na.rm = TRUE) <
          100) {
        output$deSanityWarnThresh <- renderText("")
      } else {
        output$deSanityWarnThresh <- renderText("Selected assay seems not logged (MAST style sanity check). Forcing to plot anyway. ")
      }
      output$deThreshplot <- renderPlot({
        plotMASTThresholdGenes(inSCE = vals$counts,
                               useAssay = input$deAssay,
                               check_sanity = FALSE)
      }, height = 800)
    }
  })

  # Data table
  output$deResult <- DT::renderDataTable({
    if(!is.null(input$deResSel)){
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

  observeEvent(input$dePlotVio, {
    if(!is.null(input$deResSel) &&
       !input$deResSel == ""){
      sce <- vals$counts
      useResult <- input$deResSel
      nrow <- input$deVioNRow
      ncol <- input$deVioNCol

      useAssay <- metadata(vals$counts)$diffExp[[useResult]]$useAssay
      if(input$deVioLabel == "Default ID"){
        labelBy = NULL
      } else {
        labelBy = input$deVioLabel
      }
      # MAST style sanity check for whether logged or not
      x <- assay(vals$counts, useAssay)
      if (!all(floor(x) == x, na.rm = TRUE) & max(x, na.rm = TRUE) <
          100) {
        output$deSanityWarnViolin <- renderText("")
      } else {
        output$deSanityWarnViolin <- renderText("Selected assay seems not logged (MAST style sanity check). Forcing to plot anyway. ")
      }
      output$deViolinPlot <- renderPlot({
        plotDEGViolin(inSCE = sce, useResult = useResult,
                      #threshP = input$deVioUseThresh,
                      nrow = nrow, ncol = ncol, labelBy = labelBy,
                      check_sanity = FALSE)
      })
    }
  })
  # Linear Regression Plot
  output$deRegTotalUI <- renderUI({
    topN <- input$deRegNRow * input$deRegNCol
    p(as.character(topN))
  })

  observeEvent(input$dePlotReg, {
    if(!is.null(input$deResSel) &&
       !input$deResSel == ""){
      sce <- vals$counts
      useResult <- input$deResSel
      nrow <- input$deRegNRow
      ncol <- input$deRegNCol
      useAssay <- metadata(vals$counts)$diffExp[[useResult]]$useAssay
      if(input$deRegLabel == "Default ID"){
        labelBy = NULL
      } else {
        labelBy = input$deRegLabel
      }
      # MAST style sanity check for whether logged or not
      x <- assay(vals$counts, useAssay)
      if (!all(floor(x) == x, na.rm = TRUE) & max(x, na.rm = TRUE) <
          100) {
        output$deSanityWarnReg <- renderText("")
      } else {
        output$deSanityWarnReg <- renderText("Selected assay seems not logged (MAST style sanity check). Forcing to plot anyway. ")
      }
      output$deRegPlot <- renderPlot({
        plotDEGRegression(inSCE = sce, useResult = useResult,
                          #threshP = input$deVioUseThresh,
                          nrow = nrow, ncol = ncol, labelBy = labelBy,
                          check_sanity = FALSE)
      })
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
    selectInput("deHMSplitRow", "Split columns by", multiple = TRUE,
                choices = c('regulation', otherAvail),
                selected = 'regulation')
  })

  observeEvent(input$dePlotHM, {
    if(!is.null(input$deResSel) &&
       !input$deResSel == ""){
      sce <- vals$counts
      useResult <- input$deResSel
      onlyPos <- input$deHMPosOnly
      log2fcThreshold <- input$deHMFC
      fdrThreshold <- input$deHMFDR
      rowDataName <- input$deHMrowData
      colDataName <- input$deHMcolData
      colSplitBy <- input$deHMSplitCol
      rowSplitBy <- input$deHMSplitRow
      output$deHeatmap <- renderPlot({
        plotDEGHeatmap(inSCE = sce, useResult = useResult, onlyPos = onlyPos,
                       log2fcThreshold = log2fcThreshold,
                       fdrThreshold = fdrThreshold, rowDataName = rowDataName,
                       colDataName = colDataName, colSplitBy = colSplitBy,
                       rowSplitBy = rowSplitBy)
      })
    }
  })


  #-----------------------------------------------------------------------------
  # Page 5.2: Find Marker ####
  #-----------------------------------------------------------------------------
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
                                         log2fcThreshold = input$fmLogFC,
                                         fdrThreshold = input$fmFDR)
        shinyalert::shinyalert("Success", "Find Marker completed.",
                               "success")
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
  }, filter = "top")

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
  observeEvent(input$fmUseTopN, {
    if (!isTRUE(input$fmUseTopN)) {
      shinyjs::disable("fmTopN")
    } else {
      shinyjs::enable("fmTopN")
    }
  })

  output$fmHMAssayUI <- renderUI({
    if(!is.null(vals$counts)){
      allAssay <- assayNames(vals$counts)
      selectInput('fmHMAssay', "Assay to plot", allAssay,
                  selected = input$fmAssay)
    }
  })

  observeEvent(input$plotFM, {
    if(!is.null(vals$counts) &&
       'findMarker' %in% names(metadata(vals$counts)) &&
       !is.null(input$fmHMAssay)){
      withBusyIndicatorServer("plotFM", {
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
      inSCE <- vals$counts
      useAssay <- input$fmHMAssay
      orderBy <- input$fmHMOrder
      log2fcThreshold <- input$fmHMFC
      fdrThreshold <- input$fmHMFDR
      decreasing <- input$fmHMdec
      rowDataName <- input$fmHMrowData
      colDataName <- input$fmHMcolData
      if(!isTRUE(input$fmUseTopN)) {
        topN <- NULL
      } else {
        topN <- input$fmTopN
      }
      # Take value before rendering plot, so that the plot doesnt auto re-render
      # while we tweak the parameter
      output$fmHeatmap <- renderPlot({
        plotMarkerDiffExp(inSCE = inSCE, useAssay = useAssay, orderBy = orderBy,
                          log2fcThreshold = log2fcThreshold, topN = topN,
                          fdrThreshold = fdrThreshold, decreasing = decreasing,
                          rowDataName = rowDataName, colDataName = colDataName)
      })
    })
    }
  })

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
      # updateAssayInputs()
      vals$counts <- singleCellTK:::.seuratInvalidate(inSCE = vals$counts)
    })
    updateCollapse(session = session, "SeuratUI", style = list("Normalize Data" = "danger"))
    shinyjs::enable(selector = "div[value='Scale Data']")
    showNotification("Normalization Complete")
  })

  #Perform scaling
  observeEvent(input$scale_button, {
    req(vals$counts)
    withProgress(message = "Scaling", max = 1, value = 1, {
      vals$counts <- seuratScaleData(inSCE = vals$counts,
                                     useAssay = "seuratNormData",
                                     scaledAssayName = "seuratScaledData",
                                     model = input$model.use,
                                     scale = input$do.scale,
                                     center = input$do.center,
                                     scaleMax = input$scale.max)
      # updateAssayInputs()
      vals$counts <- singleCellTK:::.seuratInvalidate(inSCE = vals$counts, scaleData = FALSE)
    })
    updateCollapse(session = session, "SeuratUI", style = list("Scale Data" = "danger"))
    shinyjs::enable(selector = "div[value='Highly Variable Genes']")
    showNotification("Scale Complete")
  })

  #Find HVG
  observeEvent(input$find_hvg_button, {
    req(vals$counts)
    withProgress(message = "Finding highly variable genes", max = 1, value = 1, {
      vals$counts <- seuratFindHVG(inSCE = vals$counts,
                                   useAssay = "seuratNormData",
                                   hvgMethod = input$hvg_method,
                                   hvgNumber = as.numeric(input$hvg_no_features))

      vals$counts <- singleCellTK:::.seuratInvalidate(inSCE = vals$counts, scaleData = FALSE, varFeatures = FALSE)
    })
    withProgress(message = "Plotting HVG", max = 1, value = 1, {
      output$plot_hvg <- renderPlotly({
        plotly::ggplotly(seuratPlotHVG(vals$counts))
      })
    })
    updateCollapse(session = session, "SeuratUI", style = list("Highly Variable Genes" = "danger"))
    shinyjs::enable(selector = "div[value='Dimensionality Reduction']")
    showNotification("Find HVG Complete")
  })

  #Display highly variable genes
  output$hvg_output <- renderText({
    if (!is.null(vals$counts)) {
      if (!is.null(vals$counts@metadata$seurat$obj)) {
        if (length(slot(vals$counts@metadata$seurat$obj, "assays")[["RNA"]]@var.features) > 0) {
          singleCellTK:::.seuratGetVariableFeatures(vals$counts, input$hvg_no_features_view)
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
        plotly::ggplotly(seuratReductionPlot(inSCE = vals$counts,
                                             useReduction = "pca",
                                             showLegend = FALSE))
      })
    })
    if (input$pca_compute_elbow) {
      appendTab(inputId = "seuratPCAPlotTabset", tabPanel(title = "Elbow Plot",
                                                          panel(heading = "Elbow Plot",
                                                                plotlyOutput(outputId = "plot_elbow_pca")
                                                          )
      ))

      withProgress(message = "Generating Elbow Plot", max = 1, value = 1, {
        updateNumericInput(session = session, inputId = "pca_significant_pc_counter", value = singleCellTK:::.computeSignificantPC(vals$counts))
        output$plot_elbow_pca <- renderPlotly({
          seuratElbowPlot(inSCE = vals$counts,
                          significantPC = singleCellTK:::.computeSignificantPC(vals$counts))
        })
        output$pca_significant_pc_output <- renderText({
          paste("<p>Number of significant components suggested by ElbowPlot: <span style='color:red'>", singleCellTK:::.computeSignificantPC(vals$counts)," </span> </p> <hr>")
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
          plotly::ggplotly(seuratJackStrawPlot(inSCE = vals$counts,
                                               dims = input$pca_no_components))
        })
      })
    }
    if (input$pca_compute_heatmap) {
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
                                                                      jqui_resizable(plotOutput(outputId = "plot_heatmap_pca"), options = list(maxWidth = 700))
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
          seuratHeatmapPlot(plotObject = vals$counts@metadata$seurat$heatmap_pca,
                            dims = input$pca_no_components,
                            ncol = 2,
                            labels = c("PC1", "PC2", "PC3", "PC4"))
        })
        updatePickerInput(session = session, inputId = "picker_dimheatmap_components_pca", choices = singleCellTK:::.getComponentNames(vals$counts@metadata$seurat$count_pc, "PC"))
      })
    }
    updateCollapse(session = session, "SeuratUI", style = list("Dimensionality Reduction" = "danger"))

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

    showNotification("PCA Complete")
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
        plotly::ggplotly(seuratReductionPlot(inSCE = vals$counts,
                                             useReduction = "ica",
                                             showLegend = FALSE))
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
                                                                      jqui_resizable(plotOutput(outputId = "plot_heatmap_ica"), options = list(maxWidth = 700))
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
          seuratHeatmapPlot(plotObject = vals$counts@metadata$seurat$heatmap_ica,
                            dims = input$ica_no_components,
                            ncol = 2,
                            labels = c("IC1", "IC2", "IC3", "IC4"))
        })
        updatePickerInput(session = session, inputId = "picker_dimheatmap_components_ica", choices = singleCellTK:::.getComponentNames(vals$counts@metadata$seurat$count_ic, "IC"))
      })
    }
    updateCollapse(session = session, "SeuratUI", style = list("Dimensionality Reduction" = "danger"))

    #Enable/Disable ICA plot panels not selected for computation (Heatmap)
    shinyjs::enable(
      selector = ".seurat_ica_plots a[data-value='ICA Plot']")

    shinyjs::toggleState(
      selector = ".seurat_ica_plots a[data-value='Heatmap Plot']",
      condition = input$ica_compute_heatmap)

    shinyjs::enable(
      selector = "div[value='tSNE/UMAP']")

    shinyjs::show(selector = ".seurat_ica_plots")

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
      updateCollapse(session = session, "SeuratUI", style = list("Clustering" = "danger"))
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
            plotly::ggplotly(seuratReductionPlot(inSCE = vals$counts,
                                                 useReduction = "pca",
                                                 showLegend = TRUE))
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
            plotly::ggplotly(seuratReductionPlot(inSCE = vals$counts,
                                                 useReduction = "ica",
                                                 showLegend = TRUE))
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
            plotly::ggplotly(seuratReductionPlot(inSCE = vals$counts,
                                                 useReduction = "tsne",
                                                 showLegend = TRUE))
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
            plotly::ggplotly(seuratReductionPlot(inSCE = vals$counts,
                                                 useReduction = "umap",
                                                 showLegend = TRUE))
          })
        })
        shinyjs::toggleState(
          selector = ".seurat_clustering_plots a[data-value='UMAP Plot']",
          condition = !is.null(slot(vals$counts@metadata$seurat$obj, "reductions")[["umap"]]))
      }

      shinyjs::show(selector = ".seurat_clustering_plots")
    }
    else{
      showNotification(paste0("'", input$reduction_clustering_method, "' reduction not found in input object"))
    }
  })

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
          plotly::ggplotly(seuratReductionPlot(inSCE = vals$counts,
                                               useReduction = "tsne",
                                               showLegend = FALSE))
        })
      })
      updateCollapse(session = session, "SeuratUI", style = list("tSNE/UMAP" = "danger"))
      shinyjs::enable(selector = "div[value='Clustering']")

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
          plotly::ggplotly(seuratReductionPlot(inSCE = vals$counts,
                                               useReduction = "umap",
                                               showLegend = FALSE))
        })
      })
      updateCollapse(session = session, "SeuratUI", style = list("tSNE/UMAP" = "danger"))
      shinyjs::enable(selector = "div[value='Clustering']")
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
        seuratHeatmapPlot(plotObject = vals$counts@metadata$seurat$heatmap_pca,
                          dims = length(input$picker_dimheatmap_components_pca),
                          ncol = input$slider_dimheatmap_pca,
                          labels = input$picker_dimheatmap_components_pca)
      })
    }
  })

  #Customize heatmap (ica) with selected options
  observeEvent(input$plot_heatmap_ica_button, {
    if (!is.null(input$picker_dimheatmap_components_ica)) {
      output$plot_heatmap_ica <- renderPlot({
        seuratHeatmapPlot(plotObject = vals$counts@metadata$seurat$heatmap_ica,
                          dims = length(input$picker_dimheatmap_components_ica),
                          ncol = input$slider_dimheatmap_ica,
                          labels = input$picker_dimheatmap_components_ica)
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
  # Page: Column Annotation (colData)
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
  output$outputColumnAnnotationTable_colData <- renderUI(
    {
      output$colOutTable <- DT::renderDataTable({ DT::datatable(vals$columnAnnotation, editable = 'cell', options = list(pageLength = 5)) })
      DT::dataTableOutput("colOutTable")
    })

  #create selectinput for selecting attribute with colnames from incoming dataset
  #create selectinput for selecting attribute value
  output$inputSelectAttribute_colData <- renderUI(
    {
      if(ncol(vals$columnAnnotation) > 0){
        selectInput("inputSelectAttribute_colData", label = "select attribute", choices = colnames(vals$columnAnnotation))
      }
    })
  output$inputSelectAttributeDelete_colData <- renderUI(
    {
      if(ncol(vals$columnAnnotation) > 0){
        selectInput("inputSelectAttributeDelete_colData", label = "select attribute to delete", choices = colnames(vals$columnAnnotation))
      }
    })

  #create selectinput for selecting column to delete
  output$inputSelectAttributeValue_colData <- renderUI(
    {
      if(ncol(vals$columnAnnotation) > 0){
        selectInput("inputSelectAttributeValue_colData", label = "select attribute value", choices = vals$columnAnnotation[, match(input$inputSelectAttribute_colData, colnames(vals$columnAnnotation))])
      }
    })

  #create selectinput for selecting merge_1 attribute
  #create selectinput for selecting merge_2 attribute
  output$inputSelectAttributeMerge1_colData <- renderUI(
    {
      if(ncol(vals$columnAnnotation) > 0){
        selectInput("inputSelectAttributeMerge1_colData", label = "select first column", choices = colnames(vals$columnAnnotation))
      }
    })
  output$inputSelectAttributeMerge2_colData <- renderUI(
    {
      if(ncol(vals$columnAnnotation) > 0){
        selectInput("inputSelectAttributeMerge2_colData", label = "select second column", choices = colnames(vals$columnAnnotation))
      }
    })

  #create selectinput for selecting fill_1 attribute
  #create selectinput for selecting fill_2 attribute
  output$inputSelectAttributeFill1_colData <- renderUI(
    {
      if(ncol(vals$columnAnnotation) > 0){
        selectInput("inputSelectAttributeFill1_colData", label = "select attribute column", choices = colnames(vals$columnAnnotation))
      }
    })
  output$inputSelectAttributeFill2_colData <- renderUI(
    {
      if(ncol(vals$columnAnnotation) > 0){
        selectInput("inputSelectAttributeFill2_colData", label = "select column to fill", choices = colnames(vals$columnAnnotation))
      }
    })

  #create selectinput for selecting attribute value for magic fill
  output$inputSelectAttributeFillvalue_colData <- renderUI(
    {
      if(ncol(vals$columnAnnotation) > 0){
        selectInput("inputSelectAttributeFillvalue_colData", label = "select attribute value", choices = vals$columnAnnotation[, match(input$inputSelectAttributeFill1_colData, colnames(vals$columnAnnotation))])
      }
    })

  #update criteria parameter text input when attribute value selectinput is changed
  observeEvent(input$inputSelectAttributeValue_colData,
               {
                 updateTextInput(session = session, "inputCriteria_colData",value = input$inputSelectAttributeValue_colData)
               })

  #create selectinput for selecting attribute for clean operation
  output$inputSelectAttributeClean_colData <- renderUI(
    {
      if(ncol(vals$columnAnnotation) > 0){
        selectInput("inputSelectAttributeClean_colData", label = "select attribute column", choices = colnames(vals$columnAnnotation))
      }
    })

  #confirm create bin button
  observeEvent(input$buttonConfirmBin_colData,
               {
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
  observeEvent(input$buttonConfirmMerge_colData,
               {
                 df <- vals$columnAnnotation
                 colname1 <- input$inputSelectAttributeMerge1_colData
                 colname2 <- input$inputSelectAttributeMerge2_colData
                 df <- unite_(df, col = colname1, c(colname1, colname2), sep = input$inputSelectSeparatorMerge_colData)

                 vals$columnAnnotation <- df

                 output$changesWarning_colData <- renderUI({
                   HTML("<h5><span style='color:red'> You have made changes to the Cell Annotation data. Select 'Save' to finalize these changes or 'Reset' to discard the changes.</span></h5></br>")
                 })
                 showNotification("You have made changes to the Cell Annotation data. Select 'Save' to finalize these changes or 'Reset' to discard the changes.", type = "error")
               })

  #fill column button
  observeEvent(input$buttonConfirmFill_colData,
               {
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
  observeEvent(input$buttonConfirmClean_colData,
               {
                 #get df from reactive input, backup column datatypes and convert factor to character
                 data <- singleCellTK:::.manageFactor(vals$columnAnnotation, operation = "backup")
                 df <- data$df

                 #perform operation
                 selected_attribute <- input$inputSelectAttributeClean_colData
                 selected_column_no <- match(selected_attribute, colnames(df))
                 selected_choice <- input$inputRemovalOperation_colData
                 selected_choice_no <- match(selected_choice, c("remove alphabets", "remove digits", "remove spaces", "remove symbols"))

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
  observeEvent(input$buttonConfirmEmptyColumnName_colData,
               {
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
    showNotification("Changes saved successfully.")
  })


  #-----------------------------------------------------------------------------
  # Page: Row Annotation (rowData)
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
  output$outputColumnAnnotationTable_rowData <- renderUI(
    {
      output$rowOutTable <- DT::renderDataTable({ DT::datatable(vals$rowAnnotation, editable = 'cell', options = list(pageLength = 5)) })
      DT::dataTableOutput("rowOutTable")
    })

  #create selectinput for selecting attribute with colnames from incoming dataset
  #create selectinput for selecting attribute value
  output$inputSelectAttribute_rowData <- renderUI(
    {
      if(ncol(vals$rowAnnotation) > 0){
        selectInput("inputSelectAttribute_rowData", label = "select attribute", choices = colnames(vals$rowAnnotation))
      }
    })
  output$inputSelectAttributeDelete_rowData <- renderUI(
    {
      if(ncol(vals$rowAnnotation) > 0){
        selectInput("inputSelectAttributeDelete_rowData", label = "select attribute to delete", choices = colnames(vals$rowAnnotation))
      }
    })

  #create selectinput for selecting column to delete
  output$inputSelectAttributeValue_rowData <- renderUI(
    {
      if(ncol(vals$rowAnnotation) > 0){
        selectInput("inputSelectAttributeValue_rowData", label = "select attribute value", choices = vals$rowAnnotation[, match(input$inputSelectAttribute_rowData, colnames(vals$rowAnnotation))])
      }
    })

  #create selectinput for selecting merge_1 attribute
  #create selectinput for selecting merge_2 attribute
  output$inputSelectAttributeMerge1_rowData <- renderUI(
    {
      if(ncol(vals$rowAnnotation) > 0){
        selectInput("inputSelectAttributeMerge1_rowData", label = "select first column", choices = colnames(vals$rowAnnotation))
      }
    })
  output$inputSelectAttributeMerge2_rowData <- renderUI(
    {
      if(ncol(vals$rowAnnotation) > 0){
        selectInput("inputSelectAttributeMerge2_rowData", label = "select second column", choices = colnames(vals$rowAnnotation))
      }
    })

  #create selectinput for selecting fill_1 attribute
  #create selectinput for selecting fill_2 attribute
  output$inputSelectAttributeFill1_rowData <- renderUI(
    {
      if(ncol(vals$rowAnnotation) > 0){
        selectInput("inputSelectAttributeFill1_rowData", label = "select attribute column", choices = colnames(vals$rowAnnotation))
      }
    })
  output$inputSelectAttributeFill2_rowData <- renderUI(
    {
      if(ncol(vals$rowAnnotation) > 0){
        selectInput("inputSelectAttributeFill2_rowData", label = "select column to fill", choices = colnames(vals$rowAnnotation))
      }
    })

  #create selectinput for selecting attribute value for magic fill
  output$inputSelectAttributeFillvalue_rowData <- renderUI(
    {
      if(ncol(vals$rowAnnotation) > 0){
        selectInput("inputSelectAttributeFillvalue_rowData", label = "select attribute value", choices = vals$rowAnnotation[, match(input$inputSelectAttributeFill1_rowData, colnames(vals$rowAnnotation))])
      }
    })

  #update criteria parameter text input when attribute value selectinput is changed
  observeEvent(input$inputSelectAttributeValue_rowData,
               {
                 updateTextInput(session = session, "inputCriteria_rowData",value = input$inputSelectAttributeValue_rowData)
               })

  #create selectinput for selecting attribute for clean operation
  output$inputSelectAttributeClean_rowData <- renderUI(
    {
      if(ncol(vals$rowAnnotation) > 0){
        selectInput("inputSelectAttributeClean_rowData", label = "select attribute column", choices = colnames(vals$rowAnnotation))
      }
    })

  #confirm create bin button
  observeEvent(input$buttonConfirmBin_rowData,
               {
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
  observeEvent(input$buttonConfirmMerge_rowData,
               {
                 df <- vals$rowAnnotation
                 colname1 <- input$inputSelectAttributeMerge1_rowData
                 colname2 <- input$inputSelectAttributeMerge2_rowData
                 df <- unite_(df, col = colname1, c(colname1, colname2), sep = input$inputSelectSeparatorMerge_rowData)

                 vals$rowAnnotation <- df

                 output$changesWarning_rowData <- renderUI({
                   HTML("<h5><span style='color:red'> You have made changes to the Cell Annotation data. Select 'Save' to finalize these changes or 'Reset' to discard the changes.</span></h5></br>")
                 })
                 showNotification("You have made changes to the Cell Annotation data. Select 'Save' to finalize these changes or 'Reset' to discard the changes.", type = "error")
               })

  #fill column button
  observeEvent(input$buttonConfirmFill_rowData,
               {
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
  observeEvent(input$buttonConfirmClean_rowData,
               {
                 #get df from reactive input, backup column datatypes and convert factor to character
                 data <- singleCellTK:::.manageFactor(vals$rowAnnotation, operation = "backup")
                 df <- data$df

                 #operations
                 selected_attribute <- input$inputSelectAttributeClean_rowData
                 selected_column_no <- match(selected_attribute, colnames(df))
                 selected_choice <- input$inputRemovalOperation_rowData
                 selected_choice_no <- match(selected_choice, c("remove alphabets", "remove digits", "remove spaces", "remove symbols"))

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
  observeEvent(input$buttonConfirmEmptyColumnName_rowData,
               {
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
    showNotification("Changes saved successfully.")
  })


  #-----------------------------------------------------------------------------
  # Page Download
  #-----------------------------------------------------------------------------

  path = '~'

  observeEvent(
    ignoreNULL = TRUE,
    eventExpr = {
      input$outputDirectory
    },
    handlerExpr = {
      if (input$outputDirectory > 0) {
        # condition prevents handler execution on initial app launch
        path <<- shinyDirectoryInput::choose.dir(default = shinyDirectoryInput::readDirectoryInput(session, 'outputDirectory'))
        shinyDirectoryInput::updateDirectoryInput(session, 'outputDirectory', value = path)
      }
    }
  )

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
        return
      }

      if (input$exportChoice == "rds") {
        filename = paste("SCE-", Sys.Date(), ".rds", sep = "")
        saveRDS(vals$counts, paste(path, "/", filename, sep = ""))
      } else if (input$exportChoice == "annData") {
        exportassay <- input$exportAssay
        compression <- input$compression
        compressionOpts = input$compressionOpts
        forceDense <- input$forceDense
        overwrite <- if(input$overwrite == 'True') TRUE else FALSE
        exportSCEtoAnnData(sce=vals$counts,
                           useAssay = exportassay,
                           outputDir=input$outputDirectory__chosen_dir,
                           prefix = paste("SCE-", Sys.Date(),sep = ""),
                           overwrite=overwrite,
                           compression = compression,
                           compressionOpts = compressionOpts,
                           forceDense = forceDense)
      } else if (input$exportChoice == "textfile") {
        overwrite <- if(input$overwrite == 'True') TRUE else FALSE
        gzipped <- if(input$gzip == 'True') TRUE else FALSE
        exportSCEtoFlatFile(sce = vals$counts,
                            outputDir=path,
                            overwrite=overwrite,
                            gzipped=gzipped,
                            sample = paste("SCE-", Sys.Date(),sep = ""))
      }
    })
  })
})

