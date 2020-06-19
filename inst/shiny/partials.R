#---------------#
# IMPORT MODALS #
#---------------#

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