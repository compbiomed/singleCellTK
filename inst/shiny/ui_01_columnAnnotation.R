shinyPanelColumnAnnotation <- fluidPage(
  includeCSS('styles.CSS'),
              panel(heading = "Options for editing and importing Cell Annotation data",
                              tabsetPanel(
                                 tabPanel("Bins",
                                          panel(
                                             uiOutput("inputSelectAttribute_colData"),
                                             uiOutput("inputSelectAttributeValue_colData"),
                                             textInput("inputCriteria_colData", "criteria parameter"),
                                             selectInput("inputOperator_colData", "select comparison", choices = c("=",">","<",">=","<=")),                 
                                             textInput("inputBinName_colData", "bin name"),
                                             actionButton("buttonConfirmBin_colData", "Confirm Bin")
                                          )
                                 ),
                                 tabPanel("Merge Columns",
                                          panel(
                                             uiOutput("inputSelectAttributeMerge1_colData"),
                                             uiOutput("inputSelectAttributeMerge2_colData"),
                                             textInput("inputSelectSeparatorMerge_colData", "add separator between merged values", value = "_"),
                                             actionButton("buttonConfirmMerge_colData","Confirm Merge")
                                          )
                                 ),
                                 tabPanel("Magic Fill",
                                          panel(
                                             uiOutput("inputSelectAttributeFill1_colData"),
                                             uiOutput("inputSelectAttributeFill2_colData"),
                                             uiOutput("inputSelectAttributeFillvalue_colData"),
                                             textInput("inputReplaceText_colData", "new value"),
                                             actionButton("buttonConfirmFill_colData","Confirm Fill")
                                          )
                                 ),
                                 tabPanel("Clean",
                                          panel(
                                             uiOutput("inputSelectAttributeClean_colData"),
                                             selectInput("inputRemovalOperation_colData", "select removal criteria", choices = c("remove alphabets", "remove digits", "remove spaces","remove symbols")),
                                             actionButton("buttonConfirmClean_colData","Confirm Clean")
                                          )
                                 ),
                                 tabPanel("Add Column",
                                          panel(
                                             textInput("inputEmptyColumnName_colData", "enter new empty column name"),
                                             textInput("inputDefaultValueAddColumn_colData", "default value to fill"),
                                             actionButton("buttonConfirmEmptyColumnName_colData","Create Column")
                                          )
                                 ),
                                 tabPanel("Delete Column",
                                          panel(
                                             uiOutput("inputSelectAttributeDelete_colData"),
                                             actionButton("buttonConfirmDeleteColumn_colData", "Delete")
                                          )
                                          ),
                                 tabPanel("Import Cell Annotation from File",
                                          panel(
                                            radioButtons(
                                              inputId = "editorChoiceRadio_colData",
                                              label = "Select source for column annotation:",
                                              choices = c("Replace Cell Annotations" = "replace", "Add to existing Cell Annotations" = "concatenate"),
                                              selected = "concatenate"
                                            ),
                                            h6("You can either replace the existing colData or you can add/merge the new colData with the existing one."),
                                            conditionalPanel(
                                              condition = "input.editorChoiceRadio_colData == 'concatenate'",
                                              HTML("<h6><span style='color:red'> Warning:</span> Adding to existing colData will override the columns with same names! </h6>")
                                            ),
                                            fileInput('uploadFile_colData', 'Choose file to upload',
                                                      accept = c(
                                                        'text/csv',
                                                        'text/comma-separated-values',
                                                        '.csv'
                                                      )),
                                            actionButton(inputId = "importDataButton_colData", label = "Import")  
                                          )
                                          )
                              ),
                              
                              br(),
                              uiOutput("changesWarning_colData"),
                              fluidRow(column(6,
                                              panel("Save",
                                                    h6("Changes made to the annotation must be saved before they can be used in other modules of the toolkit:"),
                                                    actionButton("buttonSave_colData","Save",icon = icon("save"))
                                              )
                              ),
                              column(6,
                                     panel("Reset",
                                           h6("Reset annotation to point after changes were last saved:"),
                                           actionButton("buttonRestore_colData","Reset",icon = icon("refresh"))
                                     )
                                     
                              )
                              )
                              ),
              
              panel(heading = "Table of Cell Annotations",
                              uiOutput("outputColumnAnnotationTable_colData")
              )
)

