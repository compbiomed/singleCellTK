shinyPanelRowAnnotation <- fluidPage(
  includeCSS('styles.CSS'),
             panel(heading = "Options for editing and importing Feature Annotation data",
                             tabsetPanel(
                               tabPanel("Bins",
                                        panel(
                                          uiOutput("inputSelectAttribute_rowData"),
                                          uiOutput("inputSelectAttributeValue_rowData"),
                                          textInput("inputCriteria_rowData", "criteria parameter"),
                                          selectInput("inputOperator_rowData", "select comparison", choices = c("=",">","<",">=","<=")),                 
                                          textInput("inputBinName_rowData", "bin name"),
                                          actionButton("buttonConfirmBin_rowData", "Confirm Bin")
                                        )
                               ),
                               tabPanel("Merge Columns",
                                        panel(
                                          uiOutput("inputSelectAttributeMerge1_rowData"),
                                          uiOutput("inputSelectAttributeMerge2_rowData"),
                                          textInput("inputSelectSeparatorMerge_rowData", "add separator between merged values", value = "_"),
                                          actionButton("buttonConfirmMerge_rowData","Confirm Merge")
                                        )
                               ),
                               tabPanel("Magic Fill",
                                        panel(
                                          uiOutput("inputSelectAttributeFill1_rowData"),
                                          uiOutput("inputSelectAttributeFill2_rowData"),
                                          uiOutput("inputSelectAttributeFillvalue_rowData"),
                                          textInput("inputReplaceText_rowData", "new value"),
                                          actionButton("buttonConfirmFill_rowData","Confirm Fill")
                                        )
                               ),
                               tabPanel("Clean",
                                        panel(
                                          uiOutput("inputSelectAttributeClean_rowData"),
                                          selectInput("inputRemovalOperation_rowData", "select removal criteria", choices = c("remove alphabets", "remove digits", "remove spaces","remove symbols")),
                                          actionButton("buttonConfirmClean_rowData","Confirm Clean")
                                        )
                               ),
                               tabPanel("Add Column",
                                        panel(
                                          textInput("inputEmptyColumnName_rowData", "enter new empty column name"),
                                          textInput("inputDefaultValueAddColumn_rowData", "default value to fill"),
                                          actionButton("buttonConfirmEmptyColumnName_rowData","Create Column")
                                        )
                               ),
                               tabPanel("Delete Column",
                                        panel(
                                          uiOutput("inputSelectAttributeDelete_rowData"),
                                          actionButton("buttonConfirmDeleteColumn_rowData", "Delete")
                                        )
                               ),
                               tabPanel("Import Feature Annotation from File",
                                        panel(
                                          radioButtons(
                                            inputId = "editorChoiceRadio_rowData",
                                            label = "Select source for column annotation:",
                                            choices = c("Replace Feature Annotations" = "replace", "Add to existing Feature Annotations" = "concatenate"), #update this to include values
                                            selected = "concatenate"
                                          ),
                                          h6("You can either replace the existing rowData or you can add/merge the new rowData with the existing one."),
                                          conditionalPanel(
                                            condition = "input.editorChoiceRadio_rowData == 'concatenate'",
                                            HTML("<h6><span style='color:red'> Warning:</span> Adding to existing rowData will override the columns with same names! </h6>")
                                          ),
                                          fileInput('uploadFile_rowData', 'Choose file to upload',
                                                    accept = c(
                                                      'text/csv',
                                                      'text/comma-separated-values',
                                                      '.csv'
                                                    )),
                                          actionButton(inputId = "importDataButton_rowData", label = "Import")  
                                        )
                               )
                             ),
                             
                             br(),
                             uiOutput("changesWarning_rowData"),
                             br(),
                             fluidRow(column(6,
                                             panel("Save",
                                                   h6("Changes made to the annotation must be saved before they can be used in other modules of the toolkit:"),
                                                   actionButton("buttonSave_rowData","Save",icon = icon("save"))
                                             )
                             ),
                             column(6,
                                    panel("Reset",
                                          h6("Reset annotation to point after changes were last saved:"),
                                          actionButton("buttonRestore_rowData","Reset",icon = icon("refresh"))
                                    )
                                    
                             )
                             )
             ),
             
             panel(heading = "Table of Feature Annotations",
                             uiOutput("outputColumnAnnotationTable_rowData")
             )
)

