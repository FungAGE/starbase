#' submit UI Function
#'
#' @description A shiny Module.
#'
#' @param id,input,output,session Internal parameters for {shiny}.
#'
#' @noRd
#'
#' @import DT RSQLite pool shinyjs uuid dplyr
#'
#' @importFrom shiny NS tagList

# Create the SQL database and responses table
pool <- pool::dbPool(RSQLite::SQLite(), dbname = "SQL/starbase.sqlite")

#Label mandatory fields
labelMandatory <- function(label) {
  tagList(
    label,
    span("*", class = "mandatory_star")
  )
}

mod_submit_ui <- function(id) {
  ns <- NS(id)
  # ui
  fluidPage(
    shinyjs::useShinyjs(),
    fluidRow(
      actionButton(ns("add_button"), "Add", icon("plus")),
      actionButton(ns("edit_button"), "Edit", icon("edit")),
      actionButton(ns("copy_button"), "Copy", icon("copy")),
      actionButton(ns("delete_button"), "Delete", icon("trash-alt"))
    ),
    br(),
    fluidRow(DT::dataTableOutput(ns("ship"))
    )
  )
}

#' submit Server Functions
#'
#' @noRd
mod_submit_server <- function(id) {
  moduleServer(id, function(input, output, session) {
    ns <- session$ns
      #load ship_df and make reactive to inputs  
      ship_df <- reactive({
        
        #make reactive to
        input$submit
        input$submit_edit
        input$copy_button
        input$delete_button
        
        dbReadTable(pool, "ship")
        
      })  
      
      #List of mandatory fields for submission
      fieldsMandatory <- c("fna","gff3","genus","species","evidence","uploader")
      
      #define which input fields are mandatory 
      observe({
        
        mandatoryFilled <-
          vapply(fieldsMandatory,
                 function(x) {
                   !is.null(input[[x]]) && input[[x]] != ""
                 },
                 logical(1))
        mandatoryFilled <- all(mandatoryFilled)
        
        shinyjs::toggleState(id = "submit", condition = mandatoryFilled)
      })
      
      #Form for data entry
      entry_form <- function(button_id){
        
        showModal(
          modalDialog(
            div(id=("entry_form"),
                tags$head(tags$style(".modal-dialog{ width:400px}")),
                tags$head(tags$style(HTML(".shiny-split-layout > div {overflow: visible}"))),
                fluidPage(
                  fluidRow(
                      textInput("uploader",label = labelMandatory("Name of curator")),
                      textInput("evidence",label=labelMandatory("What evidence exists for ship annotation?")),
                      textInput("genus",label = labelMandatory("Enter genus name")),
                      textInput("species",label = labelMandatory("Enter species name")),
                      
                      # TODO: store file and put path in SQL table
                      fileInput("fna",label = labelMandatory("Upload ship sequence"),accept = c(".fa",".fna",".fasta")),
                      fileInput("gff3",label = labelMandatory("Upload ship annotations (GFF3 format)"),accept = c(".gff",".gff3")),
                    textAreaInput("comment", "Comment", placeholder = "", height = 100, width = "354px"),
                    helpText(labelMandatory(""), paste("Mandatory field.")),
                    actionButton(button_id, "Submit")
                  ),
                  easyClose = TRUE
                )
            )
          )
        )
      }
      
      #
      fieldsAll <- c("fna","gff3","genus","species","evidence","uploader","comment")
      
      #save form data into data_frame format
      formData <- reactive({
        
        formData <- data.frame(row_id = UUIDgenerate(),
                               fna = input$fna,
                               gff3 = input$gff3,
                               genus = input$genus,
                               species = input$species,
                               evidence = input$evidence,
                               uploader = input$uploader,
                               comment = input$comment,
                               stringsAsFactors = FALSE)
        return(formData)
        
      })
      
      #Add data
      appendData <- function(data){
        query <- sqlAppendTable(pool, "ship", data, row.names = FALSE)
        dbExecute(pool, query)
      }
      
      observeEvent(input$add_button, priority = 20,{
        
        entry_form("submit")
        
      })
      
      observeEvent(input$submit, priority = 20,{
        
        appendData(formData())
        shinyjs::reset("entry_form")
        removeModal()
        
      })
      
      #delete data
      deleteData <- reactive({
        
        SQL_df <- dbReadTable(pool, "ship")
        row_selection <- SQL_df[input$ship_rows_selected, "row_id"]
        
        query <- lapply(row_selection, function(nr){
          
          dbExecute(pool, sprintf('DELETE FROM "ship" WHERE "row_id" == ("%s")', nr))
        })
      })
      
      observeEvent(input$delete_button, priority = 20,{
        
        if(length(input$ship_rows_selected)>=1 ){
          deleteData()
        }
        
        showModal(
          
          if(length(input$ship_rows_selected) < 1 ){
            modalDialog(
              title = "Warning",
              paste("Please select row(s)." ),easyClose = TRUE
            )
          })
      })
      
      #copy data
      unique_id <- function(data){
        replicate(nrow(data), UUIDgenerate())
      }
      
      copyData <- reactive({
        
        SQL_df <- dbReadTable(pool, "ship")
        row_selection <- SQL_df[input$ship_rows_selected, "row_id"] 
        SQL_df <- SQL_df %>% filter(row_id %in% row_selection)
        SQL_df$row_id <- unique_id(SQL_df)
        
        query <- sqlAppendTable(pool, "ship", SQL_df, row.names = FALSE)
        dbExecute(pool, query)
        
      })
      
      observeEvent(input$copy_button, priority = 20,{
        
        if(length(input$ship_rows_selected)>=1 ){
          copyData()
        }
        
        showModal(
          
          if(length(input$ship_rows_selected) < 1 ){
            modalDialog(
              title = "Warning",
              paste("Please select row(s)." ),easyClose = TRUE
            )
          })
        
      })
      
      #edit data
      observeEvent(input$edit_button, priority = 20,{
        
        SQL_df <- dbReadTable(pool, "ship")
        
        showModal(
          if(length(input$ship_rows_selected) > 1 ){
            modalDialog(
              title = "Warning",
              paste("Please select only one row." ),easyClose = TRUE)
          } else if(length(input$ship_rows_selected) < 1){
            modalDialog(
              title = "Warning",
              paste("Please select a row." ),easyClose = TRUE)
          })  
        
        if(length(input$ship_rows_selected) == 1 ){
          
          entry_form("submit_edit")
          
          # TODO: update file path text in table, not input
          updateTextInput(session, "fna", value = SQL_df[input$ship_rows_selected,"fna"])
          updateTextInput(session, "gff3", value = SQL_df[input$ship_rows_selected,"gff3"])
          updateTextInput(session, "genus", value = SQL_df[input$ship_rows_selected,"genus"])
          updateTextInput(session, "species", value = SQL_df[input$ship_rows_selected,"species"])
          updateTextInput(session, "evidence", value = SQL_df[input$ship_rows_selected,"evidence"])
          updateTextInput(session, "uploader", value = SQL_df[input$ship_rows_selected,"uploader"])
          updateTextInput(session, "comment", value = SQL_df[input$ship_rows_selected,"comment"])
          
          
        }
        
      })
      
      observeEvent(input$submit_edit, priority = 20, {
        
        SQL_df <- dbReadTable(pool, "ship")
        row_selection <- SQL_df[input$ship_row_last_clicked, "row_id"] 
        dbExecute(pool, sprintf('UPDATE "ship" SET "fna" = ?, "gff3" = ?, "genus" = ?, "species" = ?, "evidence" = ?, "uploader" = ?, "comment" = ? WHERE "row_id" = ("%s")', row_selection), 
                  param = list(
                    input$fna,
                    input$gff3,
                    input$genus,
                    input$species,
                    input$evidence,
                    input$uploader,
                    input$comment))
        removeModal()
        
      })
      
      
      output$ship <- DT::renderDataTable({
        
        table <- ship_df() # %>% select(-row_id) 
        
        # TODO: add names/labels to table columns?
        table <- datatable(table, 
                           rownames = FALSE,
                           options = list(searching = FALSE, lengthChange = FALSE)
        )
        
      })
    
  })
}