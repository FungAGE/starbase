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

pool <- pool::dbPool(RSQLite::SQLite(), dbname = "SQL/starbase.sqlite")
con <- pool::poolCheckout(pool)

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
      textInput(ns("uploader"),label = labelMandatory("Name of curator")),
      textInput(ns("evidence"),label=labelMandatory("What evidence exists for ship annotation?")),
      textInput(ns("genus"),label = labelMandatory("Enter genus name")),
      textInput(ns("species"),label = labelMandatory("Enter species name")),
      
      # TODO: store file and put path in SQL table
      fileInput(ns("fna"),label = labelMandatory("Upload ship sequence"),accept = c(".fa",".fna",".fasta")),
      fileInput(ns("gff3"),label = labelMandatory("Upload ship annotations (GFF3 format)"),accept = c(".gff",".gff3")),
      textAreaInput(ns("comment"), "Comment", placeholder = "", height = 100, width = "354px"),
      helpText(labelMandatory(""), paste("Mandatory field.")),
      actionButton(ns("submit"), "Submit")
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
        dbReadTable(con, "joined_ships")
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
        query <- sqlAppendTable(con, "joined_ships", data, row.names = FALSE)
        dbExecute(con, query)
      }
      
      observeEvent(input$submit, priority = 20,{
        
        appendData(formData())
      })
  })
}