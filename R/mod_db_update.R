#' db_update UI Function
#'
#' @description A shiny Module.
#'
#' @param id,input,output,session Internal parameters for {shiny}.
#'
#' @noRd 
#'
#' @importFrom shiny NS tagList 
library(shiny)
library(DBI)

# TODO: options to update should only appear when search returns items

mod_db_update_ui <- function(id){
  ns <- NS(id)
  tagList(
    fluidPage(
      titlePanel("SQLite Database Search and Update"),
      sidebarLayout(
        sidebarPanel(
          textInput(ns("search_text"), "Search for records:"),
          textInput(ns("update_id"), "ID to Update:"),
          textInput(ns("new_value"), "New Value:"),
          actionButton(ns("search_button"), "Search"),
          actionButton(ns("update_button"), "Update")
        ),
        mainPanel(
          tableOutput("table")
        )
      )
    )
  )
}
    
#' db_update Server Functions
#'
#' @noRd 
mod_db_update_server <- function(id){
  moduleServer( id, function(input, output, session){
    ns <- session$ns
    con <- dbConnect(RSQLite::SQLite(), "../SQL/starbase.sqlite")
    
    # Define a reactive dataset based on the user's search
    searched_data <- reactive({
      if (input$search_button > 0) {
        sql <- paste0("SELECT * FROM mytable WHERE column_name LIKE '%", input$search_text, "%'")
        dbGetQuery(con, sql)
      }
    })
    
    # Render the table based on the searched dataset
    output$table <- renderTable({
      searched_data()
    })
    
    # Update a record in the database
    observeEvent(input$update_button, {
      if (!is.null(input$update_id) && !is.null(input$new_value)) {
        sql <- paste0("UPDATE mytable SET column_name = '", input$new_value, "' WHERE id = ", input$update_id)
        dbExecute(con, sql)
        # Update the displayed table after the update
        output$table <- renderTable({
          dbGetQuery(con, "SELECT * FROM mytable")
        })
      }
    })
    
    # Close the database connection when the app is closed
    session$onSessionEnded(function() {
      dbDisconnect(con)
    })
  })
}