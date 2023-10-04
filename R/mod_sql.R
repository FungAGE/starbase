#' sql UI Function
#'
#' @description A shiny Module.
#'
#' @param id,input,output,session Internal parameters for {shiny}.
#'
#' @noRd 
#'
#' @importFrom shiny NS tagList 
mod_sql_ui <- function(id){
  ns <- NS(id)
  tagList(
 
  )
}
    
#' sql Server Functions
#'
#' @noRd 
mod_sql_server <- function(id){
  moduleServer( id, function(input, output, session){
    ns <- session$ns
 
  })
}
    
## To be copied in the UI
# mod_sql_ui("sql_1")
    
## To be copied in the server
# mod_sql_server("sql_1")
