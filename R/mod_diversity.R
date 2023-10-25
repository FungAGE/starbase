#' diversity UI Function
#'
#' @description A shiny Module.
#'
#' @param id,input,output,session Internal parameters for {shiny}.
#'
#' @noRd 
#'
#' @importFrom shiny NS tagList 
mod_diversity_ui <- function(id){
  ns <- NS(id)
  tagList(
    img(src = "img/heat_tree.png",
    width = "100%",
    style = "background-color: black;")
  )
}
    
#' diversity Server Functions
#'
#' @noRd 
mod_diversity_server <- function(id){
  moduleServer( id, function(input, output, session){
    ns <- session$ns
  })
}
    
## To be copied in the UI
# mod_diversity_ui("diversity_1")
    
## To be copied in the server
# mod_diversity_server("diversity_1")
