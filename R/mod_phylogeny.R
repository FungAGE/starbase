#' phylogeny UI Function
#'
#' @description A shiny Module.
#'
#' @param id,input,output,session Internal parameters for {shiny}.
#'
#' @noRd 
#'
#' @import tidyverse ggiraph stringi htmltools htmlwidgets crosstalk DT 
#' @importFrom shiny NS tagList 

# Load required packages
mod_phylogeny_ui <- function(id){
  ns <- NS(id)
  tagList(
    readRDS(file="data/captain-tree.RDS")
  )
}
    
#' phylogeny Server Functions
#'
#' @noRd 
mod_phylogeny_server <- function(id){
  moduleServer( id, function(input, output, session){
    ns <- session$ns
    
  })
}
