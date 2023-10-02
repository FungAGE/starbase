#' explore UI Function
#'
#' @description A shiny Module.
#'
#' @param id,input,output,session Internal parameters for {shiny}.
#'
#' @noRd 
#'
#' @importFrom shiny NS tagList 
mod_explore_ui <- function(id){
  ns <- NS(id)
  families<-c("superfam01-1","superfam01-2","superfam01-3","superfam01-4","superfam01-5","superfam02-1","superfam02-2","superfam02-3","superfam03-1","superfam03-2","superfam03-3","superfam03-4","superfam03-5")
  fluidRow(div(
    h1("Captain Phylogeny"),
    box(mod_phylogeny_ui("phylogeny_1")),
    box(mod_metadata_ui("metadata_1"))
  )
  )
}
    
#' explore Server Functions
#'
#' @noRd 
mod_explore_server <- function(id){
  moduleServer( id, function(input, output, session){
    ns <- session$ns
    mod_phylogeny_server("phylogeny_1")
    mod_metadata_server("metadata_1")
  })
}
