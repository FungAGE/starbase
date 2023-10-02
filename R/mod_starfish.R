#' starfish UI Function
#'
#' @description A shiny Module.
#'
#' @param id,input,output,session Internal parameters for {shiny}.
#'
#' @noRd 
#'
#' @importFrom shiny NS tagList 
mod_starfish_ui <- function(id){
  ns <- NS(id)
  dashboardBody(
    div(
        div(h1("Annotate Starships with starfish"),
          img(src = "img/STARFISH_LOGO.png",
              width = "15%")),
        p("starfish is a modular toolkit for giant mobile element annotation. Built primarily for annotating Starship elements in fungal genomes, it can be easily adapted to find any large mobile element (≥6kb) that shares the same basic architecture as a fungal Starship or a bacterial integrative and conjugative element: a 'captain' gene with zero or more 'cargo genes downstream of its 3' end. It is particularly well suited for annotating low-copy number elements in a content independent manner.")
    )
  )
}
    
#' starfish Server Functions
#'
#' @noRd 
mod_starfish_server <- function(id){
  moduleServer( id, function(input, output, session){
    ns <- session$ns
 
  })
}

