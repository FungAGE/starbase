#' home UI Function
#'
#' @description A shiny Module.
#'
#' @param id,input,output,session Internal parameters for {shiny}.
#'
#' @noRd 
#'
#' @importFrom shiny NS tagList 
mod_home_ui <- function(id){
  ns <- NS(id)
  dashboardBody(
    fluidRow(
      box(icon=img(src="www/favicon.ico"),
        p("starbase: Database and tools for exploring large eukaryotic transposable elements in Fungi")),
      box(div(
        img(src = "img/starship-summary.png",
            width = "75%",
            style = "background-color: white;"
          )
        )
      )
    )
  )
}
    
#' home Server Functions
#'
#' @noRd 
mod_home_server <- function(id){
  moduleServer( id, function(input, output, session){
    ns <- session$ns
 
  })
}