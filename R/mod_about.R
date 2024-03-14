#' about UI Function
#'
#' @description A shiny Module.
#'
#' @param id,input,output,session Internal parameters for {shiny}.
#'
#' @noRd 
#'
#' @importFrom shiny NS tagList 

mod_about_ui <- function(id){
  ns <- NS(id)
  fluidPage(
    fluidRow(
      box(
        title = "About starbase",
        solidHeader = TRUE,
        collapsible = TRUE,
        width = NULL,
        p("starbase was developed by the ", a("FungAGE lab", href="https://fungage.github.io/")),
        p("code for starbase is available on ", a("GitHub", href="https://github.com/FungAGE/starbase"))
      )
    ),
    fluidRow(
      userBox(title=userDescription(title = "Adrian Forsythe", subtitle = "Lead developer of starbase", type = 2, image = "img/adrian.png"),
                status = "teal",socialButton(href="mailto:adrian.e.forsythe@gmail.com", icon = icon("envelope")),collapsible=FALSE),
      userBox(title=userDescription(title = "Aaron Vogan", subtitle = "FungAGE group leader", type = 2, image = "img/aaron.png"),
                status = "purple",socialButton(href="mailto:aaron.vogan@ebc.uu.se", icon = icon("envelope")),collapsible=FALSE),
      userBox(title=userDescription(title = "Emile Gluck-Thaler", subtitle = "Developer of starfish", type = 2, image = "img/emile.png"),
                status = "orange",socialButton(href="mailto:emilegluckthaler@gmail.com", icon = icon("envelope")),collapsible=FALSE)
    )
  )
}
    
#' about Server Functions
#'
#' @noRd 
mod_about_server <- function(id){
  moduleServer( id, function(input, output, session){
    ns <- session$ns
 
  })
}
    
## To be copied in the UI
# mod_about_ui("about_1")
    
## To be copied in the server
# mod_about_server("about_1")
