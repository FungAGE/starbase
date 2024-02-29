#' home UI Function
#'
#' @description A shiny Module.
#'
#' @param id,input,output,session Internal parameters for {shiny}.
#'
#' @noRd
#'
#' @import markdown
#' @importFrom shiny NS tagList
#' @importFrom shinyjs disabled
#'

working<-c("BLAST/HMMER searches", "Catalogue/Wiki of Starship Metadata", "Submission of new Ships", "Editing of Database entries")
notworking_ul <- working_ul <- tags$ul()
working_ul$children <- purrr::map(working, function(.x) tags$li(.x))
notworking<-c("Synteny/Genome Browser","Running Starfish")
notworking_ul$children <- purrr::map(working, function(.x) tags$li(.x))

mod_home_ui <- function(id) {
  ns <- NS(id)
  fluidPage(
    column(width=8,
    box(width=8,
    tags$table(style = "width: 85%",
      tags$tr(tags$td(style = "width: 85%",
                      align = "middle",
                      img(src = "img/favicon.png", width = "50%"))),
      tags$tr(tags$td(style = "width: 75%",
                      align = "middle",
                      h1("Database and tools for exploring large eukaryotic transposable elements in Fungi"))),
      tags$tr(tags$td(style = "width: 85%",
      align = "left",
        h4("This is the development version of starbase"))),
      tags$tr(tags$td(style = "width: 85%",
      align = "left",
        h4("Working Functions:"),
        working_ul
        )),
      tags$tr(tags$td(style = "width: 85%",
      align = "left",
      h4("Non-Working Functions:"),
      notworking_ul      
      )),
  ))))
}

#' home Server Functions
#'
#' @noRd
mod_home_server <- function(id) {
  moduleServer(id, function(input, output, session) {
    ns <- session$ns
  })
}
