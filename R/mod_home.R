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
    fluidRow(
      box(width=12,
        tags$table(style = "width: 85%",
          tags$tr(tags$td(style = "width: 85%",
                          align = "middle",
                          h1("starbase: A database and toolkit for exploring large eukaryotic transposable elements in Fungi"))),
          tags$tr(tags$td(style = "width: 100%",
                          align = "middle",
                          img(src = "img/starbase-map.png", width = "75%"))),
          tags$tr(tags$td(style = "width: 85%",
          align = "left",
            h4("This version of starbase is under active development."))),
          tags$tr(tags$td(style = "width: 85%",
          align = "left",
            h4("Working Functions:"),
            working_ul
            )),
          tags$tr(tags$td(style = "width: 85%",
          align = "left",
          h4("Non-Working Functions:"),
          notworking_ul      
          ))
      ))
    )
  )
}

#' home Server Functions
#'
#' @noRd
mod_home_server <- function(id) {
  moduleServer(id, function(input, output, session) {
    ns <- session$ns
  })
}
