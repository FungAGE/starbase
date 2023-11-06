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

mod_home_ui <- function(id) {
  ns <- NS(id)
  dashboardBody(
    div(
      fluidRow(
        style = "padding-bottom: 15px;",
        column(width = 1, img(src = "img/favicon.png", width = "100%")),
        column(width = 3, h2("starbase: Database and tools for exploring large eukaryotic transposable elements in Fungi"))
      )
    ),
    fluidRow(
      box(
        title = "This is the development version of starbase", status = "warning",
        disabled(
          checkboxGroupInput(ns("functions"),
            label = "Working Functions:",
            choices = c("BLAST/HMMER searches", "Synteny/Genome Browser", "Starship Catalogue", "Running Starfish", "Submission of new Ships", "Editing of Database entries"),
            selected = "BLAST/HMMER searches"
          )
        )
      )
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
