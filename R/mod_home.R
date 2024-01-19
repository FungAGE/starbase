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
    fluidRow(
      style = "padding-bottom: 15px;",
          tags$table(style = "width: 85%",
             tags$tr(tags$td(style = "width: 20%",
                             align = "right",
                             img(src = "img/favicon.png", width = "100%")),
                     tags$td(style = "width: 5%"),
                     tags$td(style = "width: 75%",
                             align = "left",
                      h2("Database and tools for exploring large eukaryotic transposable elements in Fungi"))))),
    fluidRow(
      box(
        title = "This is the development version of starbase", 
        status = "warning",
        width = 8,
        disabled(
          checkboxGroupInput(ns("functions"),
            label = "Working Functions:",
            choices = c("BLAST/HMMER searches", "Synteny/Genome Browser", "Catalogue/Wiki of Starship Metadata", "Running Starfish", "Submission of new Ships", "Editing of Database entries"),
            selected = c("BLAST/HMMER searches","Submission of new Ships", "Editing of Database entries")
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
