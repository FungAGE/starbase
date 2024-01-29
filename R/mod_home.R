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
  shiny::tagList(
    br(),
    box(
      title = "Database and tools for exploring large eukaryotic transposable elements in Fungi", 
      status = "warning",
      width = NULL,
      tags$table(style = "width: 85%",
        tags$tr(tags$td(style = "width: 20%",
                        align = "middle",
                        img(src = "img/favicon.png", width = "100px")),
                tags$td(style = "width: 5%"),
                tags$td(style = "width: 20%",
                        align = "left",
            h4("This is the development version of starbase")))),
      footer = fluidRow(
        column(width=12,
          h2(),
          h4("Working Functions:"),
          p("BLAST/HMMER searches"),
          br(),
          p("Catalogue/Wiki of Starship Metadata"),
          br(),
          p("Submission of new Ships"),
          br(),
          p( "Editing of Database entries"),
          br(),
          h4("Non-Working Functions:"),
          p("Synteny/Genome Browser"),
          br(),
          p("Running Starfish")
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
