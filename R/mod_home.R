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

starbase_version <- "0.0.0.9"

working<-c("BLAST/HMMER searches", "Catalogue/Wiki of Starship Metadata", "Submission of new Ships")
notworking_ul <- working_ul <- tags$ul()
working_ul$children <- purrr::map(working, function(.x) tags$li(.x))
notworking<-c("Synteny/Genome Browser","Running Starfish")
notworking_ul$children <- purrr::map(notworking, function(.x) tags$li(.x))

# TODO: href features to pages?

mod_home_ui <- function(id) {
  ns <- NS(id)
  fluidPage(
    fluidRow(
      tags$table(style = "width: 85%",
        tags$tr(tags$td(style = "width: 100%",
                        align = "middle",
                        h1("starbase: A database and toolkit for exploring large eukaryotic transposable elements in Fungi"))),
        tags$tr(tags$td(style = "width: 100%",
                        align = "middle",
                        column(width=4,
                          box(width=NULL,title="What can I do with starbase?",status="primary",
                            working_ul)),
                        column(width=4,
                          box(width=NULL,title="This version of starbase is under active development",status="warning",
                            tags$table(style = "width: 85%",
                              tags$tr(tags$td(style = "width: 85%",
                                              align = "middle",
                                              h4("Non-Working Functions:"),
                                              notworking_ul))
                            )
                          )
                        ),
                        column(width=4,
                          box(width=NULL,
                              title="Contact info:",
                              status="success",
                              p("starbase was developed by the ", a("FungAGE lab", href="https://fungage.github.io/")),
                              p("code for starbase is available on GitHub ", 
                              socialButton(href="https://github.com/FungAGE/starbase",icon=icon("github")))
                        )))),
        tags$tr(tags$td(style = "width: 85%",
                        align = "middle",
                        img(src = "img/starbase-map.png", width = "100%"))),
        tags$tr(tags$td(style = "width: 85%",
                        align = "middle",
                        column(width=10,
                          box(width=NULL,title = "Data Availability",status="primary",
                            p("We have been maintaining starbase data on our GitHub repo. We are currently in the process of migrating to a new back-end, which will provide more options for data export. In the mean time, you can retrieve all Starship sequences, annotations, and more, in a single .zip file (size ~200Mb)"),
                            downloadButton(outputId="dl_package",label="Download the latest version of starbase.")
                          )
                        )
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
    filename <- paste0("starbase-v",starbase_version,".tar.gz")
    output$dl_package <- downloadHandler(
      filename = filename, 
      content = "https://github.com/FungAGE/Starships/archive/refs/heads/main.zip"
    )
  })
}
