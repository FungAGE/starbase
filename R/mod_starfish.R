#' starfish UI Function
#'
#' @description A shiny Module.
#'
#' @param id,input,output,session Internal parameters for {shiny}.
#'
#' @noRd
#'
#' @importFrom shiny NS tagList

# TODO: add functionality for multiple inputs?
mod_starfish_ui <- function(id) {
  ns <- NS(id)
  tagList(
    box(
      title = "Annotate Starships with starfish", icon = img(src = "img/STARFISH_LOGO.png", width = "15%"),
      p("starfish is a modular toolkit for giant mobile element annotation. Built primarily for annotating Starship elements in fungal genomes, it can be easily adapted to find any large mobile element (≥6kb) that shares the same basic architecture as a fungal Starship or a bacterial integrative and conjugative element: a 'captain' gene with zero or more 'cargo genes downstream of its 3' end. It is particularly well suited for annotating low-copy number elements in a content independent manner.")
    ),
    box(
      h1("Run starfish on your genome"),
      fileInput(ns("input_fasta"), "Upload a genome file in FASTA format", width = "100px"),
      fileInput(ns("input_gff"), "Upload annotations for this genome in GFF3 format", width = "100px"),
      actionButton(ns("starfish"), "Launch starfish")
    )
  )
}

#' starfish Server Functions
#'
#' @noRd
mod_starfish_server <- function(id) {
  moduleServer(id, function(input, output, session) {
    ns <- session$ns
    data <- reactive({
      validate(
        need(input$input_fasta != "", "Provide a genome in FASTA format"),
        need(input$input_gff != "", "Provide annotations in GFF3 format")
      )

      ext <- tools::file_ext(input$input_fasta$name)
      validate(need(!ext %in% c("fasta", "fa", "fna"), "Provide a genome in FASTA format"))
    })
  })
}
