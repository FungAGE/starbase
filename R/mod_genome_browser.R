#' genome_browser UI Function
#'
#' @description A shiny Module.
#'
#' @param id,input,output,session Internal parameters for {shiny}.
#'
#' @importFrom  shinyjs inlineCSS
#' @importFrom shiny NS tagList htmlOutput
#' @importFrom htmltools renderTags 
#' 
#' @noRd

mod_genome_browser_ui <- function(id) {
  ns <- NS(id)
  tagList(
    fluidPage(
      verticalLayout(
        fluidRow(
          inlineCSS(".form-group {margin-bottom: 0;}
                                .irs-with-grid {bottom: 0px;}
                                .irs-grid {height: 13px;}
                                .irs-grid-text {height: 0px;}
                                "),
          JBrowseROutput(ns("browserOutput"))
        )
      )
    )
  )
}

#' genome_browser Server Functions
#' @importFrom JBrowseR serve_data renderJBrowseR assembly track_feature tracks default_session JBrowseR JBrowseROutput
#'
#' @noRd
mod_genome_browser_server <- function(id) {
  moduleServer(id, function(input, output, session) {
    ns <- session$ns

    # TODO: make sure fasta headers match gff
    # TODO: make coordinates in gff relative to the element, rather than the chr. Use columns 4 and 5 in `mycodb.final.starships.feat` to correct. and remember, gff is 1-indexed
    path.fa <- "tmp/browser/aspcla3_s00849.fa.gz"
    path.gff <- "tmp/browser/aspcla3.gff.gz"
    g.chr <- "aspcla3_DS026990.1"
    mks.range.1 <- 1
    mks.range.2 <- 1000
    port <- httpuv::randomPort()

    if (!grepl("^http", path.fa)) {
      data_server <- serve_data(dirname(path.fa), port = port)
    } else {
      data_server <- NULL
    }

    if (!grepl("^http", path.fa)) {
      assembly <- assembly(
        paste0("http://127.0.0.1:", port, "/", basename(path.fa)),
        bgzip = TRUE
      )
    } else {
      assembly <- assembly(
        path.fa,
        bgzip = TRUE
      )
    }
    ## create configuration for a JB2 GFF FeatureTrack

    if (!is.null(path.gff)) {
      if (!grepl("^http", path.gff)) {
        annotations_track <- track_feature(
          paste0("http://127.0.0.1:", port, "/", basename(path.gff)),
          assembly
        )
      } else {
        annotations_track <- track_feature(
          path.gff,
          assembly
        )
      }
    } else {
      annotations_track <- NULL
    }

    ## create the tracks array to pass to browser
    tracks <- tracks(annotations_track)

    tracks_set <- c(annotations_track)

    theme <- JBrowseR::theme("#6c81c0", "#22284c")

    if (any(!is.null(tracks_set))) {
      default_session <- default_session(
        assembly,
        tracks_set[which(!is.null(tracks_set))]
      )
      output$browserOutput <- renderJBrowseR(JBrowseR(
        "View",
        assembly = assembly,
        tracks = tracks,
        location = paste0(g.chr, ":", mks.range.1, "..", mks.range.2),
        defaultSession = default_session,
        theme = theme
      ))
    } else {
      output$browserOutput <- renderJBrowseR(JBrowseR(
        "View",
        assembly = assembly,
        location = paste0(g.chr, ":", mks.range.1, "..", mks.range.2),
        theme = theme
      ))
    }
    # data_server$stop_server()
  })
}
