#' submit UI Function
#'
#' @description A shiny Module.
#'
#' @param id,input,output,session Internal parameters for {shiny}.
#'
#' @noRd
#'
#' @importFrom shiny NS tagList
mod_submit_ui <- function(id) {
  ns <- NS(id)
  tagList()
}

#' submit Server Functions
#'
#' @noRd
mod_submit_server <- function(id) {
  moduleServer(id, function(input, output, session) {
    ns <- session$ns
  })
}
