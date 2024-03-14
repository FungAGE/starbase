#' The application server-side
#'
#' @param input,output,session Internal parameters for {shiny}.
#'     DO NOT REMOVE.
#' @import shiny shinydashboard shinydashboardPlus tidyverse
#' @noRd 

app_server <- function(input, output, session) {
  mod_home_server("home_1")
  mod_wiki_server("wiki_1")
  mod_blast_server("blast_1")
  mod_genome_browser_server("genome_browser_1")
  mod_starfish_server("starfish_1")
  mod_submit_server("submit_1")
  mod_about_server("about_1")
}
