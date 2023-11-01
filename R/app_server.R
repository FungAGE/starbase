#' The application server-side
#'
#' @param input,output,session Internal parameters for {shiny}.
#'     DO NOT REMOVE.
#' @import shiny
#' @noRd
app_server <- function(input, output, session) {
  # Your application server logic
  mod_home_server("home_1")
  mod_wiki_server("wiki_1")
  mod_explore_server("explore_1")
  mod_blast_server("blast_1")
  # mod_genome_browser_server("genome_browser_1")
  mod_blast_syn_viz_server("blast_syn_viz_1")
  mod_starfish_server("starfish_1")
  mod_submit_server("submit_1")
  mod_user_server("user_1")
}
