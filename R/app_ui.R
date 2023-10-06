#' The application User-Interface
#'
#' @param request Internal parameter for `{shiny}`.
#'     DO NOT REMOVE.
#' @import shiny
#' @noRd
library(shinythemes)
app_ui <- function(request) {
    tagList(
        # Leave this function for adding external resources
      golem_add_external_resources(),
      # tags$head(
      #   
      # h1("starbase"),

      fluidPage(theme = shinytheme("cerulean"),
        mod_dashboard_ui("dashboard_1")
    ),
  )
}