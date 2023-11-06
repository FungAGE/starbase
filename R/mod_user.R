#' user UI Function
#'
#' @description A shiny Module.
#'
#' @param id,input,output,session Internal parameters for {shiny}.
#'
#' @noRd
#'
#' @importFrom shiny NS tagList
mod_user_ui <- function(id) {
  ns <- NS(id)
  tagList(
    dashboardBody(
      box(title = "starbase User Account: ", userOutput("user"))
    )
  )
}

#' user Server Functions
#'
#' @noRd
mod_user_server <- function(id) {
  moduleServer(id, function(input, output, session) {
    ns <- session$ns
    output$user <- renderUser({
      dashboardUser(
        name = "Adrian",
        image = "",
        title = "shinydashboardPlus",
        subtitle = "Author",
        footer = p("The footer", class = "text-center"),
      )
    })
  })
}
