#' The application User-Interface
#'
#' @param request Internal parameter for `{shiny}`.
#'     DO NOT REMOVE.
#' @import shiny shinydashboard shinydashboardPlus bslib dplyr glue shinyauthr RSQLite DBI lubridate
#' @noRd

app_ui <- function(request) {
  tagList(
    # Leave this function for adding external resources
    golem_add_external_resources(),
    shinyjs::useShinyjs(),
    dashboardPage(
      skin = "blue",
      options = list(sidebarExpandOnHover = TRUE),
      header = dashboardHeader(title = "starbase",   
                              controlbarIcon = img(app_sys("img/favicon.ico")),
                              leftUi = tagList(
                                div(textOutput("welcome"), style = "padding: 10px"),
                                # TODO: should be on right side of page
                                tags$li(
                                  class = "dropdown",
                                  style = "padding: 8px;",
                                  shinyauthr::logoutUI("logout")
                                ),
                                # tags$li(
                                #   class = "dropdown",
                                #   tags$a(
                                #     icon("github"),
                                #     href = "https://github.com/FungAGE/starbase",
                                #     title = "See the code on github"
                                #   )
                                # )
                              )),
      controlbar = dashboardControlbar(),
      title = "starbase home",
      sidebar = uiOutput("loginSidebar"),
      body = dashboardBody(
              shinyauthr::loginUI("login",
                title = "Please log in to access starbase",
                cookie_expiry = cookie_expiry,
                additional_ui=mod_home_ui("home_1")
              ),
              uiOutput("loginBody")
            )
    )
  )
}
