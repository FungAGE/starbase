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
    dashboardPage(
      skin = "midnight",
      md = TRUE,
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
      sidebar = dashboardSidebar(
        minified = FALSE, collapsed = TRUE,
        sidebarMenu(
          menuItem("Welcome to starbase",
            tabName = "home", icon = NULL, badgeLabel = NULL, badgeColor = "green",
            href = NULL, newtab = FALSE, selected = NULL,
            expandedName = as.character(gsub("[[:space:]]", "", "Welcome to Starbase")),
            startExpanded = FALSE
          ),
          menuItem("Wiki",
            tabName = "wiki", icon = NULL, badgeLabel = NULL, badgeColor = "green",
            href = NULL, newtab = FALSE, selected = NULL,
            expandedName = as.character(gsub("[[:space:]]", "", "Wiki")),
            startExpanded = FALSE
          ),
          menuItem("BLAST/HMMER Searches",
            tabName = "blast", icon = NULL, badgeLabel = NULL, badgeColor = "green",
            href = NULL, newtab = FALSE, selected = NULL,
            expandedName = as.character(gsub("[[:space:]]", "", "BLAST/HMMER Searches")),
            startExpanded = FALSE
          ),
          menuItem("Explore Starships",
            tabName = "explore", icon = NULL, badgeLabel = NULL, badgeColor = "green",
            href = NULL, newtab = FALSE, selected = NULL,
            expandedName = as.character(gsub("[[:space:]]", "", "Explore Starships")),
            startExpanded = FALSE
          ),
          menuItem("Starship Synteny",
            tabName = "synteny", icon = NULL, badgeLabel = NULL, badgeColor = "green",
            href = NULL, newtab = FALSE, selected = NULL,
            expandedName = as.character(gsub("[[:space:]]", "", "Starship Synteny")),
            startExpanded = FALSE
          ),
          menuItem("starfish",
            icon = NULL, tabName = "starfish", badgeLabel = NULL, badgeColor = "green",
            href = NULL, newtab = FALSE, selected = NULL,
            expandedName = as.character(gsub("[[:space:]]", "", "starfish")),
            startExpanded = FALSE
          ),
          menuItem("Submit Starships to starbase",
            tabName = "submit", icon = NULL, badgeLabel = NULL, badgeColor = "green",
            href = NULL, newtab = FALSE, selected = NULL,
            expandedName = as.character(gsub("[[:space:]]", "", "Submit Starships to starbase")),
            startExpanded = FALSE
          ),
          menuItem("Update starbase Entries", 
            tabName = "db_update", icon = icon("th")
            ),
          id = NULL, .list = NULL
        )
      ),
      body = dashboardBody(
        shinyjs::useShinyjs(),
          fluidPage(
            fluidRow(
              shinyauthr::loginUI("login",cookie_expiry = cookie_expiry),
              uiOutput("loginUI")
            )
          )
      )
    )
  )
}
