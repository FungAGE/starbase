#' dashboard UI Function
#'
#' @description main dash board to starbase
#'
#' @param id,input,output,session Internal parameters for {shiny}.
#'
#' @noRd 
#'
#' @importFrom shiny NS tagList 

library(shiny)
library(shinydashboard)
library(shinydashboardPlus)

mod_dashboard_ui <- function(id){
  ns <- NS(id)
    tagList(
        # Your application UI logic
        # includeMarkdown("dashboard.Rmd")

        dashboardPage(
          options = list(sidebarExpandOnHover = TRUE),
          header = dashboardHeader(title="starbase"),
          sidebar = dashboardSidebar(minified = FALSE, collapsed = FALSE,
            sidebarMenu(
              menuItem("Welcome to starbase", tabName = "home", icon = NULL, badgeLabel = NULL, badgeColor = "green",
                href = NULL, newtab = FALSE, selected = NULL,
                expandedName = as.character(gsub("[[:space:]]", "", "Welcome to Starbase")),
                startExpanded = TRUE),
              menuItem("Explore Starships", tabName = "explore", icon = NULL, badgeLabel = NULL, badgeColor = "green",
                href = NULL, newtab = FALSE, selected = NULL,
                expandedName = as.character(gsub("[[:space:]]", "", "Explore Starships")),
                startExpanded = FALSE),
              menuItem("BLAST/HMMER Searches", tabName = "blast", icon = NULL, badgeLabel = NULL, badgeColor = "green",
                href = NULL, newtab = FALSE, selected = NULL,
                expandedName = as.character(gsub("[[:space:]]", "", "BLAST/HMMER Searches")),
                startExpanded = FALSE),
              menuItem("starfish", icon = NULL, tabName = "starfish", badgeLabel = NULL, badgeColor = "green",
                href = NULL, newtab = FALSE, selected = NULL,
                expandedName = as.character(gsub("[[:space:]]", "", "starfish")),
                startExpanded = FALSE),
              menuItem("Submit Starships to starbase", tabName = "submit", icon = NULL, badgeLabel = NULL, badgeColor = "green",
                href = NULL, newtab = FALSE, selected = NULL,
                expandedName = as.character(gsub("[[:space:]]", "", "Submit Starships to starbase")),
                startExpanded = FALSE),
              menuItem("Sign in", tabName = "user", icon = NULL, badgeLabel = NULL, badgeColor = "green",
                       href = NULL, newtab = FALSE, selected = NULL,
                       expandedName = as.character(gsub("[[:space:]]", "", "Sign in")),
                       startExpanded = FALSE,
                       menuSubItem("Update starbase Entries", tabName = "db_update")
                       ),
              id = NULL, .list = NULL)

          ),
          body = dashboardBody(
                      tabItems(
                        tabItem(tabName = "home",
                          mod_home_ui("home_1")
                        ),
                        tabItem(tabName = "explore",
                          mod_explore_ui("explore_1")                        
                        ),
                        tabItem(tabName = "blast",
                          mod_blast_ui("blast_1")
                          # mod_blast_viz_ui("blast_viz_1")
                        ),
                        tabItem(tabName = "starfish",
                          mod_starfish_ui("starfish_1")
                        ),
                        tabItem(tabName = "submit",
                          mod_submit_ui("submit_1")
                        ),
                        tabItem(tabName = "user",
                          mod_user_ui("user_1")
                        ),
                        tabItem(tabName = "db_update",
                          mod_db_update_ui("db_update_1")
                        )
                      )
                    )

          ),
          controlbar = dashboardControlbar(),
          title = "starbase home")
}
    
#' dashboard Server Functions
#'
#' @noRd 
mod_dashboard_server <- function(id){
  moduleServer( id, function(input, output, session){
    ns <- session$ns
    function(input, output, session) {
      mod_home_server("home_1")
      mod_explore_server("explore_1")
      mod_blast_server("blast_1")
      # mod_blast_viz_server("blast_viz_1")
      mod_starfish_server("starfish_1")
      mod_submit_server("submit_1")
      mod_user_server("user_1")
    }
  })
}