#' The application User-Interface
#'
#' @param request Internal parameter for `{shiny}`.
#'     DO NOT REMOVE.
#' @import shiny shinydashboard shinydashboardPlus
#' @noRd

app_ui <- function(request) {
    tagList(
        # Leave this function for adding external resources
      golem_add_external_resources(),
      fluidPage(
            tagList(
        # Your application UI logic
        # includeMarkdown("dashboard.Rmd")

        dashboardPage(
          skin="blue",
          options = list(sidebarExpandOnHover = TRUE),
          header = dashboardHeader(title="starbase"),
          sidebar = dashboardSidebar(minified = FALSE, collapsed = FALSE,
            sidebarMenu(
              menuItem("Welcome to starbase", tabName = "home", icon = NULL, badgeLabel = NULL, badgeColor = "green",
                href = NULL, newtab = FALSE, selected = NULL,
                expandedName = as.character(gsub("[[:space:]]", "", "Welcome to Starbase")),
                startExpanded = TRUE),
              menuItem("BLAST/HMMER Searches", tabName = "blast", icon = NULL, badgeLabel = NULL, badgeColor = "green",
                href = NULL, newtab = FALSE, selected = NULL,
                expandedName = as.character(gsub("[[:space:]]", "", "BLAST/HMMER Searches")),
                startExpanded = FALSE),
              menuItem("Explore Starships", tabName = "explore", icon = NULL, badgeLabel = NULL, badgeColor = "green",
                href = NULL, newtab = FALSE, selected = NULL,
                expandedName = as.character(gsub("[[:space:]]", "", "Explore Starships")),
                startExpanded = FALSE),
              # menuItem("Starship Browser", tabName = "browse", icon = NULL, badgeLabel = NULL, badgeColor = "green",
              #   href = NULL, newtab = FALSE, selected = NULL,
              #   expandedName = as.character(gsub("[[:space:]]", "", "Starship Browser")),
              #   startExpanded = FALSE),
              menuItem("Starship Synteny", tabName = "synteny", icon = NULL, badgeLabel = NULL, badgeColor = "green",
                href = NULL, newtab = FALSE, selected = NULL,
                expandedName = as.character(gsub("[[:space:]]", "", "Starship Synteny")),
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
                  menuSubItem("Update starbase Entries", tabName = "db_update")),
              id = NULL, .list = NULL)

          ),
          body = dashboardBody(
                      tabItems(
                        tabItem(tabName = "home",
                          mod_home_ui("home_1")
                        ),
                        tabItem(tabName = "blast",
                          mod_blast_ui("blast_1"),
                          mod_blast_viz_ui("blast_viz_1")
                        ),
                        tabItem(tabName = "explore",
                          mod_explore_ui("explore_1")                        
                        ),
                        # tabItem(tabName = "browse",
                        #   mod_genome_browser_ui("genome_browser_1")
                        # ),
                        tabItem(tabName = "synteny",
                          mod_blast_syn_viz_ui("blast_syn_viz_1")
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
                    ),
          controlbar = dashboardControlbar(),
          title = "starbase home")
          )
    ),
  )
}