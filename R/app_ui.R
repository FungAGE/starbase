#' The application User-Interface
#'
#' @param request Internal parameter for `{shiny}`.
#'     DO NOT REMOVE.
#' @import shiny shinydashboard shinydashboardPlus bslib
#' @noRd

app_ui <- function(request) {
  tagList(
    # Leave this function for adding external resources
    golem_add_external_resources(),
    shinyjs::useShinyjs(),
    dashboardPage(
      skin = "blue",
      options = list(sidebarExpandOnHover = TRUE),
      header = dashboardHeader(title = "starbase"#,controlbarIcon = img(src="img/favicon.ico")
      ),
      controlbar = dashboardControlbar(disable=TRUE),
      title = "starbase home",
      sidebar = dashboardSidebar(
        minified = FALSE,
        collapsed = FALSE,
        sidebarMenu(
          menuItem("Home",
            tabName = "home",
            href = NULL, newtab = FALSE, selected = TRUE,
            expandedName = as.character(gsub("[[:space:]]", "", "home")),
            startExpanded = FALSE
          ),
          menuItem("Wiki",
            tabName = "wiki",
            href = NULL, newtab = FALSE, selected = NULL,
            expandedName = as.character(gsub("[[:space:]]", "", "Wiki")),
            startExpanded = TRUE
          ),
          menuItem("BLAST/HMMER Searches",
            tabName = "blast",
            href = NULL, newtab = FALSE, selected = NULL,
            expandedName = as.character(gsub("[[:space:]]", "", "BLAST/HMMER Searches")),
            startExpanded = FALSE
          ),
          menuItem("Explore Starships",
            tabName = "explore",
            href = NULL, newtab = FALSE, selected = NULL,
            expandedName = as.character(gsub("[[:space:]]", "", "Explore Starships")),
            startExpanded = FALSE
          ),
          menuItem("Starship Browser",
            tabName = "genome_browser",
            href = NULL, newtab = FALSE, selected = NULL,
            expandedName = as.character(gsub("[[:space:]]", "", "Starship Browser")),
            startExpanded = FALSE
          ),
          menuItem("starfish",
            icon = NULL, tabName = "starfish",
            href = NULL, newtab = FALSE, selected = NULL,
            expandedName = as.character(gsub("[[:space:]]", "", "starfish")),
            startExpanded = FALSE
          ),
          menuItem("Submit Starships to starbase",
            tabName = "submit",
            href = NULL, newtab = FALSE, selected = NULL,
            expandedName = as.character(gsub("[[:space:]]", "", "Submit Starships to starbase")),
            startExpanded = FALSE
          ),
          menuItem("About starbase",
            tabName = "about",
            href = NULL, newtab = FALSE, selected = NULL,
            expandedName = as.character(gsub("[[:space:]]", "", "About starbase")),
            startExpanded = FALSE
          ),

          id = NULL, .list = NULL
        )),
        body=dashboardBody(
          tabItems(
            tabItem(
              tabName = "home",
              mod_home_ui("home_1")
            ),
            tabItem(
              tabName = "wiki",
              mod_wiki_ui("wiki_1")
            ),
            tabItem(
              tabName = "blast",
              mod_blast_ui("blast_1")
            ),
            tabItem(
              tabName = "explore",
              mod_explore_ui("explore_1")
            ),
            tabItem(
              tabName = "genome_browser",
              mod_genome_browser_ui("genome_browser_1")
            ),
            tabItem(
              tabName = "starfish",
              mod_starfish_ui("starfish_1")
            ),
            tabItem(
              tabName = "submit",
              mod_submit_ui("submit_1")
            ),
            tabItem(
              tabName = "about",
              mod_about_ui("about_1")
            )
          )
        )
    )
  )
}
