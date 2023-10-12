#' blast_viz UI Function
#'
#' @description A shiny Module.
#'
#' @param id,input,output,session Internal parameters for {shiny}.
#'
#' @noRd 
#'
#' @importFrom shiny NS tagList 

# TODO: use Biostrings here and then link to modules for making alignments, circos, dot plots, and jbrowse

mod_blast_viz_ui <- function(id){
  ns <- NS(id)
  tagList(
      golem_add_external_resources(),
      fluidPage(
        h1("BlasterJS Visualization in Shiny"),
        # uiOutput("html_output")
        # includeHTML(readChar(app_sys("www/BlasterJS.html"), file.info(app_sys("www/BlasterJS.html"))$size))
        shiny::htmlTemplate("inst/app/www/BlasterJS.html")
      )

        # Div to hold the BlasterJS visualization
        # div(id = "blast-container", style = "height: 600px;",
        #     shiny::htmlTemplate(app_sys("www/test.html")))
      )
}

golem_add_external_resources <- function() {
  add_resource_path("www",app_sys("app/www"))
  add_resource_path( 'img', app_sys('app/img'))
  add_resource_path( 'js', app_sys('app/www/js'))
  tagList(
    tags$head(
      favicon(),
      bundle_resources(
        path = app_sys("app/www"),
        app_title = "starbase"
      ),
      # Add here other external resources
      # for example, you can add shinyalert::useShinyalert()
      tags$link(rel="stylesheet", type="text/css",href="www/custom.css")
    )
  )
}


#' blast_viz Server Functions
#'
#' @noRd 
mod_blast_viz_server <- function(id){
  moduleServer( id, function(input, output, session){
    ns <- session$ns
    # output$html_output <- renderUI({
    #   includeHTML(app_sys("www/BlasterJS.html"))
    # })
    }
  )
}