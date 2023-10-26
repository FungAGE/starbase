#' dotplot_syn UI Function
#'
#' @description A shiny Module.
#'
#' @param id,input,output,session Internal parameters for {shiny}.
#'
#' @noRd 
#'
#' @importFrom shiny NS tagList 
mod_dotplot_syn_ui <- function(id){
  ns <- NS(id)
  tabPanel(
    "Dot View",
    fluidRow(
      column(
        6,
        fluidRow(
          column(
            12,
            div(id="dotView")
          )
        ),
        fluidRow(
          class = "justify-content-end",
          div(class = "col-sm-auto",
              style = "padding-bottom: 7.5px;",
              shinyjs::hidden(
                actionButton(
                  inputId = ns("dotView_download"),
                  status = "secondary",
                  icon = icon("download"),
                  label = "Dot View PNG"
                )
                ##downloadButton_custom(
                ##    "dotView_download",
                ##    status = "secondary",
                ##    icon = icon("download"),
                ##    label = "Dot View SVG"
                ##)
              )
          )
        )
      ),
      column(
        6,
        div(id="dotView_geneTable",
            class = "table-responsive",
            ## div(id="dotView_tableHeader",
            ##     h4("Anchor Genes:"),
            ##     p(style="color: gray;",
            ##       "Please selected your region of interest from the dot plot on the left panel, the table below will be updated automatically.")
            ## ),
            ## DTOutput("selected_anchors"))
            uiOutput("dotviewTable"))
      )
    ),
    icon = icon("microscope")
  )
}
    
#' dotplot_syn Server Functions
#'
#' @noRd 
mod_dotplot_syn_server <- function(id){
  moduleServer( id, function(input, output, session){
    ns <- session$ns
 
  })
}
    
## To be copied in the UI
# mod_dotplot_syn_ui("dotplot_syn_1")
    
## To be copied in the server
# mod_dotplot_syn_server("dotplot_syn_1")
