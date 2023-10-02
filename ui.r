library(shinythemes)
library(DT)
library(golem)

ui <-fluidPage(theme = shinytheme("cerulean"),
            tagList(
              tags$head(
                tags$link(rel="stylesheet", type="text/css",href="style.css"),
                tags$script(type="text/javascript", src = "busy.js"),
                    # input block
                    mainPanel(
                      headerPanel('Starship BLAST/HMMER searches'),
                      textAreaInput('query', 'Input sequence:', value = "", placeholder = "", width = "600px", height="200px"),
                      # choice should be limited to 1) whole starship or 2) just a gene sequence
                      # if 2) then choose what gene or unknown
                      selectInput("input_type", "What is the input sequence?", choices=c("starship","gene"), width="120px"),
                      # BUG: where db_type when choice = gene is not overwritten from default db_type when choice = starship
                      conditionalPanel(
                        condition = "input.input_type == 'gene'",
                        selectInput("db_type", "Select a gene model:", selected = character(0), multiple = TRUE, choices = c("tyr", "freB", "nlr", "DUF3723","plp"), width="120px")
                      ),
                      conditionalPanel(
                        condition = "input.input_type == 'starship'",
                        selectInput("db_type", "Which database do you want to search?",selected = character(0), multiple = FALSE,choices = c("curated", "starfish"), width="120px"),
                        selectInput("search_ship_genes", "Search for genes in starship sequence?",choices = c("No", "Yes"), width="120px")
                      ),
                      div(style="display:inline-block",
                          selectInput("eval", "e-value:", choices=c(1,0.001,1e-4,1e-5,1e-10), width="120px")),
                      actionButton("blast", "Search")
                    ),
                    
                    # this snippet generates a progress indicator for long BLASTs
                    # FIXME: always spinning
                    # div(class = "busy",  
                    #     p("Calculation in progress.."), 
                    #     img(src="https://i.stack.imgur.com/8puiO.gif", height = 100, width = 100,align = "center")
                    # ),

                #Basic results output
                # mainPanel(
                  h4("Results"),
                  DTOutput('tbl')
                  # p("Alignment:", tableOutput("clicked") ),
                  # verbatimTextOutput("alignment")
                # )
    )
  )
)