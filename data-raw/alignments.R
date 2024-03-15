library(tidyverse)
library(Biostrings)
library(msaR)

# TODO: do this for all alignments?
# ! alignment is too large to load 

# read some sequences from a multiple sequence alignment file and display
aln<-readAAMultipleAlignment("Starships/genes/cap_tyr/alignments/YRsuperfamRefs.mafft.faa") %>%
  msaR(., seqlogo = T,menu=F, overviewbox = F,  colorscheme = "clustal")

aln %>%
  saveRDS("data/YRsuperfamRefs.mafft.faa.RDS")
  
# aln %>%
#   htmlwidgets::saveWidget(file=app_sys("html/YRsuperfamRefs.mafft.faa.html"))

# example of integrating this with shiny:
# shiny output binding for a widget named 'foo'
# fooOutput <- function(outputId, width = "100%", height = "400px") {
#   htmlwidgets::shinyWidgetOutput(outputId, "foo", width, height)
# }
# 
# # shiny render function for a widget named 'foo'
# renderFoo <- function(expr, env = parent.frame(), quoted = FALSE) {
#   if (!quoted) { expr <- substitute(expr) } # force quoted
#   htmlwidgets::shinyRenderWidget(expr, fooOutput, env, quoted = TRUE)
# }
