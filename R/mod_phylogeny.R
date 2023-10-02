#' phylogeny UI Function
#'
#' @description A shiny Module.
#'
#' @param id,input,output,session Internal parameters for {shiny}.
#'
#' @noRd 
#'
#' @importFrom shiny NS tagList 

# Load required packages
library(tidyverse)
library(ggtree)
library(ggiraph)
library(treeio)
library(stringi)
library(htmltools)
library(htmlwidgets)
library(crosstalk)
library(DT)

mod_phylogeny_ui <- function(id){
  ns <- NS(id)
  tagList(
    readRDS(file="RDS/captain-tree.RDS")
  )
}
    
#' phylogeny Server Functions
#'
#' @noRd 
mod_phylogeny_server <- function(id){
  moduleServer( id, function(input, output, session){
    ns <- session$ns
    
  })
}