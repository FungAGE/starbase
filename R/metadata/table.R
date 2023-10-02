library(tidyverse)
library(htmltools)
library(htmlwidgets)
library(crosstalk)
library(DT)

dat<-read_tsv("/home/adrian/Systematics/Starship_Database/MTDB/starships.db",col_names=c("ome","genus","species","strain","taxonomy","version","source","biosample","assembly_acc","acquisition_date","published","fna","faa","gff3")) %>%
  separate_rows(taxonomy,sep=",") %>%
  mutate(taxonomy=gsub(c("\\{|\\}|\\'"),"",trimws(taxonomy))) %>%
  separate(taxonomy,sep=": ",into=c("rank","name")) %>%
  pivot_wider(id_cols=c("ome","genus","species","strain","version","source","biosample","assembly_acc","acquisition_date","published","fna","faa","gff3"),names_from="rank",values_from="name") %>%
  relocate("kingdom","clade","phylum","subphylum","class","subclass","order","suborder","family","subfamily","genus","species") %>%
  select(-c(kingdom:subfamily))

dat.html<-datatable(dat, options = list(), class = "display",
    callback = JS("return table;"), #rownames, colnames, container,
    caption = NULL, filter = c("none", "bottom", "top"), escape = TRUE,
    style = "auto", width = NULL, height = NULL, elementId = NULL,
    fillContainer = getOption("DT.fillContainer", NULL),
    autoHideNavigation = getOption("DT.autoHideNavigation", NULL),
    selection = c("multiple", "single", "none"), extensions = list(),
    plugins = NULL, editable = FALSE)

htmlwidgets::saveWidget(dat.html, "table.html")
