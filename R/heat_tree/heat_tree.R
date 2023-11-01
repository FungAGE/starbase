remotes::install_local("/home/adrian/Systematics/bin/metacoder",force=TRUE)

library(tidyverse)
library(metacoder)
library(ggiraph)

dat<-read_tsv("MTDB/joined_ships.tsv") %>%
  # distinct(checksum,.keep_all=TRUE) %>%
  rowwise() %>%
  mutate(genus=ifelse(is.na(genus),  gsub(" (.+)","",species),genus)) %>%
  group_by(genus) %>%
  fill(c(kingdom,phylum,class,order,family,genus),.direction="downup") %>% # fill in missing taxonomic info
  group_by(species) %>%
  fill(taxid,.direction="downup") %>%
  ungroup() %>%
  mutate(across(c(kingdom,phylum,class,order,family,genus,species),~replace_na(.,"Unknown"))) %>%
  select(ome,kingdom,phylum,class,order,family,genus,species,taxid)

# fill in missing tax levels
missing_tax <- dat %>%
  filter(is.na(phylum)|is.na(class)|is.na(order)|is.na(family)) %>%
  distinct(taxid,genus,species,.keep_all=T) %>%
  mutate(tax=paste(genus,species))

# map(missing_tax$tax,~{
#   x<-system(paste0("echo 'Coccidioides immitis' | /usr/local/bin/taxonkit name2taxid | /usr/local/bin/taxonkit lineage -i 2 -r -R"),intern = TRUE) 
#   y<-x%>% read.csv(text=.,sep="\t",header = F)
#   colnames(y)<-c("tax","taxid","names","rank","ranks")
#   y %>% separate_rows(names,ranks,sep = ";") %>%
#     filter(ranks != "no rank") %>%
# })

# create taxmap object
obj<-dat %>% select(-taxid) %>% 
    parse_tax_data(., class_cols = 2:8,
                   named_by_rank = TRUE)

# Each number will produce a slightly different result for some layouts
set.seed(2) 

p<-obj %>%
  # filter_taxa(grepl(pattern = "^[a-zA-Z]+$", taxon_names)) %>% # remove "odd" taxa
  heat_tree(
    node_label = ifelse(taxon_ranks == "species", "", taxon_names),
    node_label_color = "black",
    node_size = n_obs,
    node_color = n_obs,
    # initial_layout = "re", 
    layout = "da",
    node_legend_title = "Number of Starships",
    node_size_axis_label = "",
    node_color_axis_label = "")

# save as png/svg
ggsave(plot=p,file="img/heat_tree.png")

interactive_plot  <- girafe(
  ggobj = p, width_svg = 15,
  height_svg = 15,
  options = list(
    opts_sizing(width = 0.8),
    opts_zoom(min=0.7,max=10),
    opts_tooltip(opacity = .7,
                 offx = 20, offy = -10,
                 use_fill = TRUE, use_stroke = TRUE,
                 delay_mouseout = 1000
    )
  )
)


# Save the interactive plot as an HTML file
htmlwidgets::saveWidget(interactive_plot, app_sys("html/heat_tree.html"),title="Starship Database Diversity")


saveRDS(p,"RDS/heat_tree.rds")