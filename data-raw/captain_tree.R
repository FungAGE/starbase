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

nexus_tree<-"../Starships/captain/tyr/faa/tree/funTyr50_cap25_crp3_p1-512_activeFilt.clipkit.new_colored.treefile"

# Read the Newick tree using the read.tree function
tree <- treeio::read.nexus(nexus_tree)
tree$tip.label<-gsub("'","",tree$tip.label)

# make base plot
p<-ggtree(tree)

# info on the groups we care about
family_nodes<-c(1834,1811,2032,2087,2153,1235,1384,1497,1643,1748,1778,1531,1606)
families<-read_tsv("../Starships/metadata/family/family-names.tsv",show_col_types = FALSE) %>%
  rename("groups"="longFamilyID")

group_df<-data.frame(parent=c(1834,1811,2032,2087,2153,1235,1384,1497,1643,1748,1778,1531,1606,2424,2425,2426,2442,2444,1477,1479,1487,1640,1227,1803),
  groups=c("superfam01-1","superfam01-2","superfam01-3","superfam01-4","superfam01-5","superfam02-1","superfam02-2","superfam02-3","superfam03-1","superfam03-2","superfam03-3","superfam03-4","superfam03-5",
  "extra group","extra group","extra group","extra group","extra group","extra group","extra group","extra group","extra group","extra group","extra group"),
  custom.fill=c("#8dd3c7","#ededa8","#adabc4","#33a02c","#fb8072","#80b1d3","#b3de69","#b0b0b0","#fdb45a","#fccde5","#ffed6f","#bc80bd","#ccebc5","black","black","black","black","black","black","black","black","black","black","black")) %>%
  left_join(families)

# create a annotated plot first
# layer data will be pulled out later to make interactive rects

p1<-p %<+% group_df +
  geom_hilight(type="rect",node=1834) + 
  geom_hilight(type="rect",node=1811) + 
  geom_hilight(type="rect",node=2032) + 
  geom_hilight(type="rect",node=2087) + 
  geom_hilight(type="rect",node=2153) + 
  geom_hilight(type="rect",node=1235) + 
  geom_hilight(type="rect",node=1384) + 
  geom_hilight(type="rect",node=1497) + 
  geom_hilight(type="rect",node=1643) + 
  geom_hilight(type="rect",node=1748) + 
  geom_hilight(type="rect",node=1778) + 
  geom_hilight(type="rect",node=1531) + 
  geom_hilight(type="rect",node=1606) + 
  geom_hilight(type="rect",node=2424) + 
  geom_hilight(type="rect",node=2425) + 
  geom_hilight(type="rect",node=2426) + 
  geom_hilight(type="rect",node=2442) + 
  geom_hilight(type="rect",node=2444) + 
  geom_hilight(type="rect",node=1477) + 
  geom_hilight(type="rect",node=1479) + 
  geom_hilight(type="rect",node=1487) + 
  geom_hilight(type="rect",node=1640) + 
  geom_hilight(type="rect",node=1227) + 
  geom_hilight(type="rect",node=1803)

tree_dat<-p$data %>% filter(parent %in% family_nodes) %>%
    left_join(group_df) %>%
    group_by(parent) %>%
  mutate(parentlab=ifelse(parent%in%family_nodes,parent,NA)) %>% 
  ungroup()

# hilights start at layer #3
highlighted_tree_dat<-map(3:(length(family_nodes)+2),~{layer_data(p1,.x)}) %>% list_rbind()

tree_dat<-tree_dat %>% 
    left_join(highlighted_tree_dat,by=c("parent"="clade_root_node"))

window_open = 'window.parent.open('
html_extension = '.html'
closing_quote = ')'

inter.tree.p<-p +
  geom_treescale(x=5,y=2,offset=5)+
  geom_rect_interactive(data=tree_dat,
    aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax,fill=groups,
    tooltip = groups,
    data_id = groups,
    size=1
    # onclick = sprintf('%s"%s%s"%s', !!window_open, groups, !!html_extension, !!closing_quote)
    ),
    alpha=0.5) + 
    geom_text(data=tree_dat,aes(x = x, y = y, label = groups), size = 5,nudge_x=5)+
    scale_fill_manual(values=group_df$custom.fill)+
    guides(size="none",fill="none")
    # scale_fill_manual(values=group_df$custom.fill, guide=guide_legend_interactive(title="Starship family_nodes"))
    # theme(legend.position="right")

ggsave(inter.tree.p,file="~/Downloads/captain-tree.svg")

captain_tree  <- girafe(
  ggobj = inter.tree.p,
  width_svg = 12,
  height_svg = 12,
  options = list(
    opts_hover(css = "fill:white;stroke:black;r:5pt;"),
    opts_sizing(width = 0.8),
    opts_zoom(min=0.7,max=10),
    opts_tooltip(opacity = .7,
      css = "text:black;",
      offx = 20, offy = -10,
      use_fill = TRUE, use_stroke = TRUE,
      delay_mouseout = 1000
    )
  )
)

usethis::use_data(captain_tree, overwrite = TRUE)