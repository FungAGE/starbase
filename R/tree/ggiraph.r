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

nexus_tree<-"/home/adrian/Systematics/Starship_Database/Starships/starships/starfish/phylo/funTyr50_cap25_crp3_p1-512_activeFilt.clipkit.new_colored.treefile"

# Read the Newick tree using the read.tree function
tree <- treeio::read.nexus(nexus_tree)
tree$tip.label<-gsub("'","",tree$tip.label)

# make base plot
p<-ggtree(tree)

# load taxonomic info
tax<-full_join(read_tsv("/home/adrian/Systematics/Starship_Database/MTDB/mycodb2899.ome2species.txt",col_names=c("code","species")),
               read_csv("/home/adrian/Systematics/Starship_Database/Starships/starships/manual-annotations/Starships.csv") %>% rename_all(tolower))
# info on the groups we care about
families<-c(1834,1811,2032,2087,2153,1235,1384,1497,1643,1748,1778,1531,1606)

group_df<-data.frame(parent=c(1834,1811,2032,2087,2153,1235,1384,1497,1643,1748,1778,1531,1606,2424,2425,2426,2442,2444,1477,1479,1487,1640,1227,1803),
  groups=c("superfam01-1","superfam01-2","superfam01-3","superfam01-4","superfam01-5","superfam02-1","superfam02-2","superfam02-3","superfam03-1","superfam03-2","superfam03-3","superfam03-4","superfam03-5",
  "extra group","extra group","extra group","extra group","extra group","extra group","extra group","extra group","extra group","extra group","extra group"),
  custom.fill=c("#8dd3c7","#ededa8","#adabc4","#33a02c","#fb8072","#80b1d3","#b3de69","#b0b0b0","#fdb45a","#fccde5","#ffed6f","#bc80bd","#ccebc5","black","black","black","black","black","black","black","black","black","black","black"))

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

tree_dat<-p$data %>% filter(parent %in% families) %>%
    left_join(group_df) %>%
    group_by(parent) %>%
  mutate(parentlab=ifelse(parent%in%families,parent,NA)) %>% 
  ungroup()

# hilights start at layer #3
highlighted_tree_dat<-map(3:(length(families)+2),~{layer_data(p1,.x)}) %>% list_rbind()

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
    size=1,
    onclick = sprintf('%s"%s%s"%s', !!window_open, groups, !!html_extension, !!closing_quote)),
    alpha=0.5) + 
    scale_fill_manual(values=group_df$custom.fill)+
    guides(size="none",fill="none")
    # scale_fill_manual(values=group_df$custom.fill, guide=guide_legend_interactive(title="Starship Families"))
    # theme(legend.position="right")
    # geom_text(data=tree_dat,aes(x = x, y = y, label = groups), size = 5,nudge_x=5)

interactive_plot  <- girafe(
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

saveRDS(interactive_plot,file="RDS/captain-tree.RDS")

# TODO: to make it easier, use the `SharedData` object to subset other data to make plots/tables
# Wrap data frame in SharedData
sd <- SharedData$new(tree_dat)

# Use SharedData like a dataframe with Crosstalk-enabled widgets
row1 <- div(h2("Select a superfamily to subset:"),filter_checkbox("groups","", sd, ~groups,inline=TRUE))

row2 <- div(
  class = "row",
  div(
    class = "col-md-6",
      girafe(
        ggobj =   p + geom_treescale(x=5,y=2,offset=5)+
        geom_rect_interactive(data=sd$data(),
          aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax,fill=groups,
          tooltip = groups,
          data_id = groups,
          size=1,
          onclick = sprintf('%s"%s%s"%s', !!window_open, groups, !!html_extension, !!closing_quote)),
          alpha=0.5) + 
          scale_fill_manual(values=group_df$custom.fill)+
          guides(size="none",fill="none"),
        width_svg = 10,
        height_svg = 14,
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
  ),
  div(
    class = "col-md-6",
        datatable(sd,#selection=list(mode = 'multiple', selected = c(4,10), target ='column'),
        extensions="Scroller", style="bootstrap", class="compact", width="100%",options=list(deferRender=TRUE, scrollY=300, scroller=TRUE)
              )

    # datatable(sd$data() %>% 
    #             mutate(label=str_split(label,"_",simplify=TRUE)[1]) %>%
    #             left_join(tax,by=c("label"="code")) %>% 
    #             select(label,species), 
    #           extensions="Scroller", style="bootstrap", class="compact", width="100%",options=list(deferRender=TRUE, scrollY=300, scroller=TRUE) # TODO: instead use `selection` here
    #           )
  )
)

comb_widget<-bscols(div(class="container",row1,row2)) %>% 
  htmltools::renderTags() 

# save as RDS
bscols(div(class="container",row1,row2)) %>% 
  saveRDS("RDS/captain-tree-wtable.RDS")

# Save the interactive plot as an HTML file
save_html(html=comb_widget,file="/home/adrian/Systematics/Starship_Database/sequenceserver/public/trees/captain-tree.html")

# replace headers
# headers<-read_tsv("/home/adrian/Systematics/Starship_Database/tree/headers.tsv",col_names=NULL) %>% 
#   rename("label"="X3") %>%
#   mutate(X1=stringi::stri_extract_first_regex(X1,"^([^_]*)"),
#   label=gsub("_R_|_\\d+$","",label))

# tree_dat<-data.frame(label=tree$tip.label) %>%
#   left_join(headers) %>%
#   mutate(tip_label=paste0(tolower(X1),"_",label)) %>%
#   relocate("label","X2","X1")

# replace missing labels
# p2$data<-p2$data %>% 
#   filter(!is.na(x) & !is.na(y)) %>%
#   left_join(group_df) %>% 
#   mutate(label=ifelse(is.na(label),groups,label),
#          label=ifelse(label=="extra group",node,label),
#          isTip=ifelse(!is.na(label),TRUE,FALSE)) %>% 
#          select(-groups)

# make subtrees
for(i in families){
  # Wrap data frame in SharedData
  sub_sd <- SharedData$new(tree_dat %>% filter(groups==i))

  # TODO: borrow elements from sequenceserver html and recreate here with `htmltools`
  comb_sub_widget<-bscols(
    div(class="container",
      h1(paste0("Captain gene phylogeny of ",group_df[group_df$parent==i,]$groups)),
      div(class = "col-md-6",
        girafe(
          ggobj = ggtree(tree_subset(tree,node=i,levels_back=0))+
          geom_tiplab(align=T,size=2)+
          geom_point_interactive()+
          geom_treescale()+
          labs(main=paste0("Family ",group_df[group_df$parent==i,]$groups," Phylogeny"))+
          hexpand(0.25),
          width_svg = 14,
          height_svg = 14,
          options = list(
            opts_sizing(width = 0.8),
            opts_zoom(min=0.7,max=10)
          )
        )


        ),
    div(class = "col-md-6",
          datatable(sd,extensions="Scroller", style="bootstrap", class="compact", width="100%",options=list(deferRender=TRUE, scrollY=300, scroller=TRUE))
        )
    )
  ) %>%
  htmltools::renderTags() 

  # Save the interactive plot as an HTML file
  save_html(html=comb_sub_widget,file=paste0("/home/adrian/Systematics/Starship_Database/sequenceserver/public/trees/",group_df[group_df$parent==i,]$groups,".html"))

  # htmlwidgets::saveWidget(interactive_sub_plot, 
  #   file=paste0("/home/adrian/Systematics/Starship_Database/sequenceserver/public/trees/",group_df[group_df$parent==i,]$groups,".html"),
  #   title=paste0("Phylogeny of ",group_df[group_df$parent==i,]$groups," Captain Genes"))

}

#############
# subsetted tree
#############

sub.p<-ggtree(tree_subset(tree,node=families,levels_back=0))+
  geom_tiplab(align=T,size=2)+
  geom_treescale()+
# list tips to keep
keep.tips<-p2$data %>% filter(!is.na(x) & !is.na(y) & isTip == "TRUE") %>% pull(node)

# subsetted tree
ggtree(ape::keep.tip(tip=keep.tips, phy=tree))

# Create your ggtree plot
p3<-p2 %<+% (group_df %>% mutate(groups=ifelse(groups=="extra group",node,groups)) %>% select(groups)) +
  geom_point()+
  geom_tiplab()+
  geom_treescale()+
  geom_point_interactive(data=group_df,aes(tooltip = groups,
    data_id = groups,
    size=1,
    color=groups,
    onclick = sprintf("window.open(\"", node,".html\"")
    ))

# create an interactive visualization of a ggtree with girafe
# Convert ggtree plot to an interactive girafe plot
# Customize interactivity options
interactive_plot  <- girafe(
  ggobj = p3, width_svg = 15,
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
htmlwidgets::saveWidget(interactive_plot, "ggiraph.html",title="Starship Database Phylogeny")
