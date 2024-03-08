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

# TODO: to make it easier, use the `SharedData` object to subset other data to make plots/tables
# Wrap data frame in SharedData
sd <- SharedData$new(tree_dat)

# Use SharedData like a dataframe with Crosstalk-enabled widgets
row1 <- div(h2("Select a family to subset:"),
            datatable(sd$data() %>% distinct(groups,clade,newFamilyID,familyName,`type element reference`,notes))
            # filter_checkbox("groups","", sd, ~groups,inline=TRUE)
            )

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
          size=1
          # onclick = sprintf('%s"%s%s"%s', !!window_open, groups, !!html_extension, !!closing_quote)
          ),
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
        datatable(sd %>% dplyr::select(),#selection=list(mode = 'multiple', selected = c(4,10), target ='column'),
        extensions="Scroller", style="bootstrap", class="compact", width="100%",options=list(deferRender=TRUE, scrollY=300, scroller=TRUE)
              )

    # datatable(sd$data() %>% 
    #             mutate(label=str_split(label,"_",simplify=TRUE)[1]) %>%
    #             left_join(tax,by=c("label"="code")) %>% 
    #             dplyr::select(label,species), 
    #           extensions="Scroller", style="bootstrap", class="compact", width="100%",options=list(deferRender=TRUE, scrollY=300, scroller=TRUE) # TODO: instead use `selection` here
    #           )
  )
)

comb_widget<-bscols(div(class="container",row1,row2)) %>% 
  htmltools::renderTags() 

# save as RDS
bscols(div(class="container",row1,row2)) %>% 
  saveRDS("data/captain-tree-wtable.RDS")

# Save the interactive plot as an HTML file
save_html(html=comb_widget,file="sequenceserver/public/trees/captain-tree.html")

# replace headers
# headers<-read_tsv("tree/headers.tsv",col_names=NULL) %>% 
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
#          dplyr::select(-groups)

# make subtrees
for(i in family_nodes){
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
  save_html(html=comb_sub_widget,file=paste0("sequenceserver/public/trees/",group_df[group_df$parent==i,]$groups,".html"))

  # htmlwidgets::saveWidget(interactive_sub_plot, 
  #   file=paste0("sequenceserver/public/trees/",group_df[group_df$parent==i,]$groups,".html"),
  #   title=paste0("Phylogeny of ",group_df[group_df$parent==i,]$groups," Captain Genes"))

}

#############
# subsetted tree
#############

sub.p<-ggtree(tree_subset(tree,node=family_nodes,levels_back=0))+
  geom_tiplab(align=T,size=2)+
  geom_treescale()+
# list tips to keep
keep.tips<-p2$data %>% filter(!is.na(x) & !is.na(y) & isTip == "TRUE") %>% pull(node)

# subsetted tree
ggtree(ape::keep.tip(tip=keep.tips, phy=tree))

# Create your ggtree plot
p3<-p2 %<+% (group_df %>% mutate(groups=ifelse(groups=="extra group",node,groups)) %>% dplyr::select(groups)) +
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
htmlwidgets::saveWidget(interactive_plot, "html/ggiraph.html",title="Starship Database Phylogeny")
