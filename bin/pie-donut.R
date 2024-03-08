library(tidyverse)
library(webr)
library(plotly)

load("data/joined_ships.rda")

dat<-joined_ships %>%
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

plot_dat <- dat %>% 
  filter(class != "Unknown") %>%
  group_by(class,order) %>% 
  summarise(n=n_distinct(ome))
p<-PieDonut(plot_dat,aes(class, order, count=n))

svg("inst/app/img/pie-donut.svg")
p
dev.off()

htmlwidgets::saveWidget(ggplotly(p), "inst/app/www/html/pie-donut.html",title="Starship Database Diversity")

load("data/joined_ships.rda")
plot_dat<-joined_ships %>% 
  filter(!is.na(order) & !is.na(starship_family)) %>%
  group_by(starship_family,order) %>%
  summarise(n=n_distinct(ome))

families<-read_tsv("/home/adrian/Systematics/Starship_Database/Starships/metadata/family/family-names.tsv")

# check
# plot_dat %>%
#   filter(!starship_family %in% c(families$otherFamilyID,families$familyName))

plot_dat_sub<-bind_rows(joined_ships %>% inner_join(families,by=c("starship_family"="otherFamilyID")),
    joined_ships %>% inner_join(families,by=c("starship_family"="familyName"))) %>%
  filter(!is.na(order) & !is.na(starship_family)) %>%
  group_by(longFamilyID,order) %>%
  summarise(n=n_distinct(ome))

svg("inst/app/img/pie-donut-family.svg")
PieDonut(plot_dat_sub,aes(longFamilyID,order, count=n))
dev.off()