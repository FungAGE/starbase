library(tidyverse)
taxa_list<-ships_with_anno %>%
  mutate(taxa=paste(genus,species)) %>%
  group_by(ship_code) %>%
  summarise(named_vec = list(taxa)) %>%
  deframe()

usethis::use_data(taxa_list, overwrite = TRUE)
