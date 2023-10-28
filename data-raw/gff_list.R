library(tidyverse)
gff_list<-ships_with_anno %>%
  group_by(ship_code) %>%
  summarise(named_vec = list(gff3)) %>%
  deframe()

usethis::use_data(gff_list, overwrite = TRUE)
