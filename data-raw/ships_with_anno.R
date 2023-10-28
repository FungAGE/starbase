library(tidyverse)
ships_with_anno<-joined_ships %>%
  filter(!is.na(fna) & !is.na(gff3)) %>% 
  mutate(ship_code=ifelse(is.na(ship_code),
    paste0(ome,"_s",str_sub(checksum,start=1,end=5)),ship_code)) %>% 
  distinct(ship_code,fna,gff3,genus,species) 

usethis::use_data(ships_with_anno, overwrite = TRUE)
