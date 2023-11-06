library(tidyverse)
fna_list <- ships_with_anno %>%
  group_by(ship_code) %>%
  summarise(named_vec = list(fna)) %>%
  deframe()

usethis::use_data(fna_list, overwrite = TRUE)
