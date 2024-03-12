library(tidyverse)

load("data/ship_df.rda")
ship_ids <- ship_df %>% distinct(genome_name) %>% pull(genome_name)

# for package data
usethis::use_data(ship_ids, overwrite = TRUE)

