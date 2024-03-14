library(tidyverse)

load("data/joined_ships.rda")

# Define the metadata
metadata<-joined_ships %>%
  filter(!is.na(starship_family)) %>%
  mutate(starship_family=gsub("fam","superfam0",starship_family),
    starship_family = ifelse(grepl("^fam", starship_family) & !is.na(code), code, starship_family))

# for package data
usethis::use_data(metadata, overwrite = TRUE)
