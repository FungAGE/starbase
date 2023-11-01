library(vroom)
library(tidyverse)
# populate input for ships
# BUG: there's one entry with 2 sequences, using distinct to deal with this for now
joined_ships<-vroom("MTDB/joined_ships.tsv")

# for dataspice
joined_ships %>%
  write_csv("data/csv/joined_ships.csv")

# for package data
usethis::use_data(joined_ships, overwrite = TRUE)
