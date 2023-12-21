library(tidyverse)
load("data/joined_ships.rda")

ships_with_anno <- joined_ships %>%
  filter(!is.na(fna) & !is.na(gff3)) %>%
  mutate(starshipID = ifelse(is.na(starshipID),
    paste0(ome, "_s", str_sub(checksum, start = 1, end = 5)), starshipID
  )) %>%
  distinct(starshipID, fna, gff3, genus, species)

usethis::use_data(ships_with_anno, overwrite = TRUE)

# for dataspice
ships_with_anno %>%
  write_csv("data/csv/ships_with_anno.csv")
