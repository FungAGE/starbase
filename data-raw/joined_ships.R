library(vroom)
# populate input for ships
# BUG: there's one entry with 2 sequences, using distinct to deal with this for now
joined_ships<-vroom("MTDB/joined_ships.tsv")

usethis::use_data(joined_ships, overwrite = TRUE)
