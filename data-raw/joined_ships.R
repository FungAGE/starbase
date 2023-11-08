library(vroom)
library(tidyverse)
library(RSQLite)
library(pool)

# populate input for ships
# BUG: there's one entry with 2 sequences, using distinct to deal with this for now
joined_ships <- vroom("MTDB/joined_ships.tsv")

# for dataspice
joined_ships %>%
  # TODO: replace this with row_id = UUIDgenerate() ?
  mutate(row_id=NA) %>%
  write_csv("data/csv/joined_ships.csv")

# for package data
usethis::use_data(joined_ships, overwrite = TRUE)

# connect to database file
mydb <- dbConnect(RSQLite::SQLite(), "SQL/starbase.sqlite")

row_id = UUIDgenerate()

# create tables
dbWriteTable(mydb, "ship", joined_ships)

# test
# dbGetQuery(mydb, 'SELECT * FROM starbase LIMIT 5')

# disconnect from database
dbDisconnect(mydb)