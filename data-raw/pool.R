library(pool)
library(RSQLite)

# populate input for ships
# BUG: there's one entry with 2 sequences, using distinct to deal with this for now
# joined_ships <- vroom("MTDB/joined_ships.tsv",show_col_types = FALSE)

# connect to database file
pool <- pool::dbPool(RSQLite::SQLite(), dbname = "SQL/starbase.sqlite")
con <- pool::poolCheckout(pool)
poolReturn(con)

usethis::use_data(con, overwrite = TRUE)

# create tables
# dbWriteTable(con, "ship", joined_ships)

# disconnect from database
dbDisconnect(con)
