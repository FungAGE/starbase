library(vroom)
library(tidyverse)
library(RSQLite)

joined_ships <- dbGetQuery(con, 'SELECT * FROM joined_ships')

# for package data
usethis::use_data(joined_ships, overwrite = TRUE)
