library(vroom)
library(tidyverse)

con<-load_sql()
joined_ships <- dbGetQuery(con, 'SELECT * FROM joined_ships')

# for package data
usethis::use_data(joined_ships, overwrite = TRUE)

onStop(stop_sql(con))