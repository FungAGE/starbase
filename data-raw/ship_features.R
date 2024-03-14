con <- load_sql()
ship_features<-tbl(con, "genome_features")

# TODO: reduce the size of the object that gets saved

# for package data
usethis::use_data(ship_features, overwrite = TRUE)

onStop(stop_sql(con))