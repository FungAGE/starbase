con <- load_sql()
ship_df<-tbl(con, "genome_genome")

# for package data
usethis::use_data(ship_df, overwrite = TRUE)

onStop(stop_sql(con))