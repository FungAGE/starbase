con <- load_sql()
ship_seqs<-tbl(con, "genome_genome")

# for package data
usethis::use_data(ship_seqs, overwrite = TRUE)

onStop(stop_sql(con))