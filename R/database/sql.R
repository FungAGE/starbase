library(tidyverse)
library(DBI)

# potential structure:
# - main table connecting starships with cargo
#   - contains what metadata?
# - separate tables for ship sequences (also separate for nucl/prot?) with checksums
#   - contains all sequences, all versions, etc.

mydb <- dbConnect(RSQLite::SQLite(), "/home/adrian/Systematics/Starship_Database/SQL/starbase.sqlite")

df<-read_tsv("/home/adrian/Systematics/Starship_Database/MTDB/starships.db",col_names = c("ome","genus","species","strain","taxonomy","version","source","biosample","assembly_acc","acquisition_date","published","fna","faa","gff3"))

# TODO: read in checksums

dbWriteTable(mydb, "starbase", df)

# dbGetQuery(mydb, 'SELECT * FROM starbase LIMIT 5')

dbDisconnect(mydb)
