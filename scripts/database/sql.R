library(tidyverse)
library(DBI)
library(DT)

# potential structure:
# - main table containing metadata
#   - connecting to all other tables
# - separate tables for ship and gene sequences (also separate for nucl/prot?) with checksums
#   - contains all sequences, all versions, etc.

mydb <- dbConnect(RSQLite::SQLite(), "SQL/starbase.sqlite")

# mycodb
# TODO: add ncbi taxids?

mycodb_df<-read_tsv("MTDB/20211215.pub.dedup.db",c('ome', 'genus', 'species', 'strain', 'version', 'biosample','fna', 'faa', 'gff3','taxonomy','missing1', 'missing2','source','published','assembly_acc','acquisition_date'),show_col_types = FALSE) %>% 
  rowwise() %>%
  mutate(fna=ifelse(length(Sys.glob(paste0("SQL/data/fna/ships/",ome,"*.fna")))==1,
                    Sys.glob(paste0("SQL/data/fna/ships/",ome,"*.fna")),NA),
         gff=ifelse(length(Sys.glob(paste0("SQL/data/gff/mycodb/",ome,".gff")))==1,
                    Sys.glob(paste0("SQL/data/gff/mycodb/",ome,".gff")),NA)) %>%
  ungroup() %>%
  separate_rows(taxonomy,sep=",") %>%
  mutate(taxonomy=gsub(c("\\{|\\}|\\'"),"",trimws(taxonomy)),
    evidence="starfish") %>%
  separate(taxonomy,sep=": ",into=c("rank","name")) %>%
  pivot_wider(id_cols=c("ome","genus","species","strain","version","source","biosample","assembly_acc","acquisition_date","published","fna","faa","gff3"),names_from="rank",values_from="name")

# manual
# TODO: add ome's like here: https://github.com/xonq/mycotools/blob/0b4ff2977f291e7637363c0e5592d3d73a138e67/mycotools/predb2mtdb.py#L191
manual_df<-read_csv("Starships/ships/manual-annotations/Starships.fulltaxa.csv") %>%
  rename("starship_class"="Class") %>%
  rename_all(tolower) %>%
  rowwise() %>%
  mutate(fna=ifelse(length(Sys.glob(paste0("Starships/starships/starfish/mycodb.final.starships.fna.split/",ome,"*.fna")))==1,
                    Sys.glob(paste0("Starships/starships/starfish/mycodb.final.starships.fna.split/",ome,"*.fna")),NA),
         gff=ifelse(length(Sys.glob(paste0("MTDB/mycodb.final.starships_split/",ome,".gff")))==1,
                    Sys.glob(paste0("MTDB/mycodb.final.starships_split/",ome,".gff")),NA)) %>%
  ungroup() %>%
  mutate(source="manual",
        evidence="manual")


# ships
## checksums
read_file("SQL/data/fna/ships")

# genes
## checksums
read_file("SQL/data/faa/cargo/tyr")

mycodb_genes<-mycodb_df %>%
  rowwise() %>%
  mutate(tyr_faa=ifelse(length(Sys.glob(paste0("SQL/data/faa/cargo/tyr/mycodb/*",ome,"*.faa")))==1,
                    Sys.glob(paste0("SQL/data/faa/cargo/tyr/mycodb/*",ome,"*.faa")),NA),
         plp_faa=ifelse(length(Sys.glob(paste0("SQL/data/faa/cargo/plp/mycodb/*",ome,"*.faa")))==1,
                        Sys.glob(paste0("SQL/data/faa/cargo/plp/mycodb/*",ome,"*.faa")),NA),
         nlr_faa=ifelse(length(Sys.glob(paste0("SQL/data/faa/cargo/nlr/mycodb/*",ome,"*.faa")))==1,
                        Sys.glob(paste0("SQL/data/faa/cargo/nlr/mycodb/*",ome,"*.faa")),NA),
         fre_faa=ifelse(length(Sys.glob(paste0("SQL/data/faa/cargo/fre/mycodb/*",ome,"*.faa")))==1,
                        Sys.glob(paste0("SQL/data/faa/cargo/fre/mycodb/*",ome,"*.faa")),NA),
         duf3723_faa=ifelse(length(Sys.glob(paste0("SQL/data/faa/cargo/duf3723/mycodb/*",ome,"*.faa")))==1,
                        Sys.glob(paste0("SQL/data/faa/cargo/duf3723/mycodb/*",ome,"*.faa")),NA)) %>%
  ungroup()

# dbWriteTable(mydb, "starbase", df)

# TODO: read in checksums

dbGetQuery(mydb, 'SELECT * FROM starbase') %>%
  distinct(ome,genus,species,strain) %>%
  datatable(options = list(), class = "display",
            callback = JS("return table;"), #rownames, colnames, container,
            caption = NULL, filter = c("none", "bottom", "top"), escape = TRUE,
            style = "auto", width = NULL, height = NULL, elementId = NULL,
            fillContainer = getOption("DT.fillContainer", NULL),
            autoHideNavigation = getOption("DT.autoHideNavigation", NULL),
            selection = c("multiple", "single", "none"), extensions = list(),
            plugins = NULL, editable = FALSE)

dbDisconnect(mydb)
