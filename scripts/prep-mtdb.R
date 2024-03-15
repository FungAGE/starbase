library(tidyverse)
library(DBI)

# mycodb
mycodb_df_extra<-read_tsv("MTDB/mycodb2899.ome2species.txt",col_names=c("ome","taxa"),na=c("","NA","Unknown","unknown"),show_col_types = FALSE) %>%
  mutate(species=trimws(gsub("^([^ ]+)","",taxa)),
    genus=gsub(" (.+)","",taxa)) %>%
  select(-taxa)

mycodb_df<-read_tsv("MTDB/20211215.pub.dedup.db",c('ome', 'genus', 'species', 'strain', 'version', 'biosample','fna', 'faa', 'gff3','taxonomy','missing1', 'missing2','genomeSource','published','assembly_acc','acquisition_date'),na=c("","NA","Unknown","unknown"),show_col_types = FALSE) %>% 
  separate_rows(taxonomy,sep=",") %>%
  mutate(taxonomy=gsub(c("\\{|\\}|\\'"),"",trimws(taxonomy))) %>%
  separate(taxonomy,sep=": ",into=c("rank","name")) %>%
  mutate(name=ifelse(name %in% c("","unknown","Unknown"),NA,name)) %>%
  pivot_wider(id_cols=c("ome","genus","species","strain","version","biosample","assembly_acc","acquisition_date","genomeSource","published","fna","faa","gff3"),names_from="rank",values_from="name") %>%
  full_join(mycodb_df_extra) %>%
  mutate(evidence="starfish",source="mycodb")
  
# fill in metadata using starfish features
mycodb_feat <- read_tsv("ships/starfish/output/mycodb.final.starships.feat",show_col_types = FALSE) %>%
  rename("size"="elementLength") %>%
  mutate(dr=ifelse(!is.na(upDR) & !is.na(downDR),paste0(upDR,"/",downDR),
                   ifelse(is.na(upDR),downDR,
                          ifelse(is.na(downDR),upDR,
                                 ifelse(is.na(upDR) & is.na(downDR),NA)))),
         atir=ifelse(!is.na(upTIR) & !is.na(downTIR),paste0(upTIR,"/",downTIR),
                     ifelse(is.na(upTIR),downTIR,
                            ifelse(is.na(downTIR),upTIR,
                                   ifelse(is.na(upTIR) & is.na(downTIR),NA)))))

# fill in lineage information
mycodb_lineage <- read_tsv("ships/starfish/output/mycodb.final.lineage.txt",col_names = c("starshipID","lineages"),show_col_types = FALSE) %>% 
  separate(lineages,into = c("starship_family","starship_navis","starship_haplotype"),sep = ":")

# manual
# TODO: add publication information
manual_df<-read_csv("ships/manual-annotations/Starships.fulltaxa.csv",na=c("","NA","Unknown","unknown","?"),show_col_types = FALSE) %>%
  rename("starship_navis"="Navis","starship_family"="Family") %>%
  select(-species) %>%
  rename_all(str_to_lower) %>%
  mutate(source="manual",evidence="manual",genomeSource="ncbi")

# ships
mycodb_ships <- left_join(
  tibble(fna=Sys.glob("SQL/data/fna/ships/starfish/*")) %>% 
    rowwise() %>%
    mutate(starshipID=gsub("mycodb.final.starships.part_|__.*|.fna|.fa","",basename(fna)),
          ome=str_split(starshipID,"_",simplify=TRUE)[1]),
  tibble(gff3=Sys.glob("metadata/ships/starfish/gff/starfish/*")) %>% 
    rowwise() %>%
    mutate(starshipID=gsub(".gff","",basename(gff3)),
          ome=str_split(starshipID,"_",simplify=TRUE)[1])
          ,by=c("starshipID","ome")) %>% 
  left_join(mycodb_df %>% select(-c('fna', 'faa', 'gff3'))) %>% 
  left_join(mycodb_feat) %>%
  left_join(mycodb_lineage) %>%
  filter(!is.na(genus))

manual_ships<-manual_df %>%
  select(-taxrank) %>%
  rowwise() %>%
  mutate(fna=ifelse(length(Sys.glob(paste0("SQL/data/fna/ships/manual-annotations/*",code,".*")))==1,
                    Sys.glob(paste0("SQL/data/fna/ships/manual-annotations/*",code,".*")),NA),
        gff3=NA) %>%
  ungroup()

# join ships and checksums
# also adding ome's to manual entries, like in mycotools: https://github.com/xonq/mycotools/blob/0b4ff2977f291e7637363c0e5592d3d73a138e67/mycotools/predb2mtdb.py#L191
# first three letters of genus, first three letters of species (or "sp."), unique database number, and optional MTDB version tag
# not adding version for now...
ship_checksums<- read_delim("SQL/data/fna/ships/fna.checksums.txt",col_names=c("checksum","fna"),show_col_types = FALSE) %>%
  mutate(fna=gsub("./","SQL/data/fna/ships/",fixed=TRUE,fna),
         checksum=gsub("seqkit.v0.1_DLS_k0_","",checksum))

joined_ships <- full_join(mycodb_ships,manual_ships) %>%
  left_join(ship_checksums) %>% 
  group_by(genus) %>%
  fill(c(superkingdom, clade, kingdom, subkingdom, phylum, subphylum, class, order, family,genus)) %>% # fill in missing taxonomic info
  mutate(ome=ifelse(is.na(ome),
                      paste0(str_to_lower(str_sub(string=str_extract(species,"^([^ ]+)"),start=1,end=3)),
                        ifelse(is.na(genus) | genus == "" | str_length(genus)<=3,"sp.",str_to_lower(str_sub(genus,start=1,end=3))),
                        row_number()+sum(ifelse(source=="mycodb", 1, 0))),ome)) %>%
  ungroup() %>%
  # * group_by excludes rows where the grouping column is NA
  # ? if the checksums are identical, we would assume that these should all be the same anyway?
  # * found little/no overlap with checksums...
  # fill(c(biosample,assembly_acc,acquisition_date,genomeSource,published),.direction="downup")
  group_by(checksum) %>%
  fill(colnames(mycodb_feat)) %>%
  # TODO: overwrite mycodb meta with manual annotations meta, if it exists
  fill(c(size, dr, atir, target, spok, ars, other, hgt))  %>%
  group_by(starship_haplotype) %>%
  fill(starship_family,starship_navis) %>%
  group_by(starship_navis) %>%
  fill(starship_family) %>%
  ungroup() %>%
  # ! temporary adding `starshipID` for manual annotations
  mutate(starshipID=ifelse(is.na(starshipID), ifelse(!is.na(checksum),paste0(ome,"_s",str_sub(checksum,start=1,end=5)),ome),starshipID),
         assembly_acc=ifelse(is.na(assembly_acc),ome,assembly_acc)) %>% 
  rename("tir"="atir")

# ? add version = 1 where missing?

# save table for now
joined_ships %>% 
  write_tsv("MTDB/joined_ships.tsv")

# how many present/missing
joined_ships %>% group_by(source) %>%
  summarise(fna=sum(!is.na(fna)),gff3=sum(!is.na(gff3)),checksums=sum(!is.na(checksum)))

# how many duplicates?
joined_ships %>% summarise(uniq=n_distinct(checksum),dups=nrow(.)-uniq)

# write predb file
ship_predb<-joined_ships %>%
  mutate(useRestriction="no", # assuming this should be ok for our purposes
    genomeSource=ifelse(is.na(genomeSource)| genomeSource == "slot","ncbi",genomeSource)) %>%
  select(assembly_acc,ome,genus,species,strain,version,biosample,fna,gff3,genomeSource,useRestriction,published)

# if fna/gff3 is missing, then prepdb2mtdb will fail
ship_predb %>%
  filter(!is.na(fna) & !is.na(gff3)) %>%
  write_tsv("MTDB/ships_predb.tsv",col_names=FALSE,na="")

# write mtdb file
# TODO: turn paths into absolute paths
joined_ships %>%
  filter(!is.na(fna) & !is.na(gff3)) %>%
  mutate(useRestriction="no", # assuming this should be ok for our purposes
    genomeSource=ifelse(is.na(genomeSource)| genomeSource == "slot","ncbi",genomeSource),
    faa=NA) %>%
  select(ome,genus,species,strain,version,source,biosample,assembly_acc,published,acquisition_date,fna,faa,gff3) %>%
  left_join(read_tsv("MTDB/20211215.pub.dedup.db",c('ome', 'genus', 'species', 'strain', 'version', 'biosample','fna', 'faa', 'gff3','taxonomy','missing1', 'missing2','genomeSource','published','assembly_acc','acquisition_date'),na=c("","NA","Unknown","unknown"),show_col_types = FALSE) %>% select(ome,taxonomy,assembly_acc)) %>%
  select(ome,genus,species,strain,taxonomy,version,source,biosample,assembly_acc,published,acquisition_date,fna,faa,gff3) %>%
  filter(!is.na(taxonomy)) %>% 
  write_tsv(file="MTDB/ships.mtdb",col_names=NULL)

# captain genes
mycodb_captain <- mycodb_df %>%
  select(-gff3) %>%
  rowwise() %>%
  mutate(faa=ifelse(length(Sys.glob(paste0("SQL/data/faa/captain/tyr/starfish/*",ome,"*")))==1,
              Sys.glob(paste0("SQL/data/faa/captain/tyr/starfish/*",ome,"*")),NA),
        gene="tyr") %>%
  ungroup()

manual_captain<-manual_df %>%
  rowwise() %>%
  mutate(faa=ifelse(length(Sys.glob(paste0("SQL/data/faa/captain/tyr/manual/*",code,".*")))==1,
                    Sys.glob(paste0("SQL/data/faa/captain/tyr/manual/*",code,".*")),NA),
        gene="tyr") %>%
  ungroup()

# load checksums
# join captain genes
captain_checksums<-read_delim("SQL/data/faa/captain/tyr/faa.checksums.txt",col_names=c("checksum","file"),show_col_types = FALSE) %>% 
  mutate(file=basename(file))
joined_captain <- full_join(mycodb_captain,manual_captain) %>%
  mutate(file=basename(faa)) %>% # just in case files have changed locations...
  left_join(captain_checksums) %>%
  # mutate(checksum=gsub("seqkit.v0.1_PLS_k0_","",checksum)) %>%
  group_by(species,genus) %>%
  mutate(ome=ifelse(is.na(ome),
                      paste0(str_to_lower(str_sub(string=str_extract(species,"^([^ ]+)"),start=1,end=3)),
                        ifelse(is.na(genus) | genus == "" | str_length(genus)<=3,"sp.",str_to_lower(str_sub(genus,start=1,end=3))),
                        row_number()+sum(ifelse(source=="mycodb", 1, 0))),ome)) %>%
  ungroup() %>%
  mutate(assembly_acc=ifelse(is.na(assembly_acc),ome,assembly_acc))

joined_captain %>% summarise(faa=sum(!is.na(faa)),checksums=sum(!is.na(checksum)))

joined_captain %>%
  write_tsv(., "MTDB/captain.tsv",col_names=FALSE,na="")

# load cargo genes
# captain genes
mycodb_cargo <- mycodb_df %>%
  select(-faa) %>%
  rowwise() %>%
  mutate(nlr=ifelse(length(Sys.glob(paste0("SQL/data/faa/cargo/nlr/starfish/*",ome,"*")))==1,
              Sys.glob(paste0("SQL/data/faa/cargo/nlr/starfish/*",ome,"*")),NA),
        fre=ifelse(length(Sys.glob(paste0("SQL/data/faa/cargo/fre/starfish/*",ome,"*")))==1,
              Sys.glob(paste0("SQL/data/faa/cargo/fre/starfish/*",ome,"*")),NA),
        plp=ifelse(length(Sys.glob(paste0("SQL/data/faa/cargo/plp/starfish/*",ome,"*")))==1,
              Sys.glob(paste0("SQL/data/faa/cargo/plp/starfish/*",ome,"*")),NA),
        duf3723=ifelse(length(Sys.glob(paste0("SQL/data/faa/cargo/duf3723/starfish/*",ome,"*")))==1,
              Sys.glob(paste0("SQL/data/faa/cargo/duf3723/starfish/*",ome,"*")),NA),) %>%
  ungroup() %>%
  pivot_longer(cols=c(nlr,fre,plp,duf3723),names_to="gene",values_to="faa") %>%
  mutate(file=basename(faa)) # just in case files have changed locations...

# load checksums
faa_files<-c("SQL/data/faa/cargo/nlr/faa.checksums.txt",
  "SQL/data/faa/cargo/fre/faa.checksums.txt",
  "SQL/data/faa/cargo/plp/faa.checksums.txt",
  "SQL/data/faa/cargo/duf3723/faa.checksums.txt")

cargo_checksums<-map_df(faa_files,~{read_delim(.x,col_names=c("checksum","file"),show_col_types = FALSE) %>% 
  mutate(gene=gsub("SQL/data/faa/cargo/","",gsub("/faa.checksums.txt","",.x)),
  # checksum=gsub("seqkit.v0.1_PLS_k0_","",checksum),
  file=basename(file))})

# join cargo genes
joined_cargo <- mycodb_cargo %>%
  left_join(cargo_checksums) %>%
  full_join(manual_ships %>% select(-fna,-gff3)) %>%
  group_by(species,genus) %>%
  mutate(ome=ifelse(is.na(ome),
                      paste0(str_to_lower(str_sub(string=str_extract(species,"^([^ ]+)"),start=1,end=3)),
                        ifelse(is.na(genus) | genus == "" | str_length(genus)<=3,"sp.",str_to_lower(str_sub(genus,start=1,end=3))),
                        row_number()+sum(ifelse(source=="mycodb", 1, 0))),ome)) %>%
  ungroup() %>%
  mutate(assembly_acc=ifelse(is.na(assembly_acc),ome,assembly_acc))

joined_cargo %>% summarise(faa=sum(!is.na(faa)),checksums=sum(!is.na(checksum)))

joined_cargo %>%
  write_tsv(., "MTDB/cargo.tsv",col_names=FALSE,na="")

# genome database
# fill in missing accession numbers for assemblies
library(rentrez)

joined_ships %>% mutate(
  assembly_acc=ifelse(
    is.na(fna) & is.na(assembly_acc),
      entrez_search(db="assembly",
                term=paste0('"',genus," ",species,"[ORGN]",'"'),
                retmax=0)
  )
)
