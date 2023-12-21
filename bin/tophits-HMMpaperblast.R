library(tidyverse)

files<-Sys.glob("SQL/data/PaperBLAST/*.tsv")

captain_files <- files[grepl("YR",files)]
captain_files <- files[!grepl("tophits",captain_files)]

captain_results<-map_df(captain_files,~{
  read_tsv(.x, show_col_types = FALSE) %>%
    mutate(superfamily=gsub("SQL/data/PaperBLAST/YRsuperfams.|.tsv","",.x)) %>%
    # separate(evidence, into=c("pident","coverage","bitscore"), sep=",") %>%
    # mutate(across(c(pident,coverage,bitscore),~trimws(.))) %>%
    rowwise() %>%
    mutate(palign=str_extract(evidence, "\\(\\d+\\.\\d+%\\)"),
      pcoverage=str_extract(evidence, "\\b\\d+\\.\\d+%\\b"),
      superfam=str_extract(evidence, "\\bsuperfam\\d+(?:-\\d+)?\\b"),
      bitscore=as.numeric(str_extract(evidence, "\\b\\d+\\.\\d+\\b"))
      )
})

top_captain_results<-captain_results %>% group_by(gene_id) %>% slice_max(bitscore)

write_tsv(top_captain_results,"SQL/data/PaperBLAST/YRsuperfams.tophits.tsv")
