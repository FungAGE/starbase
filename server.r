library(tidyverse)
library(reticulate)
library(XML)
library(DT)

server <- function(input, output, session){

  # example
  # dat<-read_csv("tmp/hmmer-parsed.csv",show_col_types = FALSE) %>% select(-1)
  # print(class(dat))
  # dat
  # output$tbl = renderDT(dat)
    
  blastresults <- reactive({
    # separate lists of starship and gene databases with sub lists of nucl and prot databases
    # TODO: add starship protein databases and gene hmm to folder
    db_list <- list(starship=list(curated =
                                    list(nucl = "/home/adrian/Systematics/Starship_Database/blastdb/concatenated.fa",
                                         prot = "/home/adrian/Systematics/Starship_Database/blastdb/concatenated.faa"),
                                  starfish =
                                    list(nucl = "/home/adrian/Systematics/Starship_Database/blastdb/mycodb.final.starships.headers.fna",
                                         prot = "/home/adrian/Systematics/Starship_Database/blastdb/mycodb.final.starships.headers.faa")),
                    gene=list(tyr = list(prot = "/home/adrian/Systematics/Starship_Database/blastdb/YRsuperfamRefs.faa",
                                         nucl = ""),
                              freB = list(prot = "",
                                          nucl = ""),
                              nlr = list(prot = "",
                                         nucl = ""),
                              DUF3723 = list(prot = "",
                                             nucl = ""),
                              plp = list(prot = "",
                                         nucl = "")))

    # gather input and set up temp file
    query <- input$query
    tmp_fasta <- "tmp/fasta.fa"

    # this makes sure the fasta is formatted properly
    # removes whitespace at the beginning of lines before testing for header
    if (!startsWith(gsub("^\\s+", "", query, perl = TRUE), ">")){
      query<-cat(">Query",query,sep = "\n")
    }

    # function to clean the query sequence: first remove non-letter characters, then ambiguous characters, then guess if query is nucl or protein
    clean_lines <- function(lines) {
      cleaned_lines <- NULL
      na_count <- 0
      for (line in lines) {
        if (grepl("^>",trimws(line))) {
          header <- line
          # header<-str_split(line,"\n",simplify=TRUE)[1]
          cleaned_lines <- str_c(cleaned_lines, header,"\n")
        } else {
          # Remove non-letter characters using gsub
          cleaned_line <- gsub("[NXnx]", "", gsub("[^A-Za-z]", "", line))
          na_count<-na_count+sum(table(unlist(strsplit(cleaned_line, "")))[c("N","X","n","x")],na.rm=TRUE)
          cleaned_lines <- str_c(cleaned_lines, cleaned_line)
        }
      }
      cleaned_seq<-str_split(cleaned_lines,"\n",simplify=TRUE)[2]

      # count characters for deciding type later
      nucl_count<-sum(table(unlist(strsplit(cleaned_seq, "")))[c("a","A","t","T","g","G","c","C")],na.rm=TRUE)
      prot_count<-sum(table(unlist(strsplit(cleaned_seq, "")))[c(letters,LETTERS)],na.rm=TRUE)

      # guess nucl or protein
      if(prot_count >= (0.9 * str_length(cleaned_seq)+na_count)){
        query_type<-"prot"
        print("Query is protein sequence")
      } else {
        query_type<-"nucl"
        print("Query is nucleotide sequence")
      }
      writeLines(cleaned_lines, tmp_fasta)
      return(list("header"=header,"query_type"=query_type))
    }

    query_list <- clean_lines(query)

    # force the right database to be used
    # and calls the BLAST/HMMER
    db<-db_list[[input$input_type]][[input$db_type]][[query_list$query_type]]

    # select correct BLAST/HMMER program
    if(input$input_type == "gene" | input$search_ship_genes == "Yes") {
      if(query_list$query_type=="nucl") {
        hmmer_program <- "jackhmmer"
      } else {
        hmmer_program <- "phmmer"
      }
      if(input$input_type == "starship" & input$search_ship_genes == "Yes" ) {
        gene_list <- names(db_list[["gene"]])
      } else {
        gene_list <- input$db_type
      }
      hmmer_list<-NULL
      for(gene in gene_list) {
        hmmer_cmd <- paste(hmmer_program,"--domtblout tmp/hmmer.out","--cpu 4","--domE",input$eval,tmp_fasta,db,sep=" ")
        system(hmmer_cmd, intern = FALSE)
      }
      # TODO: combine results from multiple genes into list before parser
      system("python bin/hmm.py -i tmp/hmmer.out")
      read_csv("tmp/hmmer-parsed.csv",show_col_types = FALSE) %>% select(-1)
    } else if(input$input_type == "starship" & input$search_ship_genes == "No" ) {
      if(query_list$query_type=="nucl") {
        blast_program <- "blastn"
      } else {
        blast_program <- "blastp"
      }
      # TODO: should be separated by results from different genes?
      # TODO: split with conditional statements for HMMER, BLAST, or both
      system(paste0(blast_program," -query ",tmp_fasta," -db ",db," -evalue ",input$eval," -outfmt 5 -max_hsps 1 -max_target_seqs 10 -num_threads 4"), intern = T) %>%
        xmlParse() %>%
        xpathApply(.,'//Iteration',function(row){
          query_ID <- getNodeSet(row, 'Iteration_query-def') %>% sapply(., xmlValue)
          hit_IDs <- getNodeSet(row, 'Iteration_hits//Hit//Hit_id') %>% sapply(., xmlValue)
          hit_length <- getNodeSet(row, 'Iteration_hits//Hit//Hit_len') %>% sapply(., xmlValue)
          bitscore <- getNodeSet(row, 'Iteration_hits//Hit//Hit_hsps//Hsp//Hsp_bit-score') %>% sapply(., xmlValue)
          eval <- getNodeSet(row, 'Iteration_hits//Hit//Hit_hsps//Hsp//Hsp_evalue') %>% sapply(., xmlValue)
          bind_rows(query_ID,hit_IDs,hit_length,bitscore,eval)
        }) %>%
        lapply(.,function(y){as.data.frame((y),stringsAsFactors=FALSE)}) %>%
        bind_rows()
      }
    }) %>%
    bindEvent(input$blast)
  
  output$tbl = renderDT(blastresults()) 
  
  #this chunk gets the alignemnt information from a clicked row
  # output$clicked <- renderTable({
  #   if(is.null(input$blastResults_rows_selected)){}
  #   else{
  #     xmltop = xmlRoot(blastresults())
  #     clicked = input$blastResults_rows_selected
  #     tableout<- data.frame(parsedresults()[clicked,])
  #     
  #     tableout <- t(tableout)
  #     names(tableout) <- c("")
  #     rownames(tableout) <- c("Query ID","Hit ID", "Length", "Bit Score", "e-value")
  #     colnames(tableout) <- NULL
  #     data.frame(tableout)
  #   }
  # },rownames =T,colnames =F)
  # 
  # #this chunk makes the alignments for clicked rows
  # output$alignment <- renderText({
  #   if(is.null(input$blastResults_rows_selected)){}
  #   else{
  #     xmltop = xmlRoot(blastresults())
  #     
  #     clicked = input$blastResults_rows_selected
  #     
  #     #loop over the xml to get the alignments
  #     align <- xpathApply(blastresults(), '//Iteration',function(row){
  #       top <- getNodeSet(row, 'Iteration_hits//Hit//Hit_hsps//Hsp//Hsp_qseq') %>% sapply(., xmlValue)
  #       mid <- getNodeSet(row, 'Iteration_hits//Hit//Hit_hsps//Hsp//Hsp_midline') %>% sapply(., xmlValue)
  #       bottom <- getNodeSet(row, 'Iteration_hits//Hit//Hit_hsps//Hsp//Hsp_hseq') %>% sapply(., xmlValue)
  #       rbind(top,mid,bottom)
  #     })
  #     
  #     #split the alignments every 40 carachters to get a "wrapped look"
  #     alignx <- do.call("cbind", align)
  #     splits <- strsplit(gsub("(.{40})", "\\1,", alignx[1:3,clicked]),",")
  #     
  #     #paste them together with returns '\n' on the breaks
  #     split_out <- lapply(1:length(splits[[1]]),function(i){
  #       rbind(paste0("Q-",splits[[1]][i],"\n"),paste0("M-",splits[[2]][i],"\n"),paste0("H-",splits[[3]][i],"\n"))
  #     })
  #     unlist(split_out)
  #   }
  # })
}