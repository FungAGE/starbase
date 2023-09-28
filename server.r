library(XML)
library(xml2)
library(dplyr)
library(DT)

server <- function(input, output, session){
  blastresults <- reactive({
    # separate lists of starship and gene databases with sub lists of nucl and prot databases
    # TODO: add starship protein databases and gene hmm to folder
    db_list <- list(starship=list(curated = 
                                       list(nucl = "/home/adrian/Systematics/Starship_Database/blastdb/concatenated.fa",
                                            prot = "/home/adrian/Systematics/Starship_Database/blastdb/concatenated.faa"),
                                     starfish = 
                                       list(nucl = "/home/adrian/Systematics/Starship_Database/blastdb/mycodb.final.starships.headers.fna",
                                            prot = "/home/adrian/Systematics/Starship_Database/blastdb/mycodb.final.starships.headers.faa")),
                    gene=list(tyr = list(prot = "/home/adrian/Systematics/Starship_Database/hmm/YRsuperfamRefs.nucl.hmm",
                                              nucl = "/home/adrian/Systematics/Starship_Database/hmm/YRsuperfamRefs.prot"),
                                 freB = list(prot = "/home/adrian/Systematics/Starship_Database/hmm/freB.prot.hmm",
                                              nucl = "/home/adrian/Systematics/Starship_Database/hmm/freB.nucl.hmm"),
                                 nlr = list(prot = "/home/adrian/Systematics/Starship_Database/hmm/nlr.prot.hmm",
                                             nucl = "/home/adrian/Systematics/Starship_Database/hmm/nlr.nucl.hmm"),
                                 DUF3723 = list(prot = "/home/adrian/Systematics/Starship_Database/hmm/DUF3723.prot.hmm",
                                                 nucl = "/home/adrian/Systematics/Starship_Database/hmm/DUF3723.nucl.hmm"),
                                 plp = list(prot = "/home/adrian/Systematics/Starship_Database/hmm/plp.prot.hmm",
                                             nucl = "/home/adrian/Systematics/Starship_Database/hmm/plp.nucl.hmm")))
    
    # gather input and set up temp file
    query <- input$query
    tmp <- tempfile(fileext = ".fa")
    
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
      if(nucl_count >= (0.9 * length(cleaned_seq)+na_count)){
        query_type<-"nucl"
      } else {
        query_type<-"prot"
      }
      writeLines(cleaned_lines, tmp)
      return(list("header"=header,"query_type"=query_type))
    }
    
    query_list <- clean_lines(query)
  
    # force the right database to be used
    # and calls the BLAST/HMMER    
    db<-db_list[[input$input_type]][[input$db_type]][[query_list$query_type]]

    # select correct BLAST/HMMER program
    if(input$input_type == "gene" | input$search_ship_genes == "Yes") {
      if(query_type=="nucl") {
        hmmer_program <- "jackhmmer"
      } else {
        hmmer_program <- "phmmer"
      }
      if(input$input_type == "starship" & input$search_ship_genes == "Yes" ) {
        gene_list <- names(db_list[["gene"]])
      } else {
        gene_list <- input$db_type
      }
      for(gene in gene_list) {
        gene_list <- lapply(names(db_list[["gene"]]),function(x) x[[query_list$query_type]])
        hmmer_cmd <- paste(hmmer_program,"--domtblout --cpu 4",tmp,db,sep=" ")
        # TODO: combine results from multiple genes into list
        system(hmmer_cmd, intern = T)
        # TODO: hmmer output parser
      }
    } else if(input$input_type == "starship" & input$search_ship_genes == "No" ) {
      if(query_type=="nucl") {
        blast_program <- "blastn"
      } else {
        blast_program <- "blastp"
      }
      blast_data<-system(paste0(blast_program," -query ",tmp," -db ",db," -evalue ",input$eval," -outfmt 5 -max_hsps 1 -max_target_seqs 10 -num_threads 4 -out tmp/output.xml"), intern = F)
      print("tmp/output.xml")
    }
  }) %>%
    bindEvent(input$blast)
  
  output$blastResults <- renderUI({
    tags$iframe(src = "inst/app/www/BlasterJS.html", width = "100%", height = "500px")})
}