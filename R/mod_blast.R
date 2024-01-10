#' blast UI Function
#'
#' @description A shiny Module.
#'
#' @param id,input,output,session Internal parameters for {shiny}.
#'
#' @noRd
#'
#' @import shinydashboard shinydashboardPlus
#' @import XML tidyverse stringr dplyr DT Biostrings msaR chorddiag ggiraph stringi htmltools htmlwidgets 
#' 
#' @importFrom readr read_csv read_file
#' @importFrom shiny NS tagList

# input logic:
# 1) test query type
# 2) look for captain genes
# 3) look for cargo genes


mod_blast_ui <- function(id) {
  ns <- NS(id)
  dashboardBody(
    fluidRow(
      column(width=4,
      box(
        title = "Input for BLAST/HMMER Searches",
        id = "blastbox",
        solidHeader = FALSE,
        collapsible = TRUE,
        width = NULL,
        tagList(
          textAreaInput(ns("query_text"), "Paste a sequence here in FASTA format:", value = NULL, width = "600px", height = "200px", placeholder = "Paste any nucleotide or protein sequence here"),
          fileInput(ns("query_file"), "Upload a fasta file (10 MB maximum size)",accept = c(".fa",".fna",".fasta",".faa")),

          checkboxInput(ns("search_ship_genes"),"Screen for genes in starship sequence?",value = TRUE),
          conditionalPanel(
            condition = "input.search_ship_genes == TRUE",
            selectInput(
              ns("gene_type"), 
              "Select a specific gene model? Default is to run search for tyr genes", 
              selected = "tyr", 
              multiple = TRUE, 
              width = "200px",
              choices = c("tyr", "fre", "nlr", "DUF3723", "plp")
            ),
            ns = ns
          ),
          div(
            style = "display:inline-block",
            selectInput(ns("eval"), "Filter by e-value:", choices = c(1, 0.001, 1e-4, 1e-5, 1e-10),multiple = FALSE,selected = 1)
          ),
          rep_br(1),
          actionButton(ns("blast"), "Start BLAST/HMMER Search")
        )
      )),
      column(width = 8,
        uiOutput(ns("ship_ui")),
        uiOutput(ns("gene_ui"))
      )
      # box(title="Classification",
      #   solidHeader = FALSE,
      #   collapsible = TRUE,
      #   width = NULL,
      #   uiOutput(ns("classification_ui")))
    )
  )
}

#' blast Server Functions
#'
#' @noRd

mod_blast_server <- function(id) {
  moduleServer(id, function(input, output, session) {
    ns <- session$ns

    # Create reactive values
    # rv <- reactiveValues()

    # separate lists of starship and gene databases with sub lists of nucl and prot databases
    # TODO: add profiles for other cargo genes? i.e. metal resistance/formaldeyde?

    db_list <- list(
      starship = list(
        nucl = "Starships/blastdb/concatenated.fa"
      ),
      gene = list(
        tyr = list(
          prot = "Starships/blastdb/YRsuperfamRefs.faa",
          nucl = "Starships/blastdb/YRsuperfamRefs.fa"
        ),
        fre = list(
          prot = "Starships/blastdb/fre.mycoDB.faa",
          nucl = "Starships/blastdb/fre.fa"
        ),
        nlr = list(
          prot = "Starships/blastdb/nlr.mycoDB.faa",
          nucl = "Starships/blastdb/nlr.fa"
        ),
        DUF3723 = list(
          prot = "Starships/blastdb/duf3723.mycoDB.faa",
          nucl = "Starships/blastdb/duf3723.fa"
        ),
        plp = list(
          prot = "Starships/blastdb/plp.mycoDB.faa",
          nucl = "Starships/blastdb/plp.fa"
        )
      )
    )
    
    # gather input and set up temp file
    tmp_fasta <- tempfile(fileext = ".fa")

    observeEvent(input$blast, {
      updateBox("blastbox", action = "toggle")
    })

    submitted <- eventReactive(input$blast, {
      # throw error if none or both provided
      # choose the input
      if (all(is.null(input$query_file$datapath), input$query_text == "")) {
        shinyalert("No input provided:", "Please provide either a query sequence in the text box, or upload a fasta file", type = "error")
      } else if (all(!is.null(input$query_file$datapath),input$query_text != "")) {
        shinyalert("Multiple inputs:", "Please provide either a query sequence in the text box, or upload a fasta file, not both", type = "error")
      } 
      
      if (all(is.null(input$query_file$datapath),input$query_text != "")) {
        query <- input$query_text
      } else if (all(!is.null(input$query_file$datapath), input$query_text == "")) {
        fasta_file<-seqinr::read.fasta(input$query_file$datapath)
        fasta_file_header<-ifelse(any(is.null(names(fasta_file)),names(fasta_file) == "",is.na(names(fasta_file))),"QUERY",names(fasta_file))
        query <-str_c(paste0(">",fasta_file_header),str_flatten(fasta_file),sep="\n")
      }

      if (is.null(query)){
        shinyalert("Input error","something wrong with fasta input")
      }

      # if(!(length(grep(">", query)) > 1)){ shinyalert("Multifasta provided:","Multifastas are currently not supported. Upload one sequence at a time please!")}

      # if gene scan is chosen, use this to create separate results section
      input_type <- ifelse(input$search_ship_genes == TRUE, "gene","ship")
      
      # TODO: add better error catching for format of query sequence/file here
      # function to clean the query sequence: first remove non-letter characters, then ambiguous characters, then guess if query is nucl or protein
      clean_lines <- function(query) {
        cleaned_seq <- NULL
        na_count <- 0
        na_char <- c("N", "X", "n", "x", "*")
        nucl_char <- c("a", "A", "t", "T", "g", "G", "c", "C")
        prot_char <- c(letters, LETTERS)
        prot_char <- prot_char[!prot_char %in% c(na_char, nucl_char)]

        query_lines <- str_split(query, "\n", simplify = T)
        # BUG: header is not being added correctly in some cases
        for (line in query_lines) {
          if (str_length(line) > 1) {
            if (grepl("^>", trimws(line))) {
              header <- gsub("^>", "", trimws(line))
              # header<-str_split(line,"\n",simplify=TRUE)[1]
              cleaned_seq <- header
            } else {
              # Remove non-letter characters using gsub
              cleaned_line <- gsub(paste0(na_char, collapse = "|"), "", fixed = T, gsub("[^A-Za-z]", "", trimws(line)))
              na_count <- na_count + sum(table(unlist(strsplit(cleaned_line, "")))[na_char], na.rm = TRUE)
              cleaned_seq <- str_c(cleaned_seq, cleaned_line,sep="\n")
            }
          }
        }

        # count characters for deciding type later
        nucl_count <- sum(table(unlist(strsplit(cleaned_seq, "")))[nucl_char], na.rm = TRUE)
        prot_count <- sum(table(unlist(strsplit(cleaned_seq, "")))[prot_char], na.rm = TRUE)

        # guess if sequence is nucl or protein
        if (prot_count >= (0.1 * str_length(cleaned_seq))) {
          query_type <- "prot"
          print("Query is protein sequence")
        } else {
          query_type <- "nucl"
          print("Query is nucleotide sequence")
        }
        return(list("header" = header, "query_type" = query_type, "query" = cleaned_seq))
      }

      query_list <- clean_lines(query)

      writeLines(str_c(str_c(">", query_list$header, sep = ""), query_list$query, sep = "\n"), tmp_fasta)
      return(query_list)
    })

    blastresults <- reactive({
      # select correct BLAST/HMMER program
      # force the right database to be used
      # and calls the BLAST/HMMER
      # TODO: set threshold for # of returned hits?
      ship_blast_results<-NULL
      ship_blast_out<-NULL
      withProgress(message = "Performing BLAST search of Starship elements in database...", {
        if (submitted()$query_type == "nucl" ) {
          blast_program <- "blastn"
        } else {
          blast_program <- "tblastn"
        }
        db <- db_list[["starship"]][["nucl"]]
        ship_blast_out<-system(paste0(blast_program, " -query ", tmp_fasta, " -db ", db, " -evalue ", input$eval, " -outfmt 5 -max_hsps 1 -max_target_seqs 10 -num_threads 4"), intern = T) %>%
            xmlParse()
        xmltop <- xmlRoot(ship_blast_out)

        # the first chunk is for multi-fastas
        ship_blast_results <- xpathApply(ship_blast_out, "//Iteration", function(row) {
          query_id <- getNodeSet(row, "Iteration_query-def") %>% sapply(., xmlValue)
          hit_IDs <- getNodeSet(row, "Iteration_hits//Hit//Hit_id") %>% sapply(., xmlValue)
          aln_length <- getNodeSet(row, "Iteration_hits//Hit//Hsp_align-len") %>% sapply(., xmlValue)
          gaps <- getNodeSet(row, "Iteration_hits//Hit//Hsp_gaps") %>% sapply(., xmlValue)
          bitscore <- getNodeSet(row, "Iteration_hits//Hit//Hit_hsps//Hsp//Hsp_bit-score") %>% sapply(., xmlValue)
          eval <- getNodeSet(row, "Iteration_hits//Hit//Hit_hsps//Hsp//Hsp_evalue") %>% sapply(., xmlValue)
          query_seq <- getNodeSet(row, "Iteration_hits//Hit//Hit_hsps//Hsp//Hsp_qseq") %>% sapply(., xmlValue)
          subject_seq <- getNodeSet(row, "Iteration_hits//Hit//Hit_hsps//Hsp//Hsp_hseq") %>% sapply(., xmlValue)
          bind_cols("query_id" = query_id, "hit_IDs" = hit_IDs, "aln_length" = aln_length, "gaps" = gaps, "query_seq" = query_seq, "subject_seq" = subject_seq, "bitscore" = bitscore, "eval" = eval)
        }) %>%
        lapply(., function(y) {
          as.data.frame((y), stringsAsFactors = FALSE) # this ensures that NAs get added for no hits
        }) %>% bind_rows()
      })
        
      genes_blast_out <- NULL
      if (input$search_ship_genes == TRUE) {
        if (submitted()$query_type == "nucl") {
          hmmer_program <- "jackhmmer"
        } else {
          hmmer_program <- "phmmer"
        }

        gene_list <- input$gene_type
        
        run_hmmer <- function(gene) {
          db <- db_list[["gene"]][[gene]][[submitted()$query_type]]
          
          tmp_hmmer <- tempfile(fileext = ".hmmer")
          
          hmmer_cmd <- paste(hmmer_program, "-o", tmp_hmmer, "--cpu 4", "--domE", input$eval, tmp_fasta, db, sep = " ")
          system(hmmer_cmd, intern = FALSE)
          
          tmp_hmmer_parsed <- tempfile(fileext = ".hmmer.parsed")

          system(paste0("python bin/hmm.py --hmmer_output_file ", tmp_hmmer, " --parsed ", tmp_hmmer_parsed))

          read_tsv(tmp_hmmer_parsed, show_col_types = FALSE) %>%
            dplyr::select(-c(query_start,query_end))
            # dplyr::slice(-1L)
            # dplyr::slice_min(order_by=evalue,with_ties=FALSE)

        }
        
        withProgress(message = "Performing HMMER search of cargo genes in database...", {
          # combine results from multiple genes
          genes_blast_out<-map(gene_list, ~ {
            run_hmmer(.x)
          })
          names(genes_blast_out) <- gene_list
        })
      }
      list("ship"=ship_blast_results,"genes"=genes_blast_out)
    })
        
    # these functions are really just for parsing the results list for genes
    make_blast_table<-function(data) {
      data %>% 
        dplyr::select(-c(query_id,query_seq,subject_seq)) %>% 
        datatable(data=.,selection = "single")
    }

    # this chunk gets the alignment information from a clicked row
    # TODO: this info is redundant with main blast table. include other info like taxonomy? links to pages?
    make_clicked_table <- function(data, name) {
      tableout<-data %>% dplyr::select(-c(query_id,query_seq,subject_seq)) %>%
        dplyr::slice(input[[paste0("table_",name,"_rows_selected")]]) #%>% t()

      tableout %>% pull(hit_IDs)
        # names(tableout) <- c("")
        # rownames(tableout) <- c("Query ID", "Hit ID", "Alignment Length", "Gaps", "Bit Score", "e-value")
        # colnames(tableout) <- NULL
        # data.frame(tableout)
    }

    # this function makes the alignments for clicked rows
    make_alignment<-function(data, name) {
      #? needed for when row is not clicked?
      clicked<-input[[paste0("table_",name,"_rows_selected")]]
      if (!is.null(clicked)) {
        cdat <- data %>% dplyr::slice(clicked)
        qseq <- cdat %>% pull(query_seq)
        qid <- cdat %>% pull(query_id)
        sseqs <- cdat %>% pull(subject_seq)
        sids <- cdat %>% pull(hit_IDs)

        if (submitted()$query_type == "nucl" ) {
          ss <- DNAStringSet(c(qseq, sseqs[[clicked]]))
          base_palette="nucleotide"
        } else {
          tmp_aln<-tempfile(fileext = ".aln")
          writeLines(str_c(str_c(str_c(">", qid,sep=""), qseq, sep = "\n"),
                          str_c(str_c(">", sids[[clicked]]), sseqs[[clicked]], sep = "\n"),
                          sep="\n"),tmp_aln)
          ss <- seqinr::read.alignment(file = tmp_aln, format = "fasta")
          base_palette="clustal"
        }
        # BUG: when setting custom headers for msa object:
        # names(ss) <- c("Query:", "Subject:")        
        msaR(ss, menu = T, overviewbox = F, seqlogo = FALSE, leftheader = FALSE, labelname = TRUE, labelid = FALSE, labelNameLength = 200,colorscheme=base_palette)
      } else {}
    }

    make_chord<-function(data){     
        if (is.null(data)) {
          print("input error for chord")
        } else {
        
        # from, to, length (relative?)
        # "query_id" = query_id, "hit_IDs" = hit_IDs, "aln_length" = aln_length, "gaps" = gaps, "bitscore" = bitscore, "eval" = eval
        if(nrow(data)<10){
          nhits<-nrow(data)
        } else {
          nhits<-10
        }
          rexp <- "^(\\w+)\\s?(.*)$"
          df<-data %>%
            mutate(across(c(query_id,hit_IDs),~sub(rexp,"\\1",.)),
                  aln_length=ifelse(query_id == hit_IDs,NA,as.numeric(aln_length))) %>% # exclude self-hits
            arrange(desc(eval)) %>%
            dplyr::slice_head(n=nhits) %>%
            select(query_id,hit_IDs,aln_length)
          
        subNames<-unique(df$hit_IDs)
        queryName<-"QUERY"
        # queryName<-unique(df$query_id)

        groupColors<-RColorBrewer::brewer.pal(nhits,"Paired")
        
        m <- matrix(df$aln_length,
                    byrow = FALSE,
                    nrow = nhits, ncol = nhits)

        dimnames(m) <- list(query = rep(queryName,nhits),
                            subject = subNames)

        # Build the chord diagram:
        chorddiag(m, 
                  groupnameFontsize = 12,
                  groupColors = groupColors, 
                  groupnamePadding = 50)
    }}
    
    render_output <- function(name, plot.chord = FALSE) {
      req(blastresults())  # Check if blastresults() is available

      output[[paste0("table_",name)]] <- DT::renderDT({ make_blast_table(blastresults()[[name]]) })
      output[[paste0("clicked_table_",name)]] <- renderTable({ make_clicked_table(blastresults()[[name]], name) }, rownames = TRUE, colnames = FALSE)
      output[[paste0("alignment_",name)]] <- renderMsaR({ make_alignment(blastresults()[[name]], name) })

      tab_content <- list(
        DT::DTOutput(ns(paste0("table_", name))),
        tableOutput(ns(paste0("clicked_table_", name))),
        msaROutput(ns(paste0("alignment_", name)), width = "85%", height = "120%")
      )
      if (plot.chord) {
        output[[paste0("chord_",name)]] <- renderChorddiag({ make_chord(blastresults()[[name]]) })
        tab_content <- c(chorddiagOutput(ns(paste0("chord_", name)), width = "100%", height = "600px"), tab_content)
      }
      return(tab_content)
    }

    output$ship_ui <- renderUI({
      box(title = "Ship BLAST/HMMER Results",
          solidHeader = FALSE,
          collapsible = TRUE,
          width = NULL,
          render_output("ship", TRUE)
      )
    })

    output$gene_ui<-renderUI({
        box(title="Gene BLAST/HMMER Results",
        solidHeader = FALSE,
        collapsible = TRUE,
        width = NULL,
        # map(names(blastresults()[["genes"]]), ~ {
        #   render_output(.x, FALSE)})
        render_output("tyr", FALSE))
    })

    output$captain_tree<-renderGirafe(readRDS("data/captain-tree.RDS"))
    output$classification_ui <- renderUI({girafeOutput(ns("captain_tree"))})
  })
}