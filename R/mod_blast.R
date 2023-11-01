#' blast UI Function
#'
#' @description A shiny Module.
#'
#' @param id,input,output,session Internal parameters for {shiny}.
#'
#' @noRd
#'
#' @import XML tidyverse stringr dplyr DT Biostrings msaR
#' @importFrom readr read_csv
#' @importFrom shiny NS tagList

mod_blast_ui <- function(id) {
  ns <- NS(id)
  dashboardBody(
    fluidRow(
      # input block
      box(
        title = "BLAST/HMMER Searches",
        solidHeader = FALSE,
        collapsible = TRUE,
        width = NULL,
        # TODO: add logic to only accept one query type
        tagList(
          textAreaInput(ns("query"), "Paste your sequence here:", value = "", width = "600px", height = "200px", placeholder = "Paste any string of nucleotide or protein sequence here"),
          fileInput(ns("query_file"), "Upload a fasta file", width = "100px"),
          
          # choice should be limited to 1) whole starship or 2) just a gene sequence
          # if 2) then choose what gene or unknown
          selectInput(ns("input_type"), "What is the input sequence?", choices = c("starship", "gene"), width = "120px"),
          # BUG: where db_type when choice = gene is not overwritten from default db_type when choice = starship
          conditionalPanel(
            condition = "input.input_type == 'gene'",
            selectInput(ns("gene_type"), "Select a gene model:", selected = character(0), multiple = TRUE, choices = c("tyr", "freB", "nlr", "DUF3723", "plp"), width = "120px"),
            ns = ns
          ),
          conditionalPanel(
            condition = "input.input_type == 'starship'",
            selectInput(ns("search_ship_genes"), "Search for genes in starship sequence?", choices = c("No", "Yes"), width = "120px"),
            ns = ns
          ),
          div(
            style = "display:inline-block",
            selectInput(ns("eval"), "e-value:", choices = c(1, 0.001, 1e-4, 1e-5, 1e-10), width = "120px")
          ),
          actionButton(ns("blast"), "Search")
        )
      )
    ),
    
    # TODO: should be separated by results from different genes?
    # Basic results output
    fluidRow(
      box(
        title = "Results",
        width = NULL,
        status = "success",
        closable = FALSE,
        solidHeader = FALSE,
        collapsible = FALSE,
        uiOutput(ns("tabs"))
      )
    )
  )
}

#' blast Server Functions
#'
#' @noRd

mod_blast_server <- function(id) {
  moduleServer(id, function(input, output, session) {
    ns <- session$ns
    blastresults <- eventReactive(input$blast, {
      # separate lists of starship and gene databases with sub lists of nucl and prot databases
      
      # TODO: add profiles for other cargo genes? i.e. metal resistance/formaldeyde?
      
      db_list <- list(
        starship = list(
          nucl = "blastdb/concatenated.fa",
          prot = ""
        ),
        gene = list(
          tyr = list(
            prot = "blastdb/YRsuperfamRefs.faa",
            nucl = "blastdb/YRsuperfamRefs.fa"
          ),
          freB = list(
            prot = "blastdb/fre.mycoDB.faa",
            nucl = "blastdb/fre.fa"
          ),
          nlr = list(
            prot = "blastdb/nlr.mycoDB.faa",
            nucl = "blastdb/nlr.fa"
          ),
          DUF3723 = list(
            prot = "blastdb/duf3723.mycoDB.faa",
            nucl = "blastdb/duf3723.fa"
          ),
          plp = list(
            prot = "blastdb/plp.mycoDB.faa",
            nucl = "blastdb/plp.fa"
          )
        )
      )
      
      # gather input and set up temp file
      tmp_fasta <- tempfile(fileext = ".fa")
      
      # TODO: add better error catching for format of query sequence/file here
      
      # gather input and set up temp file
      # validate(need(input$query_file || input$query, message = TRUE))
      
      if (is.null(input$query_file)) {
        query <- input$query
      } else {
        query <- read_file(input$query_file, show_col_types = FALSE)
      }
      
      # function to clean the query sequence: first remove non-letter characters, then ambiguous characters, then guess if query is nucl or protein
      clean_lines <- function(query) {
        cleaned_seq <- NULL
        na_count <- 0
        na_char <- c("N", "X", "n", "x", "*")
        nucl_char <- c("a", "A", "t", "T", "g", "G", "c", "C")
        prot_char <- c(letters, LETTERS)
        prot_char <- prot_char[!prot_char %in% c(na_char, nucl_char)]
        
        query_lines <- str_split(query, "\n", simplify = T)
        
        for (line in query_lines) {
          if (str_length(line) > 1) {
            if (grepl("^>", trimws(line))) {
              header <- gsub("^>", "", trimws(line))
              # header<-str_split(line,"\n",simplify=TRUE)[1]
              # cleaned_seq <- str_c(cleaned_seq, header, "\n")
            } else {
              # Remove non-letter characters using gsub
              cleaned_line <- gsub(paste0(na_char, collapse = "|"), "", fixed = T, gsub("[^A-Za-z]", "", trimws(line)))
              na_count <- na_count + sum(table(unlist(strsplit(cleaned_line, "")))[na_char], na.rm = TRUE)
              cleaned_seq <- str_c(cleaned_seq, cleaned_line)
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
      
      # select correct BLAST/HMMER program
      # force the right database to be used
      # and calls the BLAST/HMMER
      
      if (input$input_type == "gene" | input$search_ship_genes == "Yes") {
        if (query_list$query_type == "nucl") {
          hmmer_program <- "jackhmmer"
        } else {
          hmmer_program <- "phmmer"
        }
        if (!is.null(input$gene_type)) {
          gene_list <- input$gene_type
        } else {
          gene_list <- names(db_list[["gene"]])
        }
        
        run_hmmer <- function(gene) {
          db <- db_list[[input$input_type]][[gene]][[query_list$query_type]]
          tmp_hmmer <- tempfile(fileext = ".hmmer")
          tmp_hmmer_parsed <- tempfile(fileext = ".hmmer.parsed")
          hmmer_cmd <- paste(hmmer_program, "-o", tmp_hmmer, "--cpu 4", "--domE", input$eval, tmp_fasta, db, sep = " ")
          system(hmmer_cmd, intern = FALSE)
          system(paste0("python bin/hmm.py -i ", tmp_hmmer, "-o ", tmp_hmmer_parsed))
          # {record.id}\t{hit.id}\t{aln_length}\t{query_seq}\t{query_start}\t{query_end}\t{subject_seq}\t{evalue}\t{bitscore}
          read_csv(tmp_hmmer_parsed, show_col_types = FALSE, col_names = c("query_id", "hit_IDs", "aln_length", "query_start", "query_end", "gaps", "query_seq", "subject_seq", "evalue", "bitscore"))
        }
        
        withProgress(message = "Performing HMMER search of cargo genes in database...", {
          # combine results from multiple genes
          map(gene_list, ~ {
            run_hmmer(.x)
          })
        })
      }
      
      if (input$input_type == "starship") {
        # catch for protein queries in ship blastdb
        if (query_list$query_type == "prot") {
          shinyalert("Error:", "Cannot search starship database with protein sequence", type = "error")
        }
        
        if (query_list$query_type == "nucl") {
          blast_program <- "blastn"
        } else {
          blast_program <- "blastp"
        }
        db <- db_list[[input$input_type]][[query_list$query_type]]
        withProgress(message = "Performing BLAST search of Starship elements in database...", {
          system(paste0(blast_program, " -query ", tmp_fasta, " -db ", db, " -evalue ", input$eval, " -outfmt 5 -max_hsps 1 -max_target_seqs 10 -num_threads 4"), intern = T) %>%
            xmlParse()
        })
      }
    })
    
    parsedresults <- reactive({
      if (is.null(blastresults())) {} else {
        if (input$input_type == "starship") {
          xmltop <- xmlRoot(blastresults())
          
          # the first chunk is for multi-fastas
          results <- xpathApply(blastresults(), "//Iteration", function(row) {
            query_ID <- getNodeSet(row, "Iteration_query-def") %>% sapply(., xmlValue)
            hit_IDs <- getNodeSet(row, "Iteration_hits//Hit//Hit_id") %>% sapply(., xmlValue)
            aln_length <- getNodeSet(row, "Iteration_hits//Hit//Hsp_align-len") %>% sapply(., xmlValue)
            gaps <- getNodeSet(row, "Iteration_hits//Hit//Hsp_gaps") %>% sapply(., xmlValue)
            bitscore <- getNodeSet(row, "Iteration_hits//Hit//Hit_hsps//Hsp//Hsp_bit-score") %>% sapply(., xmlValue)
            eval <- getNodeSet(row, "Iteration_hits//Hit//Hit_hsps//Hsp//Hsp_evalue") %>% sapply(., xmlValue)
            bind_cols("query_ID" = query_ID, "hit_IDs" = hit_IDs, "aln_length" = aln_length, "gaps" = gaps, "bitscore" = bitscore, "eval" = eval)
          })
          
          # this ensures that NAs get added for no hits
          lapply(results, function(y) {
            as.data.frame((y), stringsAsFactors = FALSE)
          }) %>% bind_rows()
        }
        
        if (input$input_type == "gene" | input$search_ship_genes == "Yes") {
          # contains a named list of DFs
          map(blastresults(), ~ {
            .x %>% select(query_id, hit_IDs, aln_length, gaps, evalue, bitscore)
          })
        }
      }
    })
    
    # these functions are really just for parsing the results list for genes
    make_blast_table <- function(data) {
      renderDT(
        {
          data
        },
        selection = "single"
      )
    }
    
    # this chunk gets the alignment information from a clicked row
    # TODO: how to sort out which tab is being clicked?
    make_clicked_table <- function(data, clicked) {
      renderTable(
        {
          tableout <- data %>%
            slice(clicked) %>%
            t()
          names(tableout) <- c("")
          rownames(tableout) <- c("Query ID", "Hit ID", "Alignment Length", "Gaps", "Bit Score", "e-value")
          colnames(tableout) <- NULL
          data.frame(tableout)
        },
        rownames = T,
        colnames = F
      )
    }
    
    make_alignment <- function(data, clicked) {
      renderMsaR({
        cdat <- data %>% slice(clicked)
        qseq <- cdat %>% pull(query_seq)
        qid <- cdat %>% pull(query_ID)
        sseq <- cdat %>% pull(subject_seq)
        sid <- cdat %>% pull(subject_IDs)
        
        ss <- DNAStringSet(c(qseq, sseq))
        names(ss) <- c(paste("Query:", qid), paste("Subject:", sid))
        msaR(ss, menu = T, overviewbox = F, seqlogo = FALSE, leftheader = FALSE, labelname = TRUE, labelid = FALSE, labelNameLength = 200)
      })
    }
    
    # creates tabs for each gene
    tabify <- function(dataset) {
      tabPanel(
        dataset,
        DT::DTOutput(paste0("table_", dataset))
        # tableOutput(ns(paste0("clicked_table_", dataset))),
        # msaROutput(ns(paste0("alignment_", dataset)), width = "100%", height = "100%")
      )
    }
    
    output$tabs <- renderUI({
      if (input$input_type == "starship" & input$search_ship_genes != "Yes") {
        selected_datasets <- input$input_type
      } else if (input$input_type == "gene") {
        selected_datasets <- input$gene_type
      } else if (input$input_type == "starship" & input$search_ship_genes == "Yes") {
        selected_datasets <- c(input$input_type, input$gene_type)
      }
      # tab_list <- map(selected_datasets, ~ {
      #   tabify(.x)
      tabs <- lapply(selected_datasets, function(dataset) {
        tabPanel(dataset, DT::dataTableOutput(ns(paste0("table_", dataset))))
      })
      do.call(tabsetPanel, tabs)
    })
    
    
    observe({
      if (is.null(blastresults())) {
      } else {
        if (input$input_type == "starship") {
          output[[paste0("table_", input$input_type)]] <- make_blast_table(parsedresults())
          clicked_ship <- input[[paste0("table_", input$input_type, "_rows_selected")]]
          
          if (is.null(clicked_ship)) {
          } else {
            output[[paste0("clicked_table_", input$input_type)]] <- renderTable(
              {
                xmltop <- xmlRoot(blastresults())
                tableout <- parsedresults() %>%
                  data.frame() %>%
                  slice(clicked_ship)
                
                tableout <- t(tableout)
                names(tableout) <- c("")
                rownames(tableout) <- c("Query ID", "Hit ID", "Alignment Length", "Gaps", "Bit Score", "e-value")
                colnames(tableout) <- NULL
                data.frame(tableout)
              },
              rownames = T,
              colnames = F
            )
            
            # this chunk makes the alignments for clicked rows
            output[[paste0("alignment_", input$input_type)]] <- renderMsaR({
              xmltop <- xmlRoot(blastresults())
              # loop over the xml to get the alignments
              align <- xpathApply(blastresults(), "//Iteration", function(row) {
                query <- getNodeSet(row, "Iteration_hits//Hit//Hit_hsps//Hsp//Hsp_qseq") %>% sapply(., xmlValue)
                subject <- getNodeSet(row, "Iteration_hits//Hit//Hit_hsps//Hsp//Hsp_hseq") %>% sapply(., xmlValue)
                rbind(query, subject)
              })
              
              # add query and subject IDs
              qid <- getNodeSet(blastresults(), "//BlastOutput_query-def") %>% sapply(., xmlValue)
              sids <- xpathApply(blastresults(), "//Iteration", function(row) {
                sid <- getNodeSet(row, "Iteration_hits//Hit//Hit_id") %>% sapply(., xmlValue)
              })
              
              alignx <- do.call("cbind", align)
              ss <- DNAStringSet(c(alignx[1, clicked], alignx[2, clicked]))
              names(ss) <- c(paste("Query:", qid), paste("Subject:", sids[[clicked]]))
              
              msaR(ss, menu = T, overviewbox = F, seqlogo = FALSE, leftheader = FALSE, labelname = TRUE, labelid = FALSE, labelNameLength = 200)
            })
          }
        }
        if (input$input_type == "gene" | input$search_ship_genes == "Yes") {
          for (dataset in names(parsedresults())) {
            output[[paste0("table_", dataset)]] <- make_blast_table(parsedresults()[[dataset]])
            clicked_gene <- input[[paste0("table_", dataset, "_rows_selected")]]
            output[[paste0("clicked_table_", dataset)]] <- make_clicked_table(parsedresults()[[dataset]], clicked_gene)
            output[[paste0("alignment_", dataset)]] <- make_alignment(parsedresults()[[dataset]], clicked_gene)
          }
        }
      }
    })
  })
}
