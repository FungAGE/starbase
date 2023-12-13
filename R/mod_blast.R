#' blast UI Function
#'
#' @description A shiny Module.
#'
#' @param id,input,output,session Internal parameters for {shiny}.
#'
#' @noRd
#'
#' @import XML tidyverse stringr dplyr DT Biostrings msaR chorddiag ggiraph stringi htmltools htmlwidgets 
#' 
#' @importFrom readr read_csv
#' @importFrom shiny NS tagList

mod_blast_ui <- function(id) {
  ns <- NS(id)
  dashboardBody(
    fluidRow(
      column(width=3,
      # input block
      box(
        title = "Input for BLAST/HMMER Searches",
        solidHeader = FALSE,
        collapsible = TRUE,
        width = NULL,
        tagList(
          textAreaInput(ns("query"), "Paste a sequence here in FASTA format:", value = NULL, width = "600px", height = "200px", placeholder = "Paste any nucleotide or protein sequence here"),
          fileInput(ns("query_file"), "Upload a fasta file (10 MB maximum size)",accept = c(".fa",".fna",".fasta",".faa")),

          # input logic:
          # 1) test query type
          # 2) look for captain genes
          # 3) look for cargo genes

          checkboxInput(ns("search_ship_genes"),"Screen for genes in starship sequence?",value = TRUE),
          conditionalPanel(
            condition = "input.search_ship_genes == TRUE",
            # ? is this even needed anymore?
            selectInput(ns("gene_type"), "Select a specific gene model? Default is to run search for tyr genes", 
                          selected = "tyr", 
                          multiple = TRUE, 
                          width = "200px",
                          choices = c("tyr", "fre", "nlr", "DUF3723", "plp")
                          ),
            ns = ns),
          div(
            style = "display:inline-block",
            selectInput(ns("eval"), "Filter by e-value:", choices = c(1, 0.001, 1e-4, 1e-5, 1e-10),multiple = FALSE,selected = 1)
          ),
          actionButton(ns("blast"), "Start BLAST/HMMER Search")
        )
      )),
      column(width=9,
        box(title="BLAST/HMMER Results",
            width = NULL,
            status = "success",
            closable = FALSE,
            solidHeader = FALSE,
            collapsible = FALSE,
          box(title = "Comparison of BLAST results",
            chorddiagOutput(ns("chord_ship"),width="100%",height = "600px")),
          box(title="Alignments",
            tabsetPanel(tabPanel("Starship Elements",
              DT::dataTableOutput(ns("table_ship")),
              tableOutput(ns("clicked_table_ship")),
              msaROutput(ns("alignment_ship"), width = "100%", height = "120%")
              ),
              tabPanel(title="Extracted Genes",
                tabsetPanel(
                  tabPanel("tyr",
                   DT::dataTableOutput(ns("table_tyr")),
                   tableOutput(ns("clicked_table_tyr")),
                   msaROutput(ns("alignment_tyr"))),
                  tabPanel("fre",
                   DT::dataTableOutput(ns("table_fre")),
                   tableOutput(ns("clicked_table_fre")),
                   msaROutput(ns("alignment_fre"))),
                  tabPanel("nlr",
                   DT::dataTableOutput(ns("table_nlr")),
                   tableOutput(ns("clicked_table_nlr")),
                   msaROutput(ns("alignment_nlr"))),
                  tabPanel("DUF3723",
                   DT::dataTableOutput(ns("table_DUF3723")),
                   tableOutput(ns("clicked_table_DUF3723")),
                   msaROutput(ns("alignment_DUF3723"))),
                  tabPanel("plp",
                   DT::dataTableOutput(ns("table_plp")),
                   tableOutput(ns("clicked_table_plp")),
                   msaROutput(ns("alignment_plp")))
                )))),
          box("Captain Phylogeney",girafeOutput(ns("captain_tree")))
        )
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

    # Create reactive values
    # rv <- reactiveValues()

    # separate lists of starship and gene databases with sub lists of nucl and prot databases
    # TODO: add profiles for other cargo genes? i.e. metal resistance/formaldeyde?

    db_list <- list(
      starship = list(
        nucl = "blastdb/concatenated.fa"
      ),
      gene = list(
        tyr = list(
          prot = "blastdb/YRsuperfamRefs.faa",
          nucl = "blastdb/YRsuperfamRefs.fa"
        ),
        fre = list(
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
    submitted <- eventReactive(input$blast, {

      # gather input and set up temp file
      # validate(need(input$query_file || input$query, message = TRUE))

      # throw error if none or both provided
      # BUG: file upload is always an error
      if (is.null(input$query_file) & is.null(input$query)) {
        shinyalert("Oops!", "Please provide a query sequence or fasta file", type = "error")
      }
      
      if (!is.null(input$query_file) & !is.null(input$query)) {
        shinyalert("Oops!", "Please provide either a query sequence or fasta file, not both", type = "error")
      }
      
      # choose the input
      if (is.null(input$query_file) & !is.null(input$query)) {
        query <- input$query
      } else if (!is.null(input$query_file) & is.null(input$query)) {
        query <- read_file(input$query_file, show_col_types = FALSE)
      }

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
      return(query_list)
    })

    blastresults <- reactive({
      # select correct BLAST/HMMER program
      # force the right database to be used
      # and calls the BLAST/HMMER
      # TODO: set threshold for # of returned hits?
        if (submitted()$query_type == "nucl" ) {
          blast_program <- "blastn"
        } else {
          blast_program <- "tblastn"
        }
        db <- db_list[["starship"]][["nucl"]]
        ship_blast_out<-withProgress(message = "Performing BLAST search of Starship elements in database...", {
          system(paste0(blast_program, " -query ", tmp_fasta, " -db ", db, " -evalue ", input$eval, " -outfmt 5 -max_hsps 1 -max_target_seqs 10 -num_threads 4"), intern = T) %>%
            xmlParse()
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

          read_tsv(tmp_hmmer_parsed, show_col_types = FALSE, col_names = c("query_id", "hit_IDs", "aln_length", "query_start", "query_end", "gaps", "query_seq", "subject_seq", "evalue", "bitscore"))
        }
        
        withProgress(message = "Performing HMMER search of cargo genes in database...", {
          # combine results from multiple genes
          genes_blast_out<-map(gene_list, ~ {
            run_hmmer(.x)
          })
          names(genes_blast_out) <- gene_list
        })
      }
      list("ship"=ship_blast_out,"genes"=genes_blast_out)
    })

    parsedresults_ship <- reactive({
      blast_out <- blastresults()[["ship"]]
      if (is.null(blast_out)) {} else {
          xmltop <- xmlRoot(blast_out)

          # the first chunk is for multi-fastas
          results <- xpathApply(blast_out, "//Iteration", function(row) {
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
    })
    
    parsedresults_genes <- reactive({
      if (input$search_ship_genes == TRUE) {
          blast_out <- blastresults()[["genes"]]
          if (is.null(blast_out)) {} else {
            # contains a named list of DFs
            map(blast_out,~{.x %>% 
              as.data.frame(., stringsAsFactors = FALSE) %>% 
              select(query_id, hit_IDs, aln_length, gaps, bitscore,evalue) %>% 
                slice(-1L)
              }) 
          }
      }
    })
    
    # these functions are really just for parsing the results list for genes
    make_blast_table<-function(data) {
      renderDT(
        {
          data 
        },
        selection = "single"
      )
    }

    # this chunk gets the alignment information from a clicked row
    # TODO: how to sort out which tab is being clicked?
    make_clicked_table <- function(data, name) {
      renderTable(
        {
          clicked<-input[[paste0("table_",name,"_rows_selected")]]
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

    # this function makes the alignments for clicked rows
    make_alignment<-function(type, data, name) {
        output_name<-paste0("alignment_",name)
        output[[output_name]]<-renderMsaR({
          #? needed for when row is not clicked?
          clicked<-input[[paste0("table_",name,"_rows_selected")]]
          if(!is.null(clicked)){
            if(type=="gene"){
              cdat <- data %>% slice(clicked)
              qseq <- cdat %>% pull(query_seq)
              qid <- cdat %>% pull(query_ID)
              sseq <- cdat %>% pull(subject_seq)
              sid <- cdat %>% pull(subject_IDs)
            } else if(type=="ship") {
              xmltop <- xmlRoot(data)
              
              align <- xpathApply(data, "//Iteration", function(row) {
                query <- getNodeSet(row, "Iteration_hits//Hit//Hit_hsps//Hsp//Hsp_qseq") %>% sapply(., xmlValue)
                subject <- getNodeSet(row, "Iteration_hits//Hit//Hit_hsps//Hsp//Hsp_hseq") %>% sapply(., xmlValue)
                rbind(query, subject)
              })

              # add query and subject IDs
              qid <- getNodeSet(data, "//BlastOutput_query-def") %>% sapply(., xmlValue)
              
              sids <- xpathApply(data, "//Iteration", function(row) {
                getNodeSet(row, "Iteration_hits//Hit//Hit_id") %>% sapply(., xmlValue)
              })

              alignx <- do.call("bind_cols", align)
              qseq<-alignx[1, clicked]
              sseq<-alignx[2, clicked]
            }

            if (submitted()$query_type == "nucl" ) {
              ss <- DNAStringSet(c(qseq, sseq))
              base_palette="nucleotide"
            } else {
              tmp_aln<-tempfile(fileext = ".aln")
              writeLines(str_c(str_c(str_c(">", qid,sep=""), qseq, sep = "\n"),
                              str_c(str_c(">", sids[[clicked]]), sseq, sep = "\n"),
                              sep="\n"),tmp_aln)
              ss <- seqinr::read.alignment(file = tmp_aln, format = "fasta")
              base_palette="clustal"
            }
            # BUG: when setting custom headers for msa object:
            # names(ss) <- c(paste("Query:", unique(qid[[clicked]])), paste("Subject:", sids[[clicked]]))
            # names(ss) <- c("Query:", "Subject:")        
            msaR(ss, menu = T, overviewbox = F, seqlogo = FALSE, leftheader = FALSE, labelname = TRUE, labelid = FALSE, labelNameLength = 200,colorscheme=base_palette)
          } else {}
        })
    }

    # creates tabs for each gene
    tabify<-function(name) {
      tabPanel(
        name,
        DT::DTOutput(ns(paste0("table_", name))),
        tableOutput(ns(paste0("clicked_table_", name))),
        msaROutput(ns(paste0("alignment_", name)), width = "85%", height = "120%")
      )
    }

    output$chord_ship <- renderChorddiag({
      if (is.null(parsedresults_ship())) {
      } else {
      
      # from, to, length (relative?)
      # "query_ID" = query_ID, "hit_IDs" = hit_IDs, "aln_length" = aln_length, "gaps" = gaps, "bitscore" = bitscore, "eval" = eval
      if(nrow(parsedresults_ship())<10){
        nhits<-nrow(parsedresults_ship())
      } else {
        nhits<-10
      }
        rexp <- "^(\\w+)\\s?(.*)$"
        df<-parsedresults_ship() %>%
          mutate(across(c(query_ID,hit_IDs),~sub(rexp,"\\1",.)),
                 aln_length=as.numeric(aln_length)) %>%
          mutate(aln_length=ifelse(query_ID == hit_IDs,NA,aln_length)) %>% # exclude self-hits
          arrange(desc(eval)) %>%
          slice_head(n=nhits) %>%
          select(query_ID,hit_IDs,aln_length)
        
      subNames<-unique(df$hit_IDs)
      queryName<-"QUERY"
      # queryName<-unique(df$query_ID)

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
    }})
    
    output$table_ship <- DT::renderDT(
      {
        if (is.null(parsedresults_ship())) {
        } else {
          parsedresults_ship() %>% as.data.frame(., stringsAsFactors = FALSE)
        } 
      },
      selection = "single"
    )
    
    output$clicked_table_ship <- renderTable(
      {
        blast_out <- blastresults()[["ship"]]
        
          if (is.null(input$table_ship_rows_selected)) {
          } else {
          xmltop <- xmlRoot(blast_out)
          clicked <- input$table_ship_rows_selected
          tableout <- parsedresults_ship() %>%
            data.frame() %>%
            slice(clicked)
  
          tableout <- t(tableout)
          names(tableout) <- c("")
          rownames(tableout) <- c("Query ID", "Hit ID", "Alignment Length", "Gaps", "Bit Score", "e-value")
          colnames(tableout) <- NULL
          data.frame(tableout)
          }
          },
      rownames = T,
      colnames = F
    )
    
    make_alignment("ship",blastresults()[["ship"]],"ship")
    
    output$table_tyr<-make_blast_table(parsedresults_genes()[["tyr"]])
    output$clicked_table_tyr<-make_clicked_table(parsedresults_genes()[["tyr"]],"tyr")
    make_alignment("gene",parsedresults_genes()[["tyr"]],"tyr")

    output$table_fre<-make_blast_table(parsedresults_genes()[["fre"]])
    output$clicked_table_fre<-make_clicked_table(parsedresults_genes()[["fre"]],"fre")
    output$alignment_fre<-make_alignment("gene",parsedresults_genes()[["fre"]],"fre")

    output$table_nlr<-make_blast_table(parsedresults_genes()[["nlr"]])
    output$clicked_table_nlr<-make_clicked_table(parsedresults_genes()[["nlr"]],"nlr")
    output$alignment_nlr<-make_alignment("gene",parsedresults_genes()[["nlr"]],"nlr")

    output$table_DUF3723<-make_blast_table(parsedresults_genes()[["DUF3723"]])
    output$clicked_table_DUF3723<-make_clicked_table(parsedresults_genes()[["DUF3723"]],"DUF3723")
    output$alignment_DUF3723<-make_alignment("gene",parsedresults_genes()[["DUF3723"]],"DUF3723")

    output$table_plp<-make_blast_table(parsedresults_genes()[["plp"]])
    output$clicked_table_plp<-make_clicked_table(parsedresults_genes()[["plp"]],"plp")
    output$alignment_plp<-make_alignment("gene",parsedresults_genes()[["plp"]],"plp")

  output$gene_tabs<- renderUI({
    gene_tabs <- map(names(parsedresults_genes()), ~ {
      tabify(.x)
    })
  do.call(tabsetPanel, gene_tabs)
  })

  output$captain_tree<-eventReactive(input$blast,
  renderGirafe({
    # if(input$search_ship_genes == TRUE & !is.null(parsedresults_genes()[["tyr"]]) & submitted()$query_type == "prot"){
    readRDS(file = "data/captain-tree.RDS")
      # } else {
        # shinyalert("Error!", "Can't run shoot with this combination (yet)", type = "error")
      # }
    })
  )
})
}