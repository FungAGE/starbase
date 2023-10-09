#' blast UI Function
#'
#' @description A shiny Module.
#'
#' @param id,input,output,session Internal parameters for {shiny}.
#'
#' @noRd
#'
#' @import XML
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
      width=NULL,
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
          selectInput(ns("db_type"), "Which database do you want to search?", selected = character(0), multiple = FALSE, choices = c("curated", "starfish"), width = "120px"),
          selectInput(ns("search_ship_genes"), "Search for genes in starship sequence?", choices = c("No", "Yes"), width = "120px"),
          ns = ns
        ),
        div(
          style = "display:inline-block",
          selectInput(ns("eval"), "e-value:", choices = c(1, 0.001, 1e-4, 1e-5, 1e-10), width = "120px")
        ),
        actionButton(ns("blast"), "Search")
      )
    )),

    # Basic results output
    fluidRow(
      box(
      title = "Results",
      width = NULL,
      status = "success",
      closable = FALSE,
      solidHeader = FALSE,
      collapsible = TRUE,
      DT::DTOutput(ns("tbl")),
      p("Alignment:", tableOutput(ns("clicked"))),
      verbatimTextOutput(ns("alignment"))
    )
  ))
}

#' blast Server Functions
#'
#' @noRd

mod_blast_server <- function(id) {
  moduleServer(id, function(input, output, session) {
    # ns <- session$ns
    blastresults <- eventReactive(input$blast, {
      # separate lists of starship and gene databases with sub lists of nucl and prot databases
      # TODO: add starship protein databases and gene hmm to folder
      db_list <- list(
        starship = list(
          curated =
            list(
              nucl = "/home/adrian/Systematics/Starship_Database/blastdb/concatenated.fa",
              prot = ""
            ),
          starfish =
            list(
              nucl = "/home/adrian/Systematics/Starship_Database/blastdb/mycodb.final.starships.headers.fna",
              prot = ""
            )
        ),
        gene = list(
          tyr = list(
            prot = "/home/adrian/Systematics/Starship_Database/blastdb/YRsuperfamRefs.faa",
            nucl = ""
          ),
          freB = list(
            prot = "",
            nucl = ""
          ),
          nlr = list(
            prot = "",
            nucl = ""
          ),
          DUF3723 = list(
            prot = "",
            nucl = ""
          ),
          plp = list(
            prot = "",
            nucl = ""
          )
        )
      )

      # gather input and set up temp file
      tmp_fasta <- "tmp/queryfile"

      # TODO: add better error catching for format of query sequence/file here

      # gather input and set up temp file
      # validate(need(input$query_file || input$query, message = TRUE))

      if (is_null(input$query_file)) {
        query <- input$query
      } else {
        query <- read_file(input$query_file,show_col_types = FALSE)
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

        # guess nucl or protein
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

      # force the right database to be used
      # and calls the BLAST/HMMER
      if (input$input_type == "starship") {
        db <- db_list[[input$input_type]][[input$db_type]][[query_list$query_type]]
      } else {
        db <- db_list[[input$input_type]][[input$gene_type]][[query_list$query_type]]
      }
      # select correct BLAST/HMMER program
      if (input$input_type == "gene" | input$search_ship_genes == "Yes") {
        if (query_list$query_type == "nucl") {
          hmmer_program <- "jackhmmer"
        } else {
          hmmer_program <- "phmmer"
        }
        if (input$input_type == "starship" & input$search_ship_genes == "Yes") {
          gene_list <- names(db_list[["gene"]])
        } else {
          gene_list <- input$gene_type
        }
        hmmer_list <- NULL
        for (gene in gene_list) {
          hmmer_cmd <- paste(hmmer_program, "-A tmp/alignment.sto --domtblout tmp/hmmer.out", "--cpu 4", "--domE", input$eval, tmp_fasta, db, sep = " ")
          system(hmmer_cmd, intern = FALSE)
        }
        # TODO: combine results from multiple genes into list before parser
        withProgress(message = "Performing HMMER search of carge genes in database...", {
          system("python bin/hmm.py -i tmp/hmmer.out -a tmp/alignment.sto")
          # TODO: have parser print to  terminal so that we don't have to load from an intermediate file every time
          read_csv("tmp/hmmer-parsed.csv", show_col_types = FALSE) %>% select(-1)
        })
      } else if (input$input_type == "starship" & input$search_ship_genes == "No") {
        if (query_list$query_type == "nucl") {
          blast_program <- "blastn"
        } else {
          blast_program <- "blastp"
        }
        # TODO: should be separated by results from different genes?
        # TODO: split with conditional statements for HMMER, BLAST, or both
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
            hit_length <- getNodeSet(row, "Iteration_hits//Hit//Hit_len") %>% sapply(., xmlValue)
            bitscore <- getNodeSet(row, "Iteration_hits//Hit//Hit_hsps//Hsp//Hsp_bit-score") %>% sapply(., xmlValue)
            eval <- getNodeSet(row, "Iteration_hits//Hit//Hit_hsps//Hsp//Hsp_evalue") %>% sapply(., xmlValue)
            bind_cols("query_ID" = query_ID, "hit_IDs" = hit_IDs, "hit_length" = hit_length, "bitscore" = bitscore, "eval" = eval)
          })

          # this ensures that NAs get added for no hits
          lapply(results, function(y) {
            as.data.frame((y), stringsAsFactors = FALSE)
          }) %>% bind_rows()
        } else if (input$input_type == "gene" | input$search_ship_genes == "Yes") {
          # TODO: iterate over genes here
          read_csv("tmp/hmmer-parsed.csv",col_names = c("X1", "query_ID", "hit_IDs", "hit_length", "bitscore", "eval"),show_col_types = FALSE) %>% select(-X1)
        }
      }
    })

    output$tbl <- renderDT(
      {
        if (is.null(blastresults())) {
        } else if (input$input_type == "starship") {
          parsedresults() %>% as.data.frame(., stringsAsFactors = FALSE)
        } else if (input$input_type == "gene") {
          parsedresults()
        }
      },
      selection = "single"
    )


    # this chunk gets the alignemnt information from a clicked row
    output$clicked <- renderTable(
      {
        if (is.null(input$tbl_rows_selected)) {} else {
          xmltop <- xmlRoot(blastresults())
          clicked <- input$tbl_rows_selected
          tableout <- parsedresults() %>%
            data.frame() %>%
            slice(clicked)

          tableout <- t(tableout)
          names(tableout) <- c("")
          rownames(tableout) <- c("Query ID", "Hit ID", "Length", "Bit Score", "e-value")
          colnames(tableout) <- NULL
          data.frame(tableout)
        }
      },
      rownames = T,
      colnames = F
    )

    # this chunk makes the alignments for clicked rows
    output$alignment <- renderText({
      if (is.null(input$tbl_rows_selected)) {} else {
        xmltop <- xmlRoot(blastresults())

        clicked <- input$tbl_rows_selected

        # loop over the xml to get the alignments
        align <- xpathApply(blastresults(), "//Iteration", function(row) {
          top <- getNodeSet(row, "Iteration_hits//Hit//Hit_hsps//Hsp//Hsp_qseq") %>% sapply(., xmlValue)
          mid <- getNodeSet(row, "Iteration_hits//Hit//Hit_hsps//Hsp//Hsp_midline") %>% sapply(., xmlValue)
          bottom <- getNodeSet(row, "Iteration_hits//Hit//Hit_hsps//Hsp//Hsp_hseq") %>% sapply(., xmlValue)
          rbind(top, mid, bottom)
        })

        # split the alignments every 40 carachters to get a "wrapped look"
        alignx <- do.call("cbind", align)
        splits <- strsplit(gsub("(.{40})", "\\1,", alignx[1:3, clicked]), ",")

        # paste them together with returns '\n' on the breaks
        split_out <- lapply(1:length(splits[[1]]), function(i) {
          rbind(paste0("Q-", splits[[1]][i], "\n"), paste0("M-", splits[[2]][i], "\n"), paste0("H-", splits[[3]][i], "\n"))
        })
        unlist(split_out)
      }
    })
  })
}
