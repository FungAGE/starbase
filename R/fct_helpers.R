#' helpers
#'
#' @description function for generating an entry form
#'
#' @return The return value, if any, from executing the function.
#'
#' @noRd

# Login
# This function must return a data.frame with columns user and sessionid.  Other columns are also okay
# and will be made available to the app after log in.
get_sessions_from_db <- function(conn = db, expiry = cookie_expiry) {
  dbReadTable(conn, "sessions") %>%
    mutate(login_time = ymd_hms(login_time)) %>%
    as_tibble() %>%
    filter(login_time > now() - days(expiry))
}

# This function must accept two parameters: user and sessionid. It will be called whenever the user
# successfully logs in with a password.
add_session_to_db <- function(user, sessionid, conn = db) {
  tibble(user = user, sessionid = sessionid, login_time = as.character(now())) %>%
    dbWriteTable(conn, "sessions", ., append = TRUE)
}

load_sql <- function(db){
  pool <- pool::dbPool(RSQLite::SQLite(), dbname = db)
  con <- pool::poolCheckout(pool)
  return(con)
}

stop_sql<-function(con){
  pool::poolReturn(con)
  dbDisconnect(con)
}


#Label mandatory fields
labelMandatory <- function(label) {
  tagList(
    label,
    span("*", class = "mandatory_star")
  )
}

# append data to SQL table
sql_append <- function(con,table,data){
  query <- sqlAppendTable(con, table, data, row.names = FALSE)
  dbExecute(con, query)
  shinyalert::shinyalert(title="Starship sequence added to database", type = "success",closeOnClickOutside = TRUE, closeOnEsc = TRUE)
}

sql_update <- function(con,table,data,row){
  sql_names <- colnames(data)
  # Generate the SET part of the SQL query dynamically
  set_clause <- paste(sql_names, "= ?", collapse = ", ")

  rowid <- dbGetQuery(con, sprintf('SELECT "rowid" FROM "%s" WHERE "rowid" = %d', table, row))

  # Construct and execute the SQL update statement
  sql_statement <- sprintf("INSERT OR REPLACE INTO %s(row_id, %s) VALUES (?, %s)", "%s", set_clause, set_clause, table)
  params <- c(as.list(unname(data[1,])), rowid)
  dbExecute(con, sql_statement, params)
}

sql_delete <- function(con,table,row){
  dbExecute(con, sprintf('DELETE FROM "%s" WHERE "rowid" == ("%s")', table, row))
}

# Form for data entry
entry_form <- function(button_id){
  showModal(
    modalDialog(
      div(id=("entry_form"),
          tags$head(tags$style(".modal-dialog{ width:400px}")),
          tags$head(tags$style(HTML(".shiny-split-layout > div {overflow: visible}"))),
          fluidPage(
            fluidRow(
              textInput("uploader",label = labelMandatory("Name of curator")),
              textInput("evidence",label=labelMandatory("What evidence exists for Starship annotation?")),
              textInput("genus",label = labelMandatory("Enter genus name")),
              textInput("species",label = labelMandatory("Enter species name")),
              
              # TODO: store file and put path in SQL table
              fileInput("fna",label = labelMandatory("Upload Starship sequence"),accept = c(".fa",".fna",".fasta")),
              fileInput("gff3",label = labelMandatory("Upload Starship annotations (GFF3 format)"),accept = c(".gff",".gff3")),
              textAreaInput("comment", "Comment", placeholder = "", height = 100, width = "354px"),
              helpText(labelMandatory(""), paste("Mandatory field.")),
              actionButton(button_id, "Submit")
            ),
            easyClose = TRUE
          )
      )
    )
  )
}
    
# save form data into data_frame format
formData <- reactive({
  formData <- data.frame(
    fna = input$fna,
    starshipID = input$starshipID,
    ome = input$ome,
    gff3 = input$gff3,
    genus = input$genus,
    species = input$species,
    strain = input$strain,
    version = input$version,
    biosample = input$biosample,
    assembly_acc = input$assembly_acc,
    acquisition_date = input$acquisition_date,
    genomeSource = input$genomeSource,
    published = input$published,
    clade = input$clade,
    kingdom = input$kingdom,
    phylum = input$phylum,
    subphylum = input$subphylum,
    class = input$class,
    subclass = input$subclass,
    order = input$order,
    family = input$family,
    subfamily = input$subfamily,
    suborder = input$suborder,
    superfamily = input$superfamily,
    tribe = input$tribe,
    evidence = input$evidence,
    source = input$source,
    X.contigID = input$X.contigID,
    captainID = input$captainID,
    elementBegin = input$elementBegin,
    elementEnd = input$elementEnd,
    size = input$size,
    strand = input$strand,
    boundaryType = input$boundaryType,
    emptySiteID = input$emptySiteID,
    emptyContig = input$emptyContig,
    emptyBegin = input$emptyBegin,
    emptyEnd = input$emptyEnd,
    emptySeq = input$emptySeq,
    upDR = input$upDR,
    downDR = input$downDR,
    DRedit = input$DRedit,
    upTIR = input$upTIR,
    downTIR = input$downTIR,
    TIRedit = input$TIRedit,
    nestedInside = input$nestedInside,
    containNested = input$containNested,
    dr = input$dr,
    tir = input$tir,
    starship_family = input$starship_family,
    starship_navis = input$starship_navis,
    starship_haplotype = input$starship_haplotype,
    code = input$code,
    target = input$target,
    spok = input$spok,
    ars = input$ars,
    other = input$other,
    hgt = input$hgt,
    taxid = input$taxid,
    superkingdom = input$superkingdom,
    subkingdom = input$subkingdom,
    checksum = input$checksum,
    rowid = UUIDgenerate(),
    stringsAsFactors = FALSE)
  return(formData)
  
})

# BLAST functions

parse_fasta_input <- function(input,tmp_fasta) {
  if (is.null(tmp_fasta)){
    tmp_fasta <- tempfile(fileext = ".fa")
  }

  if (is.null(input$query_file$datapath) & (input$query_text == "" | is.null(input$query_text))) {
    shinyalert("No input provided:", "Please provide either a query sequence in the text box, or upload a fasta file", type = "error")
  } else if (!is.null(input$query_file$datapath) & (input$query_text != "" | !is.null(input$query_text))) {
    shinyalert("Multiple inputs:", "Please provide either a query sequence in the text box, or upload a fasta file, not both", type = "error")
  } 

  if (!is.null(input$query_text) & input$query_text!="") {
    input_type<-"text"
    query<-input$query_text
  } else if (!is.null(input$query_file$datapath)) {
    input_type<-"file"
    query<-read_lines(input$query_file$datapath)
  } else {
    shinyalert("Error with reading input:", "This was not supposed to happen... Try again?", type = "error")
  }
  
  # throw error if none or both provided
  # choose the input
  # BUG: header is not being added correctly in some cases
  header_grep<-grep(">",query,fixed=TRUE)
  if(length(header_grep)>1) {
    shinyalert("Multi-fastas not supported at this time","Please provide one sequence at a time.",type="error")
  } else if (length(header_grep)==0) {
    query_header<-"QUERY"
  } else if (length(header_grep)==1) {
    query_header <- gsub(">","",query[header_grep])
  } else { 
    shinyalert("Error reading header from input",type="error")
  }
  if (length(query_header)==0){
    query_header<-"QUERY"
  }
  print(paste0("Input type is: ",input_type))
  writeLines(query, tmp_fasta)
  return(list(tmp_file=tmp_fasta,input_type=input_type,query_header=query_header))
}


# TODO: add better error catching for format of query sequence/file here
# function to clean the query sequence: first remove non-letter characters, then ambiguous characters, then guess if query is nucl or protein
clean_lines <- function(query_list) {
  cleaned_seq <- NULL
  nucl_char <- c("A", "T", "G", "C")
  nucl_char <- c(nucl_char,tolower(nucl_char))
  prot_char <- c("A","R","N","P","C","Q","U","G","H","I","L","K","M","F","P","S","T","W","Y","V")
  prot_char <- c(prot_char,tolower(prot_char))
  
  query<-readLines(query_list$tmp_file)

  query_lines <- query[which(!grepl(">",query))]
  for (line in query_lines) {
    if (str_length(line) >= 1) {
      # Remove non-letter characters using gsub
      cleaned_seq <- paste0(cleaned_seq, gsub("[^A-Za-z]", "", trimws(line)),collapse="")
    }
  }

  # count characters for deciding type later
  nucl_count <- sum(table(strsplit(cleaned_seq, ""))[nucl_char], na.rm = TRUE)
  prot_count <- sum(table(strsplit(cleaned_seq, ""))[prot_char], na.rm = TRUE)

  # guess if sequence is nucl or protein
  if (prot_count >= (0.1 * str_length(cleaned_seq))) {
    query_type <- "prot"
    print("Query is protein sequence")
  } else {
    query_type <- "nucl"
    print("Query is nucleotide sequence")
  }
  return(c(query_list,query_type = query_type, cleaned_query = cleaned_seq))
}

run_blast<-function(seq_type=NULL,blast_type=NULL,tmp_fasta=NULL,input.eval=NULL,threads=NULL,stitch=FALSE){
  db_list <- list(
    ship = list(
      nucl = "../Starships/ships/fna/blastdb/concatenated.fa"
    ),
    gene = list(
      tyr = list(
        prot = "../Starships/captain/tyr/faa/blastdb/YRsuperfamRefs.dd.faa",
        nucl = "../Starships/captain/tyr/fna/blastdb/YRsuperfamRefs.fa"
      ),
      fre = list(
        prot = "../Starships/cargo/fre/faa/blastdb/fre.mycoDB.dd.faa",
        nucl = "../Starships/cargo/fre/fna/blastdb/fre.fa"
      ),
      nlr = list(
        prot = "../Starships/cargo/nlr/faa/blastdb/nlr.mycoDB.dd.faa",
        nucl = "../Starships/cargo/nlr/fna/blastdb/nlr.fa"
      ),
      DUF3723 = list(
        prot = "../Starships/cargo/duf3723/faa/blastdb/duf3723.mycoDB.dd.faa",
        nucl = "../Starships/cargo/duf3723/fna/blastdb/duf3723.fa"
      ),
      plp = list(
        prot = "../Starships/cargo/plp/faa/blastdb/plp.mycoDB.dd.faa",
        nucl = "../Starships/cargo/plp/fna/blastdb/plp.fa"
      )
    )
  )
  if(seq_type == "nucl") {
    blast_program = "blastn"
    if (blast_type == "ship" ) {
      blastdb <- db_list[["ship"]][["nucl"]]
    } else {
      blastdb <- db_list[["gene"]][[blast_type]][["nucl"]]
    }
  } else {
    if (blast_type == "ship" ) {
      blast_program = "tblastn"
      blastdb <- db_list[["ship"]][["nucl"]]
    } else {
      blast_program = "blastp"
      blastdb <- db_list[["gene"]][[blast_type]][["prot"]]
    }
  }

  if (length(blastdb) != 1) {
     stop("Issue accessing BLAST database")
  }

  blast_tmp<-tempfile(fileext = ".blast")    
  blast_cmd<-paste(blast_program, "-query", tmp_fasta, "-db", blastdb, "-evalue", input.eval,"-out",blast_tmp,"-outfmt '6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qseq sseq' -max_hsps 1 -max_target_seqs 10 -num_threads",threads,sep=" ")
  system(blast_cmd, intern = FALSE)
  if(stitch == TRUE){
    # BUG: blast_tmp file does not get written with system call?
    stitched_blast_tmp<-tempfile(fileext = ".stitch")
    stitch_blast_cmd<-paste("python bin/BLASTstitcher.py -i",blast_tmp, "-o",stitched_blast_tmp, sep = " ")
    system(stitch_blast_cmd,intern=FALSE)
    ship_blast_out<-stitched_blast_tmp
  } else {
  }
  readr::read_delim(blast_tmp,col_names=c("query_id","hit_IDs","pident","aln_length","mismatches","gaps","query_start","query_end","subject_start","subject_end","evalue","bitscore","query_seq","subject_seq"))
}

run_hmmer <- function(seq_type=NULL,tmp_hmmer=NULL,input.eval=NULL,tmp_fasta=NULL) {
  hmmer_program <- "hmmsearch"
  hmmer_db <- "../Starships/captain/tyr/hmm/YRsuperfams.p1-512.hmm"
  hmmer_cmd <- paste(hmmer_program, "-o", tmp_hmmer, "--cpu 4", "--domE", input.eval, hmmer_db, tmp_fasta, sep = " ")
  system(hmmer_cmd, intern = FALSE)
  tmp_hmmer_parsed <- tempfile(fileext = ".hmmer.parsed.txt")
  system(paste0("python bin/hmm.py --hmmer_output_file ", tmp_hmmer, " --parsed ", tmp_hmmer_parsed),intern=FALSE)
  read_tsv(tmp_hmmer_parsed)
}

# these functions are really just for parsing the results list for genes
make_blast_table<-function(data) {
  if (!is.null(data)) {
    data %>% as.data.frame() %>% 
      dplyr::select(hit_IDs,aln_length,query_start,query_end,gaps,evalue,bitscore) %>%
      DT::datatable(
        options = list(), class = "display", rownames = FALSE,
        callback = JS("return table;"), # rownames, colnames, container,
        caption = NULL, filter = c("none", "bottom", "top"), escape = TRUE,
        style = "auto", width = NULL, height = NULL, elementId = NULL,
        fillContainer = getOption("DT.fillContainer", NULL),
        autoHideNavigation = getOption("DT.autoHideNavigation", NULL),
        selection = "single", extensions = list(),
        plugins = NULL, editable = FALSE)
  }
}

# this function makes the alignments for clicked rows
make_alignment<-function(input,data, name,query_type) {
  #? needed for when row is not clicked?
  clicked<-input[[paste0("table_",name,"_rows_selected")]]
  if (!is.null(clicked)) {
    cdat <- data[input[[paste0("table_",name,"_rows_selected")]],]
    qseq <- cdat$query_seq
    qid <- cdat$query_id
    sseqs <- cdat$subject_seq
    sids <- cdat$hit_IDs

    if (query_type == "nucl" ) {
      ss <- Biostrings::DNAStringSet(c(qseq, sseqs))
      base_palette="nucleotide"
    } else {
      ss <- Biostrings::AAStringSet(c(qseq, sseqs))
      base_palette="clustal"
    }
    # BUG: when setting custom headers for msa object:
    names(ss) <- c(paste0(qid,":"), paste0(sids,":"))        
    msaR::msaR(ss, menu = T, overviewbox = F, seqlogo = FALSE, leftheader = FALSE, labelname = TRUE, labelid = FALSE, labelNameLength = 200,colorscheme=base_palette)
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
        dplyr::select(query_id,hit_IDs,aln_length)
      
    subNames<-unique(df$hit_IDs)
    queryName<-unique(df$query_id)

    groupColors<-RColorBrewer::brewer.pal(nhits,"Paired")
    
    m <- matrix(df$aln_length,
                byrow = FALSE,
                nrow = nhits, ncol = nhits)

    dimnames(m) <- list(query = rep(queryName,nhits),
                        subject = subNames)

    # Build the chord diagram:
    chorddiag::chorddiag(m, 
              groupnameFontsize = 12,
              groupColors = groupColors, 
              groupnamePadding = 50)
}}

render_output <- function(results,input,name,output,query_type, plot_chord) {
  if (name != "hmm") {
    output[[paste0("table_",name)]] <- DT::renderDT({ make_blast_table(results) })
    output[[paste0("alignment_",name)]] <- renderMsaR({ make_alignment(input,results, name,query_type) })
    if (plot_chord) {
      output[[paste0("chord_",name)]] <- chorddiag::renderChorddiag({ make_chord(results) })
    }
  } else {
    superfam<-results %>% slice_min(evalue) %>% pull(hit_IDs)
    output[["superfam_id"]]<-renderValueBox({
      valueBox(
      subtitle = "Likely captain superfamily",
      value = superfam,
      icon = icon("dna")
    )})
    
    output[["superfam_tree"]]<-renderGirafe(captain_tree)
    

    # a custom table container
    sketch <- htmltools::withTags(table(
      class = "display",
      thead(
        tr(
          th(colspan = 4, "hmmersearch of captain genes"),
          th(colspan = 4, "Captain gene phylogeny"),
          th(colspan = 4, "Other Starships in captain superfamily")
        )
      )
    ))

    table_dat<-joined_ships %>%
      filter(!is.na(starship_family)) %>%
      mutate(starship_family=gsub("fam","superfam0",starship_family),
        starship_family = ifelse(grepl("^fam", starship_family) & !is.na(code), code, starship_family))

    reactive({
      req(superfam)

      output[["superfam_table"]]<-renderDT({
        tab_dat<-table_dat %>% filter(starship_family %in% superfam)
        if( nrow(tab_dat)==0) return(NULL)
        
        tab_dat %>%
          DT::datatable(
            options = list(), class = "display", rownames = FALSE, #container = sketch,
            callback = JS("return table;"), 
            caption = NULL, filter = c("none", "bottom", "top"), escape = TRUE,
            style = "auto", width = NULL, height = NULL, elementId = NULL,
            fillContainer = getOption("DT.fillContainer", NULL),
            autoHideNavigation = getOption("DT.autoHideNavigation", NULL),
            selection = "none", extensions = list(),
            plugins = NULL, editable = FALSE
          )
      })
      session$sendCustomMessage(type = "superfam_tree_selected", message = superfam)
    })

    observeEvent(input$reset, {
      session$sendCustomMessage(type = 'plot_set', message = character(0))
    })
      
    }
}