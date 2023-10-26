#' blast_syn_viz UI Function
#'
#' @description A shiny Module.
#'
#' @param id,input,output,session Internal parameters for {shiny}.
#'
#' @noRd
#'
#' @import tidyverse shinyalert vroom shinyWidgets colourpicker waiter purrr
#' @importFrom shiny NS tagList
#' @importFrom tibble deframe
mod_blast_syn_viz_ui <- function(id) {
  ns <- NS(id)
  tabPanel(
    "Main View",
    fluidRow(
      column(
        id = "settingPanel",
        width = 3,
        div(
          class = "boxLike",
          style = "background-color: #FAF9F6;",
          h4(icon("cog"), "Settings"),
          awesomeCheckbox(
            inputId = ns("generateDotPlot"),
            label = "Generate Dot Plot",
            value = TRUE,
            status = "primary"
          ),
          hr(class = "setting"),
          h5("Macro Synteny"),
          radioGroupButtons(
            inputId = ns("macroPlotMode"),
            label = "Choose macro synteny layout",
            choices = c("Circular", "Parallel"),
            selected = "Parallel",
            width = "100%",
            status = "primary"
          ),
          colourInput(
            inputId = ns("macroQueryColor"),
            label = "Query Chr Color",
            value = "#69a3b2"
          ),
          colourInput(
            inputId = ns("macroSubjectColor"),
            label = "Subject Chr Color",
            value = "#B27869"
          ),
          colourInput(
            inputId = ns("macroRibbonColor"),
            label = "Macro Ribbon Color",
            value = "#808080"
          ),
          sliderTextInput(
            inputId = ns("macroChrFontSize"),
            label = "Chromosomes Label Size",
            choices = paste0(seq(0.2, 2, by = 0.1), "rem"),
            selected = "1rem"
          ),
          hr(class = "setting"),
          h5("Micro Synteny"),
          awesomeCheckbox(
            inputId = ns("oneBestSubject"),
            label = HTML("Extract <strong>one best</strong> Subject"),
            value = TRUE,
            status = "primary"
          ),
          colourInput(
            inputId = ns("forwardGeneColor"),
            label = "Forward Gene's Color ",
            value = "#af8dc3"
          ),
          colourInput(
            inputId = ns("reverseGeneColor"),
            label = "Reverse Gene's Color ",
            value = "#7fbf7b"
          ),
          colourInput(
            inputId = ns("microRibbonColor"),
            label = "Micro Ribbon Color",
            value = "#10218b"
          )
        )
      ),
      column(
        id = "mainPanel",
        width = 9,
        div(
          class = "boxLike",
          fluidRow(
            column(
              12,
              h3(icon("file-upload"), "Input")
            )
          ),
          fluidRow(
            column(
              12,
              h5("Hide Setting Panel on the left")
            )
          ),
          fluidRow(
            column(
              12,
              materialSwitch(
                inputId = ns("hide_setting"),
                label = "",
                value = FALSE,
                right = TRUE,
                status = "primary",
                width = NULL
              )
            )
          ),
          fluidRow(
            column(
              12,
              h5("Compare genome to Starships in starbase")
            )
          ),
          fluidRow(
            column(
              6,
              selectizeInput(ns("query_ship_main"), choices = NULL, 
                label = "Query Ship",
                multiple = FALSE,
                width = "100%",
                selected = "Hephaestus"
              )
            ),
            column(
              6,
              selectInput(ns("subject_ship_main"), choices = NULL, 
                label = "Subject Ship",
                multiple = FALSE,
                width = "100%",
                selected = "Pegasus"
              )
            )
          ),
          fluidRow(style="padding-bottom: 15px;",
          column(
              12,
              div(class="float-left",
                  actionButton(
                      inputId = ns("macroSynteny"),
                      status = "secondary",
                      icon = icon("pagelines"),
                      label = "View Macro Synteny"
                  )
                  )
            )
          )
    ),
    icon = icon("binoculars")
  )
  ))
}

#' blast_syn_viz Server Functions
#'
#' @noRd
mod_blast_syn_viz_server <- function(id) {
  moduleServer(id, function(input, output, session) {
    ns <- session$ns

    # populate input for ships
    # BUG: there's one entry with 2 sequences, using distinct to deal with this for now
    ships<-vroom("MTDB/joined_ships.tsv") %>%
      filter(!is.na(fna) & !is.na(gff3)) %>% 
      mutate(ship_code=ifelse(is.na(ship_code),
        paste0(ome,"_s",str_sub(checksum,start=1,end=5)),ship_code)) %>% 
      distinct(ship_code,fna,gff3,genus,species) 
      
    fna_list<-ships %>%
      group_by(ship_code) %>%
      summarise(named_vec = list(fna)) %>%
      deframe()

    gff_list<-ships %>%
      group_by(ship_code) %>%
      summarise(named_vec = list(gff3)) %>%
      deframe()

    taxa_list<-ships %>%
      mutate(taxa=paste(genus,species)) %>%
      group_by(ship_code) %>%
      summarise(named_vec = list(taxa)) %>%
      deframe()

    updateSelectizeInput(session, ns("subject_ship_main"), choices = names(fna_list), server = TRUE)
    updateSelectizeInput(session, ns("query_ship_main"), choices = names(fna_list), server = TRUE)

    if (identical(ns("subject_ship_main"),ns("query_ship_main"))) {
      shinyalert("Oops!", "Please choose different ships to compare", type = "error")
    }

    ## Init synteny data as `reactiveValues`
    # list object for storing reactive values
    synteny <- reactiveValues()

    ## create waiter
    mainView_waiter <- Waiter$new()

    # hide settings
    # observe({
    #   if (input$hide_setting) {
    #     js$hideSetting()
    #   } else {
    #     js$showSetting()
    #   }
    # })

    ## Setup file paths
    # ## query Ship file path
    # queryGffFile <- reactive({gff_list[[ns("query_ship_main")]]})

    queryGffFile <- reactive({"tmp/Alyrata384.bed"})

    # ## subject Ship file path
    # subjectGffFile <- reactive({gff_list[[ns("subject_ship_main")]]})

    subjectGffFile <- reactive({"tmp/Athaliana167.bed"})

    # ## anchor file path
    # anchorFile <- reactive({
    # })

    anchorFile <- reactive({"tmp/Athaliana167.Alyrata384.anchors"})

    # ## anchor lifted file path
    # anchorLiftedFile <- reactive({
    # })
    
    anchorLiftedFile <- reactive({"tmp/Athaliana167.Alyrata384.lifted.anchors"})

    queryGff <- 
    reactive({
      req(queryGffFile())
      if (file.exists(queryGffFile())) {
        # TODO: parse gff
        vroom(queryGffFile(),
              col_names = c("chr", "start","end","gene", "score", "strand")) %>%
          mutate(chr = factor(chr, levels = unique(chr))) %>%
          arrange(chr, start, end)
      } else {
        NULL
      }
    })

    subjectGff <- 
    reactive({
      req(subjectGffFile())
      if (file.exists(subjectGffFile())) {
        vroom(subjectGffFile(),
              col_names = c("chr", "start","end","gene", "score", "strand")) %>%
          mutate(chr = factor(chr, levels = unique(chr))) %>%
          arrange(chr, start, end)
      } else {
        NULL
      }
    })

    # querySpecies <- reactive({taxa_list[[ns("query_ship_main")]]
    # })

    querySpecies <- reactive({"grape"})

    # subjectSpecies <- reactive({taxa_list[[ns("subject_ship_main")]]
    # })

    subjectSpecies <- reactive({"peach"})

    observeEvent(input$macroPlotMode, {
      ## setup different default color scheme for two macro synteny layout
      if (input$macroPlotMode == "Circular") {
        updateColourInput(session, "macroQueryColor", value = "#ca0020")
        updateColourInput(session, "macroSubjectColor", value = "#0571b0")
      } else if (input$macroPlotMode == "Parallel") {
        updateColourInput(session, "macroQueryColor", value = "#69a3b2")
        updateColourInput(session, "macroSubjectColor", value = "#B27869")
      } else {
        updateColourInput(session, "macroQueryColor", value = "#69a3b2")
        updateColourInput(session, "macroSubjectColor", value = "#B27869")
      }
    })

    observeEvent(input$macroSynteny, {
      output$microAnchor_out <- NULL
      shinyjs::hide("microSynteny_download")

      if (is.null(queryGffFile()) || is.null(subjectGffFile())) {
        shinyalert("Oops!", "Query or subject Ship file doesn't exist", type = "error")
      } else if (is.null(anchorFile()) || !file.exists(anchorFile())) {
        shinyalert("Oops!", "Anchor file doesn't exist", type = "error")
      } else if (is.null(anchorLiftedFile()) || !file.exists(anchorLiftedFile())) {
        shinyalert("Oops!", "Anchor lifted file  doesn't exist", type = "error")
      } else {
        ## show spinner
        mainView_waiter$show()
        anchorFile_new <- tempfile()
        system(paste0("awk 'BEGIN{blockID=0}{if($1~/^##/){blockID+=1;}else{print blockID\"\\t\"$0}}' ", anchorFile(), " > ", anchorFile_new))
        anchorNew <- vroom(
          anchorFile_new,
          col_names = c("blockID", "queryGene", "subjectGene", "score")
        ) %>%
          inner_join(queryGff(),
            by = c("queryGene" = "gene"),
            suffix = c(".anchor", ".gff")
          ) %>%
          dplyr::rename(
            "queryChr" = chr,
            "queryStart" = start,
            "queryEnd" = end,
            "queryStrand" = strand
          ) %>%
          select(
            blockID,
            queryGene, queryChr,
            queryStart, queryEnd,
            subjectGene
          ) %>%
          inner_join(subjectGff(),
            by = c("subjectGene" = "gene"),
            suffix = c(".anchor", ".gff")
          ) %>%
          dplyr::rename(
            "subjectChr" = chr,
            "subjectStart" = start,
            "subjectEnd" = end,
            "subjectStrand" = strand
          ) %>%
          select(
            blockID,
            queryGene, queryChr,
            queryStart, queryEnd,
            subjectGene, subjectChr,
            subjectStart, subjectEnd
          )

        get_anchorGene <- function(gene) {
          return(paste(gene, collapse = ","))
        }

        get_spanGenes <- function(chr, start, end, bed) {
          bed %>%
            filter(
              chr == {{ chr }},
              start >= {{ start }},
              end <= {{ end }}
            ) %>%
            arrange(chr, start, end) %>%
            pull(gene)
        }

        anchorSimple <- anchorNew %>%
          group_by(blockID) %>%
          summarize(
            summarized_queryChr = unique(queryChr),
            summarized_queryStart = min(queryStart,na.rm=TRUE),
            summarized_queryEnd = max(queryEnd,na.rm=TRUE),
            summarized_queryAnchorList = get_anchorGene(queryGene),
            summarized_subjectChr = unique(subjectChr),
            summarized_subjectStart = min(subjectStart,na.rm=TRUE),
            summarized_subjectEnd = max(subjectEnd,na.rm=TRUE),
            summarized_subjectAnchorList = get_anchorGene(subjectGene)
          ) %>%
          mutate(
            summarized_queryAnchorList = str_split(summarized_queryAnchorList, ","),
            summarized_subjectAnchorList = str_split(summarized_subjectAnchorList, ",")
          )

        querySpanList <- pmap(
          .l = list(
            anchorSimple$summarized_queryChr,
            anchorSimple$summarized_queryStart,
            anchorSimple$summarized_queryEnd
          ),
          .f = get_spanGenes, queryGff()
        )

        subjectSpanList <- pmap(
          .l = list(
            anchorSimple$summarized_subjectChr,
            anchorSimple$summarized_subjectStart,
            anchorSimple$summarized_subjectEnd
          ),
          .f = get_spanGenes, subjectGff()
        )

        anchorSimple <- anchorSimple %>%
          mutate(
            summarized_querySpanList = querySpanList,
            summarized_subjectSpanList = subjectSpanList
          ) %>%
          mutate(
            aspan = lengths(summarized_querySpanList),
            bspan = lengths(summarized_subjectSpanList),
            asize = lengths(summarized_queryAnchorList),
            bsize = lengths(summarized_subjectAnchorList)
          ) %>%
          mutate(
            a_idx = map2(summarized_queryAnchorList,
              summarized_querySpanList,
              .f = match
            ),
            b_idx = map2(summarized_subjectAnchorList,
              summarized_subjectSpanList,
              .f = match
            )
          ) %>%
          mutate(
            orientation = map2_chr(a_idx, b_idx,
              .f = function(x, y) {
                slope <- coefficients(lm(y ~ x))[2]
                if (slope < 0) {
                  return("-")
                } else {
                  return("+")
                }
              }
            )
          ) %>%
          mutate(
            summarized_q_startGene = map_chr(summarized_querySpanList, 1),
            summarized_q_endGene = map_chr(summarized_querySpanList, .f = function(x) {
              x[length(x)]
            }),
            summarized_s_startGene = map_chr(summarized_subjectSpanList, 1),
            summarized_s_endGene = map_chr(summarized_subjectSpanList, .f = function(x) {
              x[length(x)]
            })
          ) %>%
          dplyr::rename(
            queryChr = summarized_queryChr,
            queryStart = summarized_queryStart,
            queryEnd = summarized_queryEnd,
            q_startGene = summarized_q_startGene,
            q_endGene = summarized_q_endGene,
            subjectChr = summarized_subjectChr,
            subjectStart = summarized_subjectStart,
            subjectEnd = summarized_subjectEnd,
            s_startGene = summarized_s_startGene,
            s_endGene = summarized_s_endGene
          ) %>%
          mutate(
            score = as.integer(sqrt(aspan * bspan))
          ) %>%
          select(
            -contains("List"),
            -c(a_idx, b_idx)
          )

        ## Filter anchor file by minsize and minspan
        minsize <- 0
        minspan <- 30
        anchorSimple <- anchorSimple %>%
          filter(asize >= minsize & bsize >= minsize) %>%
          filter(aspan >= minspan & bspan >= minspan)
        ## write_tsv(anchorSimple, file = paste0(tempdir(), "/1111"))

        ## Process anchor_full
        anchor_full <- vroom(
          anchorLiftedFile(),
          col_names = c("q_Gene", "s_Gene", "score"),
          comment = "#",
          delim = "\t"
        )

        synteny$queryGff <- queryGff() %>%
          filter(chr %in% input$synteny_query_chr) %>%
          arrange(factor(chr, levels = input$synteny_query_chr), start)

        synteny$subjectGff <- subjectGff() %>%
          filter(chr %in% input$synteny_subject_chr) %>%
          arrange(factor(chr, levels = input$synteny_subject_chr), start)

        anchorSimple <- anchorSimple %>%
          filter(queryChr %in% input$synteny_query_chr & subjectChr %in% input$synteny_subject_chr)

        ## Filter and order anchor_full
        synteny$anchor_full <- anchor_full ## %>%
        ## filter(q_GeneChr %in% input$synteny_query_chr) %>%
        ## arrange(q_GeneChr, q_GeneStart)

        summarizeChrInfo <- function(inputBed, chrOrder){
            ## summarize input bed tibble to generate
            ## chr start, end and length
            summarizedBed <- inputBed %>%
                group_by(chr) %>%
                summarise(start = min(start),
                          end = max(end)) %>%
                mutate(chrLength = end - start + 1) %>%
                arrange(match(chr, chrOrder))
            return(summarizedBed)
        }

        queryChrInfo <- summarizeChrInfo(synteny$queryGff, input$synteny_query_chr)
        subjectChrInfo <- summarizeChrInfo(synteny$subjectGff, input$synteny_subject_chr)
        ## Define macro plot mode
        if (input$macroPlotMode == "Circular") {
          plotMode <- "circular"
        } else {
          plotMode <- "parallel"
        }
        ## Generate macro_synteny_data
        macro_synteny_data <- list(
          "querySpecies" = querySpecies(),
          "queryChrInfo" = queryChrInfo,
          "subjectSpecies" = subjectSpecies(),
          "subjectChrInfo" = subjectChrInfo,
          "ribbon" = anchorSimple,
          "plotMode" = plotMode,
          "queryChrColor" = input$macroQueryColor,
          "subjectChrColor" = input$macroSubjectColor,
          "ribbonColor" = input$macroRibbonColor
        )
        session$sendCustomMessage(type = "plotMacroSynteny", macro_synteny_data)

        ## Generate dot_view_data
        if (input$generateDotPlot) {
          anchorSeed <- vroom(anchorFile(),
            comment = "#", delim = "\t",
            col_names = c("queryGene", "subjectGene", "mcscan_score")
          ) %>%
            inner_join(synteny$queryGff, by = c("queryGene" = "gene")) %>%
            inner_join(synteny$subjectGff,
              by = c("subjectGene" = "gene"),
              suffix = c("_query", "_subject")
            )

          dot_view_data <- list(
            "querySpecies" = querySpecies(),
            "queryChrInfo" = queryChrInfo,
            "subjectSpecies" = subjectSpecies(),
            "subjectChrInfo" = subjectChrInfo,
            "anchorSeed" = anchorSeed
          )

          session$sendCustomMessage(type = "plotDotView", dot_view_data)
          shinyjs::show("dotView_download")
        }
        ## hide spinner
        mainView_waiter$hide()
      }
    })


    observeEvent(input$selected_macroRegion, {
      ## The start and end gene were from 5' to 3' order on the genome
      ## no matter the orientation relationship between query and subject region

      macroQueryChr <- input$selected_macroRegion[["macroQueryChr"]]
      macroQueryStart <- input$selected_macroRegion[["macroQueryStart"]]
      macroQueryEnd <- input$selected_macroRegion[["macroQueryEnd"]]
      macroSubjectChr <- input$selected_macroRegion[["macroSubjectChr"]]
      macroSubjectStart <- input$selected_macroRegion[["macroSubjectStart"]]
      macroSubjectEnd <- input$selected_macroRegion[["macroSubjectEnd"]]

      synteny$selectedQueryRegion <- synteny$queryGff %>%
        filter(
          chr == macroQueryChr,
          start >= macroQueryStart,
          end <= macroQueryEnd
        )

      synteny$selectedSubjectRegion <- synteny$subjectGff %>%
        filter(
          chr == macroSubjectChr,
          start >= macroSubjectStart,
          end <= macroSubjectEnd
        )

      synteny$selectedAnchors <- synteny$anchor_full %>%
        filter(
          q_Gene %in% synteny$selectedQueryRegion$gene,
          s_Gene %in% synteny$selectedSubjectRegion$gene
        )

      ## firstly filter anchors, retain the one with highest score value
      ## this behaviour should be the same as generating i1 blocks in mcscan doc
      if (input$oneBestSubject) {
        synteny$selectedAnchors <- synteny$selectedAnchors %>%
          group_by(q_Gene) %>%
          summarise(s_Gene = first(s_Gene, order_by = desc(score)))
      }

      ## Added genomics coordinates to anchors
      synteny$selectedAnchors <- synteny$selectedAnchors %>%
        full_join(synteny$selectedQueryRegion,
          by = c("q_Gene" = "gene"),
          suffix = c(".anchor", ".gff")
        ) %>%
        select(
          q_Gene,
          chr, start, end, strand,
          s_Gene
        ) %>%
        dplyr::rename(
          q_GeneChr = chr,
          q_GeneStart = start,
          q_GeneEnd = end,
          q_GeneStrand = strand
        ) %>%
        full_join(synteny$selectedSubjectRegion,
          by = c("s_Gene" = "gene"),
          suffix = c(".anchor", ".gff")
        ) %>%
        select(
          q_Gene,
          q_GeneChr, q_GeneStart, q_GeneEnd, q_GeneStrand,
          s_Gene,
          chr, start, end, strand
        ) %>%
        dplyr::rename(
          s_GeneChr = chr,
          s_GeneStart = start,
          s_GeneEnd = end,
          s_GeneStrand = strand
        ) %>%
        arrange(q_GeneChr, q_GeneStart, q_GeneEnd, s_GeneChr, s_GeneStart, s_GeneEnd)

      ## Put all anchor infor into result table
      output$microAnchor_out <- DT::renderDataTable(
        {
          synteny$selectedAnchors
        },
        selection = "single",
        rownames = FALSE,
        server = TRUE
      )

      micro_synteny_data <- list(
        "microQueryRegion" = synteny$selectedQueryRegion,
        "microSubjectRegion" = synteny$selectedSubjectRegion,
        "microAnchors" = synteny$selectedAnchors,
        "microForwardColor" = input$forwardGeneColor,
        "microReverseColor" = input$reverseGeneColor,
        "microRibbonColor" = input$microRibbonColor
      )

      session$sendCustomMessage(type = "plotSelectedMicroSynteny", micro_synteny_data)

      shinyjs::show("microSynteny_download")
    })

    observeEvent(input$microAnchor_out_rows_selected, {
      selectedQueryGene <- synteny$selectedAnchors[input$microAnchor_out_rows_selected, ] %>%
        pull(q_Gene)
      session$sendCustomMessage(type = "center_microSynteny", selectedQueryGene)
    })

    observeEvent(input$selected_anchors, {
      ## shiny transferred data are nested lists
      ## not friendly
      output$selected_anchors <- DT::renderDataTable(
        {
          input$selected_anchors %>%
            as_tibble() %>%
            dplyr::rename(
              "QueryGene" = queryGene,
              "QueryChr" = chr_query,
              "QueryStart" = start_query,
              "QueryEnd" = end_query,
              "QueryStrand" = strand_query,
              "SubjectGene" = subjectGene,
              "SubjectChr" = chr_subject,
              "SubjectStart" = start_subject,
              "SubjectEnd" = end_subject,
              "SubjectStrand" = strand_subject
            ) %>%
            unnest(cols = c(QueryGene, QueryChr, QueryStart, QueryEnd, QueryStrand, SubjectGene, SubjectChr, SubjectStart, SubjectEnd, SubjectStrand))
        },
        selection = "single",
        rownames = FALSE,
        server = TRUE
      )

      output$dotviewTable <- renderUI({
        tagList(
          h4("Anchor Genes:"),
          p(
            style = "color: gray;",
            "Please select your region of interest from the dot plot on the left panel, the table below will be updated automatically."
          ),
          DTOutput("selected_anchors")
        )
      })
    })

    observeEvent(input$macroChrFontSize, {
      session$sendCustomMessage(type = "updateChrFontSize", input$macroChrFontSize)
    })
  })
}