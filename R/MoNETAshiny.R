#' MoNETA Shiny app
#'
#' @import shiny
#' @import shinydashboard
#' @import shinyalert
#' @import shinyFiles
#' @import shinyMatrix
#' @import readr
#' @import glue
#' @import conflicted
#' @import fpc
#' @import plotly
#' @import ggplot2
#' @import visNetwork
#' @import tidyverse
#' @import dendextend
#' @import magrittr
#' @import dplyr
#' @import igraph
#' @return Shiny app
#' @export


MoNETAshiny = function() {
    options(shiny.maxRequestSize = 10000 * 1024^2)
    shiny::shinyApp(ui, server)

}


netSummary <- function(network){
    # -------------------------------------------------------------------------
    #  This function outputs a summary for a given network, like the network
    #  dimension, the number of genes and edges, the nodes with the highest degree
    #  and their degree.
    # -------------------------------------------------------------------------
    nodes <- c(network$source, network$dest)
    degree <- table(nodes) %>%
        tidyr::as_tibble(.) %>%
        dplyr::arrange(desc(n))
    top10 <- c(degree[c(1:10), 'nodes'])
    netSummary_list <- list( 'dim'= dim(network),
                             'edges'= dim(network)[1],
                             'nodes'= length(unique(nodes)),
                             'maxConnections'= degree[[1, 'n']],
                             'maxDegreeNodes'= top10$nodes
    )
    return(netSummary_list)
}


ui <- shinydashboard ::dashboardPage(
    shinydashboard::dashboardHeader(title = shiny::span("MoNETA ", style = "color: white; font-size: 28px")),

    shinydashboard::dashboardSidebar(
        shinydashboard::sidebarMenu(
            shiny::tags$head(shiny::tags$link(rel = "stylesheet", type = "text/css", href = "custom.css")),
            style = "position: fixed; height: 90vh; overflow-y: auto;",
            shinydashboard::menuItem("Omics Data", tabName = "mat_tab", icon = shiny::icon("table"),
                                     shinydashboard::menuItem("Loading", tabName = "mat_sub_1"),
                                     shinydashboard::menuItem("Pre-processing", tabName = "mat_sub_2")),
            shinydashboard::menuItem("Network", tabName = "net_tab", icon = shiny::icon("circle-nodes"),
                                     shinydashboard::menuItem("Network inference", tabName = "net_sub_1"),
                                     shinydashboard::menuItem("Network Filtering", tabName =  "net_sub_2")),
            shinydashboard::menuItem("RWR", tabName = "rwr_tab", icon = shiny::icon("person-running"),
                                     shinydashboard::menuItem("Loading networks", tabName = "rwr_tab_1"),
                                     shinydashboard::menuItem("Multiplex network", tabName = "rwr_tab_2"),
                                     shinydashboard::menuItem("RWR Parameters", tabName = "rwr_tab_3")),
            shinydashboard::menuItem("Dimensionality Reduction", tabName = "dr_tab", icon = shiny::icon("filter")),
            shinydashboard::menuItem("Clustering", tabName = "cl_tab", icon = shiny::icon("layer-group"))

        )),

    shinydashboard::dashboardBody(
        shinydashboard::tabItems(
            shinydashboard::tabItem(tabName = "mat_sub_1",
                                    shiny::fluidRow(
                                        shinydashboard::box(
                                            shiny::h4(shiny::span("Upload omics matrices", style = "font-weight: bold")),
                                            shiny::radioButtons(inputId = 'omics_example_opt', label = 'Load the example dataset',
                                                                choices = c('Yes', 'No'), selected = 'No'),
                                            conditionalPanel(
                                                condition = 'input.omics_example_opt == "Yes"',
                                                shiny::checkboxGroupInput(inputId = 'omics_example_files', label = 'Select one or more omics matrices',
                                                                          choices = c("GliomaCNV_norm","GliomaMethylation_norm", "GliomaExpression_norm" ),
                                                                          selected = c("GliomaCNV_norm","GliomaMethylation_norm", "GliomaExpression_norm" ))
                                                             ),
                                            conditionalPanel(
                                                condition = 'input.omics_example_opt == "No"',
                                                shiny::fileInput(inputId = "omics_files", label = 'Select one or more files' ,
                                                                 multiple = TRUE),
                                                shiny::uiOutput('omics_names')
                                            ),
                                            shiny::actionButton('load_mat_button', 'Load')
                                        ),
                                        shiny::htmlOutput(outputId = "omics_sum_info")
                                    ),
                                    shiny::fluidRow(
                                        shinydashboard::box(
                                            shiny::h4(shiny::span("Upload annotation file", style = "font-weight: bold")),
                                            shiny::radioButtons(inputId = 'anno_example_opt', label = 'Load the example annotation file',
                                                                choices = c('Yes', 'No'), selected = 'No'),
                                            conditionalPanel(
                                                condition = 'input.anno_example_opt == "No"',
                                                 shiny::fileInput(inputId = "anno_file", label = 'Select a file',
                                                             multiple = FALSE)
                                            ),
                                            shiny::actionButton('load_anno_button', 'Load')
                                        ),
                                        shiny::htmlOutput(outputId = "anno_info"),
                                        shiny::htmlOutput(outputId = "anno_info_extra")
                                    )
            ),
            shinydashboard::tabItem(tabName = "mat_sub_2",
                                    shiny::h2(shiny::span("Omics matrices Pre-processing", style = "font-weight: bold")),
                                    shiny::fluidRow(
                                        shiny::uiOutput("process_omics_mat"),
                                        shiny::uiOutput('pro_matrices_sum_info'),
                                        shiny::uiOutput('download_proc_mat_box')
                                    ),
                                    shiny::fluidRow(
                                        shiny::uiOutput('intersection'),
                                        shiny::uiOutput('intersection_info')
                                    )
            ),

            shinydashboard::tabItem(tabName = 'net_sub_1',
                                    shiny::h2(shiny::span("Omics Network Inference", style = "font-weight: bold")),
                                    shiny::fluidRow(
                                        shiny::uiOutput("omics_net_arguments"),
                                        shiny::uiOutput("download_net_box")
                                    )
            ),
            shinydashboard::tabItem(tabName = 'net_sub_2',
                                    shiny::h2(shiny::span("Omics Network Filtering", style = "font-weight: bold")),
                                    shiny::fluidRow(
                                        shiny::uiOutput("plot_net_box"),
                                        shiny::uiOutput("f_net_sum_info"),
                                        shiny::uiOutput("download_fnet_box")
                                    )
            ),

            shinydashboard::tabItem(tabName = "rwr_tab_1",
                                    shiny::fluidRow(
                                        shinydashboard::box(
                                            shiny::fileInput("omics_net_files", label = shiny::h4(shiny::span("Upload omics networks", style = "font-weight: bold")),
                                                             multiple = TRUE,
                                                             accept = c("text/csv", '.RDS', "text/comma-separated-values,text/plain", ".csv")
                                            ),
                                            shiny::uiOutput('omics_net_names'),
                                            shiny::actionButton('load_omics_net_button', 'Load')
                                        ),
                                        shiny::htmlOutput(outputId = "omics_net_sum_info")
                                    ),
                                    shiny::fluidRow(
                                        shinydashboard::box(
                                            shiny::fileInput(inputId = "anno_file1",
                                                             label = shiny::h4(shiny::span("Upload annotation file", style = "font-weight: bold")),
                                                             multiple = FALSE),
                                            shiny::actionButton('load_anno_button1', 'Load')
                                        ),
                                        shiny::htmlOutput(outputId = "anno_info1"),
                                        shiny::htmlOutput(outputId = "anno_info_extra1"),
                                        shiny::htmlOutput(outputId = "anno_info_extra2")
                                    )
            ),
            shinydashboard::tabItem(tabName = "rwr_tab_2",
                                    shiny::h2(shiny::span("Construction of Multiplex network", style = "font-weight: bold")),
                                    shiny::uiOutput("multiplex_func")
            ),
            shinydashboard::tabItem(tabName = "rwr_tab_3",
                                    shiny::h2(shiny::span("Parameters setting for RWR-(M)", style = "font-weight: bold")),
                                    shiny::fluidRow(
                                        shinydashboard::box(width = 12,
                                                            shiny::tabPanel(shiny::h4(shiny::span('Parameters', style = "font-weight: bold")),
                                                                            shiny::fluidRow(
                                                                                shiny::column( width = 4,
                                                                                               shiny::h4(shiny::span('General Parameters', style = "font-weight: bold")),
                                                                                               shiny::sliderInput('restart', 'r', min = 0, max = 1, value = 0.7, step = 0.01),
                                                                                               shiny::sliderInput('delta',  shiny::HTML("&delta;"),  min = 0, max = 1, value = 0.5, step = 0.01)
                                                                                ),
                                                                                shiny::column(width = 4,
                                                                                              shiny::h4(shiny::span('Restarting Probabilities', style = "font-weight: bold")),
                                                                                              shiny::uiOutput(outputId = "tauBIO")
                                                                                ),
                                                                                shiny::column(width = 4,
                                                                                              shiny::h4(shiny::span('RWR-NF', style = "font-weight: bold")),
                                                                                              shiny::radioButtons('omics_jump_neigh', 'Connect nodes of different layers by neighborhood', c('YES', 'NO'), selected = 'NO'),
                                                                                              shiny::conditionalPanel(
                                                                                                  condition = 'input.omics_jump_neigh == "YES"',
                                                                                                  shiny::radioButtons('omics_weight_jump_neigh', 'Weight inter-layers edges', c('YES', 'NO'), selected = 'YES')
                                                                                              )
                                                                                )
                                                                            ),
                                                                            shiny::fluidRow(
                                                                                shiny::column(width = 8,
                                                                                              shiny::hr(),
                                                                                              shiny::h4(shiny::span('Transition layer matrices', style = "font-weight: bold")),
                                                                                              shiny::radioButtons(inputId = 'bioInf_transition',
                                                                                                                  label = 'Do you want to create a biologically informed layer transition matrix?',
                                                                                                                  choices = c('YES', 'NO'), selected = 'NO'),
                                                                                              shiny::conditionalPanel(
                                                                                                  condition = 'input.bioInf_transition == "NO"',
                                                                                                  shiny::uiOutput('omics_trans_mat')
                                                                                              )
                                                                                ),
                                                                                shiny::column(width = 4,
                                                                                              shiny::hr(),
                                                                                              shiny::h4(shiny::span('Cores', style = "font-weight: bold")),
                                                                                              shiny::numericInput('cores_rwr', 'Select the number of cores to run the RWR ', min = 1, max = 190, value = 1)
                                                                                )
                                                                            ),
                                                                            shiny::fluidRow(
                                                                                shiny::tags$hr(),
                                                                                shiny::column(width = 6, shiny::actionButton('rwr_button', 'Run RWR')),
                                                                                shiny::column(width = 6, shiny::uiOutput('download_rwr_mat')),
                                                                            )
                                                            )
                                        )
                                    )
            ),

            shinydashboard::tabItem(tabName = "dr_tab",
                                    shiny::h2(shiny::span("Dimensionality Reduction Methods", style = "font-weight: bold")),
                                    shiny::fluidRow(
                                        shiny::uiOutput("dr_opt"),
                                        shiny::uiOutput('dr_plot_box'),
                                        shiny::uiOutput('download_dr_box')
                                    ),
            ),

            shinydashboard::tabItem(tabName = "cl_tab",
                                    shiny::h2(shiny::span("Clustering Analysis", style = "font-weight: bold")),
                                    shiny::fluidRow(
                                        shiny::uiOutput("cl_opt"),
                                        shiny::uiOutput(outputId = "cl_plot_box"), #%>% withSpinner(. , image = spinner_image)
                                    ),
                                    shiny::fluidRow(
                                        shiny::uiOutput('cluster_table_box'),
                                        shiny::conditionalPanel(
                                            condition = "input.cluster_method == 'hclust'",
                                            shiny::uiOutput(outputId = "plotHclust_box")
                                        )
                                    )

            )
        )
    ),

    shiny::tags$head(shiny::tags$style(shiny::HTML('* {font-family: "LM Roman 10"};')))

)


server <- function(input, output, session) {

    ############################################################################################
    #                            OMICS DATA: LOADING DATA (INPUT)                       #
    ############################################################################################

    ############## inputs loading #############

    ### Load omics matrices
    omics_files <- shiny::eventReactive(input$load_mat_button, {
        omicsFiles <- list()

        if (input$omics_example_opt == 'No'){
            inFiles <- input$omics_files
            if (is.null(inFiles))
                return(NULL)
            for (i in 1:nrow(inFiles)) {
                name <- tools::file_path_sans_ext(inFiles$name[i])
                input_name <- input[[paste0('name', i)]]
                new_name <- ifelse(nchar(input_name) == 0, name, input_name)

                ext <- tools::file_ext(inFiles$name[i])
                file <- switch(ext,
                               csv = vroom::vroom(inFiles$datapath[i], delim = ","),
                               tsv = vroom::vroom(inFiles$datapath[i], delim = "\t"),
                               RDS = readRDS(inFiles$datapath[i]),
                               validate("Invalid file; Please upload a .csv, .tsv or .RDS file")
                )
                if (is.numeric(file)){
                    sorted_file <- file[, sort(colnames(file))]
                    omicsFiles[[new_name]] <- sorted_file
                }else{
                    shinyalert::shinyalert("Type Error", "Uploaded Data is not a numeric matrix", type = "error")
                    returnValue()
                }
            }

        } else {
            req(input$omics_example_files)
            example_data <- MoNETA::GBM_mtx[input$omics_example_files]
            omicsFiles <- example_data
        }

        omicsFiles
    })

    ### Load annotation file
    annotation <- shiny::eventReactive(input$load_anno_button, {
        listFiles <- list()

        if (input$anno_example_opt == 'No'){
            inFiles <- input$anno_file
            if (is.null(inFiles)){
                return(NULL)
            }else{
                ext <- tools::file_ext(inFiles$name)
                file <- switch(ext,
                               csv = vroom::vroom(inFiles$datapath, delim = ","),
                               RDS = readRDS(inFiles$datapath),
                               validate("Invalid file; Please upload a .csv or .RDS file")
                )
                if (is.data.frame(file)){
                    listFiles[['annotation']] <- file
                    return(listFiles)
                }else {
                    shinyalert::shinyalert("Type Error", "Uploaded Data is not a dataframe", type = "error")
                    returnValue()
                }
            }
        }else{
           listFiles[['annotation']] <- MoNETA::GBM_pdata
           return(listFiles)
        }
    })


    ############## set files names #############

    output$omics_names <- shiny::renderUI({
        omicsNames <- list()
        inFiles <- input$omics_files
        if (is.null(inFiles))
            return(NULL)
        for (i in 1:nrow(inFiles)) {
            name <- tools::file_path_sans_ext(inFiles$name[i])
            omicsNames[[i]] <- shiny::textInput(inputId = paste0('name', i),
                                                label = paste('Rename the ', name, ' file'),
                                                value = NULL)
        }
        omicsNames
    })


    ############################################################################################
    #                        OMICS DATA: PRE-PROCESSING (INPUT)                      #
    ############################################################################################

    intersected_omics_mat <- shiny::reactiveValues()
    shiny::observeEvent(input$intersect_btn, {
        processed_mats <- shiny::reactiveValuesToList(proc_matrices)
        #shiny::isolate({
        if (input$intersect_opt == 'YES'){
            intersected_omics_mat$matrices <- get_intersection_matrices(processed_mats)
        }else {
            intersected_omics_mat$matrices <- processed_mats
        }
        #})
    })


    ############################################################################################
    #                        NETWORK : NETWORK INFERENCE (INPUT)                      #
    ############################################################################################

    ############################################################################################
    #                        NETWORK : NETWORK FILTERING (INPUT)                      #
    ############################################################################################

    ##### See  NETWORK INFERENCE: NETWORK FILTERING (OUTPUT)



    ############################################################################################
    #                                 RWR: LOADING NETWORKS (INPUT)                            #
    ############################################################################################

    #################### load networks #####################

    loaded_omics_net_list <- shiny::eventReactive(input$load_omics_net_button, {
        req(input$omics_net_files)
        listFiles <- list()
        inFiles <- input$omics_net_files
        if (is.null(inFiles)){
            return(NULL)
        } else {
            for (i in 1:nrow(inFiles)) {
                name <- tools::file_path_sans_ext(inFiles$name[i])
                input_name <- input[[paste0('omics_net_name', i)]]
                new_name <- ifelse(nchar(input_name) == 0, name, input_name)

                ext <- tools::file_ext(inFiles$name[i])
                file <- switch(ext,
                               csv = vroom::vroom(inFiles$datapath[i], delim = ","),
                               RDS = readRDS(inFiles$datapath[i]),
                               validate("Invalid file; Please upload a .csv or .RDS file")
                )
                if (is.data.frame(file) & ncol(file) == 3){
                    listFiles[[new_name]] <- file
                }else if (!is.data.frame(file) ){
                    shinyalert::shinyalert("Type Error", "Uploaded Data is not a dataframe", type = "error")
                    returnValue()
                }else if (ncol(file) != 3){
                    shinyalert::shinyalert("Column Error", "Uploaded Data has not 3 columns", type = "error")
                    returnValue()
                }
            }
        }
        return(listFiles)
    })

    ############## set files names #############

    output$omics_net_names <- shiny::renderUI({
        inFiles <- input$omics_net_files
        if (is.null(inFiles)){
            return(NULL)
        }else{
            omicsNames <- list()
            for (i in 1:nrow(inFiles)) {
                name <- tools::file_path_sans_ext(inFiles$name[i])
                omicsNames[[i]] <- shiny::textInput(paste0('omics_net_name', i), paste('Rename the ', name, ' file'))
            }
            omicsNames
        }
    })

    ### annotation1
    ### Load annotation file

    annotation1 <- shiny::eventReactive(input$load_anno_button1, {
        req(input$anno_file1)
        listFiles <- list()
        inFiles <- input$anno_file1
        if (is.null(inFiles)){
            return(NULL)
        } else {
            ext <- tools::file_ext(inFiles$name)
            new_name <- 'annotation1'
            file <- switch(ext,
                           csv = vroom::vroom(inFiles$datapath, delim = ","),
                           RDS = readRDS(inFiles$datapath),
                           validate("Invalid file; Please upload a .csv or .RDS file")
            )
            if (is.data.frame(file)){
                listFiles[[new_name]] <- file
            }else {
                shinyalert::shinyalert("Type Error", "Uploaded Data is not a dataframe", type = "error")
                returnValue()
            }
        }
        return(listFiles)
    })


    ############################################################################################
    #                                RWR: MULTIPLEX NETWORK (INPUT)                            #
    ############################################################################################

    g_net_list <- shiny::reactive({
        if (length(shiny::reactiveValuesToList(filteredOmicsNetworks)) != 0) {
            return(shiny::reactiveValuesToList(filteredOmicsNetworks))
        }else{
            return(NULL)
        }
    })

    l_net_list <-  shiny::reactive({
        if (!is.null(loaded_omics_net_list())) {
            return(loaded_omics_net_list())
        } else {
            return(NULL)
        }
    })

    multiplex_network <- shiny::reactiveValues()
    shiny::observeEvent(input$gen_multiplex_btn, {
        print('Generating the multiplex network...')
        net_list <- if (!is.null(g_net_list())) g_net_list() else l_net_list()
        multiplex <- create_multiplex(net_list, weighted = ifelse(input$weightMultiplex == 'YES', TRUE, FALSE))
        multiplex_network$multiplex <- multiplex %>% dplyr::ungroup(.)
        if(input$weightMultiplex == 'NO')
            multiplex_network$pruned_multiplex <- multiplex
    })
    shiny::observeEvent(input$prune_multiplex_btn, {
        print('Pruning the multiplex network...')
        multiplex <- multiplex_network$multiplex
        if(input$pruneMultiplex == 'YES')
            multiplex_network$pruned_multiplex <- prune_multiplex_network(multiplex, input$pruneMultiplex_th)
        else
            multiplex_network$pruned_multiplex <- multiplex
    })


    ############################################################################################
    #                                   RWR: PARAMETERS (INPUT)                                #
    ############################################################################################

    ######################## RWR-SIM ###########################

    RWR_output <- shiny::reactiveValues()
    shiny::observeEvent(input$rwr_button, {
        if (is.null(shiny::reactiveValuesToList(multiplex_network)$pruned_multiplex)){
            return(NULL)
        } else{
            req(multiplex_network)
            multiplex_list <- shiny::reactiveValuesToList(multiplex_network)
            omics_multiplex <- if (is.null(multiplex_list$pruned_multiplex)) multiplex_list$multiplex else multiplex_list$pruned_multiplex

            if (input$bioInf_transition == 'YES') {
                net_list <- if (!is.null(g_net_list())) g_net_list() else l_net_list()
                trans_mat <- create_layer_transition_matrix(net_list)
            }else {
                trans_mat <- input$editable_trans_mat
            }

            progress <- shiny::Progress$new()
            on.exit(progress$close())

            progress$set(message = "MoNETA", detail = paste("Doing RWR"), value = 0)

            rwr_simMat <- gen_sim_mat_M(network = omics_multiplex,
                                        tau = unlist(omics_tau_list()),
                                        restart = input$restart,
                                        delta = input$delta,
                                        layer_transition = trans_mat,
                                        jump_neighborhood = ifelse(input$omics_jump_neigh == 'YES', TRUE, FALSE),
                                        weighted_multiplex  = ifelse(input$omics_weight_jump_neigh == 'YES', TRUE, FALSE),
                                        cores = input$cores_rwr)

            progress$inc(0.8, detail = paste("Creating summary file!"))

            Summary <- list(
                'parameters' = c(
                    '################################################################################',
                    '#                             SUMMARY: RWR-M parameters                        #',
                    '################################################################################',
                    '', '',
                    paste('- Number of layers:                     ', length(unique(omics_multiplex$EdgeType))), '',
                    paste('- Restart probability (r):              ', input$restart), '',
                    paste('- Cross-jumping probability (delta):    ', input$delta), '',
                    paste('- Restart probability per layer (tau):  ', paste(omics_tau_list(), collapse = ', ')), '',
                    paste('- Weighted Multiplex Network:           ', input$weightMultiplex),'',
                    paste('- Jumping Neighborhood permission:      ', input$omics_jump_neigh),'',
                    paste('- Weighted jump_neigh edges:            ', input$omics_weight_jump_neigh),''),
                'Transition layers matrix' = trans_mat
            )

            RWR_output$rwr_simMat <- rwr_simMat
            RWR_output$Summary <- Summary
        }
    })

    ######################## download RWR-SIM mat ###########################

    output$download_rwr_mat <- shiny::renderUI({
        if (length(shiny::reactiveValuesToList(RWR_output)) == 0) {
            return(NULL)
        } else {
            shiny::downloadButton(
                outputId = 'download_rwr_btn',
                label = "Download",
                icon = shiny::icon("download")
            )
        }
    })

    output$download_rwr_btn <- shiny::downloadHandler(
        filename = function() {
            paste("MoNETAshiny_RWR_mat_", Sys.Date(), ".zip", sep = "")
        },
        content = function(file){
            temp_directory <- file.path(tempdir(), as.integer(Sys.time()))
            dir.create(temp_directory)

            to_download <- shiny::reactiveValuesToList(RWR_output)
            for (obj in names(to_download)) {
                if (obj == 'Summary'){
                    sink(file = glue(temp_directory, "/{obj}.txt"), append =  TRUE)
                    sink(NULL)
                } else {
                    file_name <- "RWR_sim_mat.csv"
                    readr::write_csv(as.data.frame(to_download[[obj]]), file.path(temp_directory, file_name))
                }
            }
            zip::zip(
                zipfile = file,
                files = dir(temp_directory),
                root = temp_directory
            )
        }
    )


    ############################################################################################
    #                                DIMENSIONALITY REDUCTION (INPUT)                          #
    ############################################################################################

    ######################## DR matrices ###########################

    to_create <-  shiny::reactiveVal(NULL)

    dr_output <- shiny::reactiveValues()
    shiny::observeEvent(input$dr_btn, {
        if (length(shiny::reactiveValuesToList(RWR_output)) == 0){
            return(NULL)
        }else{
            to_create(NULL)
            rwr_mat <- shiny::reactiveValuesToList(RWR_output)$rwr_simMat
            req(input$dr_method)
            progress <- shiny::Progress$new()
            on.exit(progress$close())
            progress$set(message = "MoNETA", detail = paste("Doing Embedding"), value = 0)

            if (input$emb_opt == 'YES'){
                emb_mat <- get_embedding(matrix = rwr_mat, embedding_size = input$dim_emb, cores = input$cores_emb)
                dr_output$emb_mat <- emb_mat
            }else{
                emb_mat <- rwr_mat
            }

            if (input$dr_method == 'UMAP'){
                progress$inc(0.8, detail = paste("Doing UMAP"))
                dr_mat <- get_parallel_umap_embedding(matrix = emb_mat,
                                                      embedding_size = 2,
                                                      n_threads = input$threads)
            }else if (input$dr_method == 'PCA'){
                progress$inc(0.8, detail = paste("Doing PCA"))
                dr_mat <- get_pca_embedding(matrix = emb_mat, embedding_size = input$emb_size)
            }else if (input$dr_method == 'tSNE'){
                progress$inc(0.8, detail = paste("Doing tSNE"))
                dr_mat <- get_tsne_embedding(matrix = emb_mat,
                                             embedding_size = input$emb_size_tsne,
                                             perplexity = input$perplexity, max_iter = input$max_iter,
                                             num_threads = input$threads)
            }

            progress$inc(1, detail = paste("Done!"))
            dr_output$dr_mat <- dr_mat
        }
    })


    dr_plot <- shiny::reactive({
        method <- ifelse(input$dr_method == 'UMAP', 'UMAP of',
                         ifelse(input$dr_method == 'PCA', 'PCA of', 'tSNE of'))
        emb <- ifelse(input$emb_opt == 'YES', 'the Embedded RWR matrix', 'the RWR matrix')
        dr_mat <- shiny::reactiveValuesToList(dr_output)$dr_mat
        dr_plot <- MoNETA::plot_2D_matrix(coord = dr_mat[1:2,], nodes_anno = data.frame(id = colnames(dr_mat)),
                                          id_name = 'id', interactive = FALSE, wo_legend = FALSE) +
            ggtitle(paste(method, emb)) +
            xlab(paste0(input$dr_method, '1')) + ylab( paste0(input$dr_method, '2')) +
            theme(text = ggplot2::element_text(family="Times New Roman"),
                  plot.title = ggplot2::element_text(size = 16, hjust = 0.5, face = 'bold'),
                  axis.title.x = ggplot2::element_text(size = 10),
                  axis.title.y = ggplot2::element_text(size = 10))
        return(dr_plot)
    })


    ############################################################################################
    #                                CLUSTERING: LOADING DATA (INPUT)                          #
    ############################################################################################

    anno <-  shiny::reactive({
        if (!is.null(annotation())) {
            return(annotation()$annotation)
        } else {
            return(NULL)
        }
    })

    anno1 <- shiny::reactive({
        if (!is.null(annotation1())) {
            return(annotation1()$annotation1)
        }else{
            return(NULL)
        }
    })


    clusters <- shiny::reactiveValues()
    shiny::observeEvent(input$cl_btn, {
        if (length(shiny::reactiveValuesToList(dr_output)) == 0){
            return(NULL)
        }else{
            if (input$cl_on_emb == 'YES'){
                mat <- shiny::reactiveValuesToList(dr_output)$emb_mat
            } else {
                mat <- shiny::reactiveValuesToList(dr_output)$dr_mat
            }

            clusters$dr_mat <- shiny::reactiveValuesToList(dr_output)$dr_mat
            t_mat <- t(mat)
            gPlot_mat <- tidyr::tibble("id" = rownames(t_mat))

            if (input$cluster_method == "kmeans"){
                to_create(NULL)

                progress <- shiny::Progress$new()
                on.exit(progress$close())
                progress$set(message = "MoNETA", detail = paste("Doing kmeans"), value = 0)

                kmeans_clus <- kmeans(t_mat, input$k)
                cluster <- gPlot_mat
                cluster$clust <- factor(kmeans_clus$cluster)
                cluster$shape <- NULL

                progress$inc(1, detail = paste("Done!"))

                if (( shiny::isTruthy(annotation()) |  shiny::isTruthy(input$anno_file1)) && ifelse(input$shape_bool == 'YES', TRUE, FALSE)){
                    anno <- if ( shiny::isTruthy(annotation())) anno() else anno1()
                    samples <- colnames(mat)
                    id_col_check <- sapply(colnames(anno), FUN = function(x) sum(samples %in% anno[[x]]))
                    id <- colnames(anno)[which.max(id_col_check)]
                    f_anno <- anno[anno[[id]] %in% samples, c(id, input$select_shape)]
                    colnames(f_anno) <- c('id', 'shape')
                    cluster <- cluster %>% dplyr::inner_join(f_anno, by = 'id')
                }
                clusters$cluster <- cluster

            } else if (input$cluster_method == "dbscan") {
                to_create(NULL)
                progress <- shiny::Progress$new()
                on.exit(progress$close())
                progress$set(message = "MoNETA", detail = paste("Doing DBSCAN"), value = 0)

                db_clust <- fpc::dbscan(t_mat, input$eps)
                cluster <- gPlot_mat
                cluster$clust <- factor(db_clust$cluster)
                cluster$shape <- NULL
                progress$inc(1, detail = paste("Done!"))

                if (( shiny::isTruthy(annotation()) |  shiny::isTruthy(input$anno_file1)) && ifelse(input$shape_bool == 'YES', TRUE, FALSE)){
                    anno <- if ( shiny::isTruthy(annotation())) anno() else anno1()
                    samples <- colnames(mat)
                    id_col_check <- sapply(colnames(anno), FUN = function(x) sum(samples %in% anno[[x]]))
                    id <- colnames(anno)[which.max(id_col_check)]
                    f_anno <- anno[anno[[id]] %in% samples, c(id, input$select_shape)]
                    colnames(f_anno) <- c('id', 'shape')
                    cluster <- cluster %>% dplyr::inner_join(f_anno, by = 'id')
                }
                clusters$cluster <- cluster

            } else if (input$cluster_method == "hclust")  {
                cluster <- NULL
                num_clust <- NULL
                progress <- shiny::Progress$new()
                on.exit(progress$close())

                if (is.null(to_create())){
                    progress$set(message = "MoNETA", detail = paste("Doing hclust"), value = 0)
                    message('hclust')
                    dendro <- t_mat %>%
                        dist() %>%
                        hclust()
                    progress$inc(1, detail = paste("Done!"))
                    message('done!')
                    message('update to_create!')
                    to_create(dendro)
                    message('done!')
                }else {
                    dendro <- to_create()
                    progress$set(message = "MoNETA", detail = paste("Doing Cut tree"), value = 0)
                    message('cut tree')
                    cutree_res <- cutree(dendro, h = input$h)
                    num_clust <- length(unique(cutree_res))
                    cluster <- gPlot_mat
                    cluster$clust <- factor(cutree_res)
                    if (( shiny::isTruthy(annotation()) |  shiny::isTruthy(input$anno_file1)) && ifelse(input$shape_bool == 'YES', TRUE, FALSE)){
                        anno <- if ( shiny::isTruthy(annotation())) anno() else anno1()
                        samples <- colnames(mat)
                        id_col_check <- sapply(colnames(anno), FUN = function(x) sum(samples %in% anno[[x]]))
                        id <- colnames(anno)[which.max(id_col_check)]
                        f_anno <- anno[anno[[id]] %in% samples, c(id, input$select_shape)]
                        colnames(f_anno) <- c('id', 'shape')
                        cluster <- cluster %>% dplyr::inner_join(f_anno, by = 'id')
                    }
                }

                message('done!')
                progress$inc(1, detail = paste("Done!"))

                clusters$cluster <-  cluster
                clusters$dendro <- dendro
                clusters$tree <- as.dendrogram(dendro)
                clusters$num_clust <- num_clust

            } else if (input$cluster_method == "annotation" | input$cluster_method == "annotation1"){
                #if (is.null(anno()) & is.null(anno1())){
                #  return(NULL)
                #}else {
                anno <- if (input$cluster_method == 'annotation') anno() else anno1()
                samples <- colnames(mat)
                id_col_check <- sapply(colnames(anno), FUN = function(x) sum(samples %in% anno[[x]]))
                id <- colnames(anno)[which.max(id_col_check)]
                f_anno <- anno[anno[[id]] %in% samples, c(id, input$custom_anno_selector, input$custom_anno_shape_selector)]
                colnames(f_anno) <- c('id', 'clust', 'shape')
                clusters$cluster <- f_anno
                #}
            }
        }
    })

    shiny::observeEvent(input$cl_btn, {
        shiny::isolate(
            if (input$cluster_method == 'hclust' & !(is.null(shiny::reactiveValuesToList(clusters)$tree))){
                shiny::updateSliderInput(inputId = "h", max = round(base::attr(shiny::reactiveValuesToList(clusters)$tree, "height"), digits = 2)
                )
            }else {
                shiny::updateSliderInput(inputId = "h", max = 100)
            }
        )
    })

    title <-  shiny::eventReactive(input$cl_btn, {
        if (input$cluster_method == "kmeans"){
            return(paste0('Kmeans clustering (k: ', input$k, ')' ))
        } else if (input$cluster_method == "dbscan") {
            return(paste0('DBSCAN clustering (epsilon: ', input$eps, ')' ))
        } else if (input$cluster_method == "hclust") {
            return(paste0('Hierarchical clustering (h: ', input$h, ', #cluster: ', shiny::reactiveValuesToList(clusters)$num_clust, ')' ))
        } else {
            method <- ifelse(input$dr_method == 'UMAP', 'UMAP of',
                             ifelse(input$dr_method == 'PCA', 'PCA of', 'tSNE of'))
            emb <- ifelse(input$emb_opt == 'YES', 'the Embedded RWR matrix', 'the RWR matrix')
            return(paste(method, emb, '(color:', input$custom_anno_selector, ', shape:', input$custom_anno_shape_selector, ')'))
        }
    })


    cl_plot <- shiny::reactive({
        clusters <- shiny::reactiveValuesToList(clusters)
        clu_anno <- clusters$cluster %>% dplyr::arrange(clust)
        if (is.null(clu_anno$shape)){
            cl_plot <- MoNETA::plot_2D_matrix(coord = clusters$dr_mat[1:2,], nodes_anno = clu_anno,
                                              id_name = 'id', id_anno_color = 'clust',
                                              interactive = FALSE, wo_legend = FALSE) +
                ggtitle(title()) +
                xlab(paste0(input$dr_method, '1')) + ylab( paste0(input$dr_method, '2')) +
                theme(text = ggplot2::element_text(family="Times New Roman"),
                      plot.title = ggplot2::element_text(size = 16, hjust = 0.5, face = 'bold'),
                      axis.title.x = ggplot2::element_text(size = 10),
                      axis.title.y = ggplot2::element_text(size = 10))
        } else{
            cl_plot <- MoNETA::plot_2D_matrix(coord = clusters$dr_mat[1:2,], nodes_anno = clu_anno,
                                              id_name = 'id', id_anno_color = 'clust', id_anno_shape = 'shape',
                                              interactive = FALSE, wo_legend = FALSE) +
                ggtitle(title()) +
                xlab(paste0(input$dr_method, '1')) + ylab( paste0(input$dr_method, '2')) +
                theme(text = ggplot2::element_text(family="Times New Roman"),
                      plot.title = ggplot2::element_text(size = 16, hjust = 0.5, face = 'bold'),
                      axis.title.x = ggplot2::element_text(size = 10),
                      axis.title.y = ggplot2::element_text(size = 10))
        }
    })



    ############################################################################################
    #                            OMICS DATA : LOADING DATA (OUTPUT)                      #
    ############################################################################################

    ############## outputs loading info #############

    output$omics_sum_info <- shiny::renderUI({
        input$load_mat_button
        if (length(omics_files()) != 0){
            nsamples <- sapply(omics_files(), ncol)
            nfeatures <- sapply(omics_files(), nrow)
            shinydashboard::box(shiny::HTML(paste0(paste("<b>", names(omics_files()), "</b>"), " loaded ",  paste("<b>", nfeatures, "</b>"), " features and ", paste("<b>", nsamples, "</b>"), " samples <br/>")))
        } else {
            shiny::HTML("<br/>")
        }
    })

    output$anno_info <- shiny::renderUI({
        input$load_anno_button
        if ( shiny::isTruthy(annotation()))
        {
            ncolumns <- ncol(annotation()$annotation)
            nrows <- nrow(annotation()$annotation)
            shinydashboard::box(
                shiny::HTML(paste0(paste("<b>", 'Annotation file', "</b>"), " loaded: ",  paste("<b>", ncolumns, "</b>"), " columns and ", paste("<b>", nrows, "</b>"), " rows <br/>")),
                shiny::HTML(paste0(paste("<b> Column names: </b>"), "<br/>")),
                shiny::HTML(paste('&ensp;', colnames(annotation()$annotation), "<br/>"))
            )
        } else {
            shiny::HTML("<br/>")
        }
    })


    output$anno_info_extra <- shiny::renderUI({
        #input$load_anno_button
        if ( shiny::isTruthy(annotation()) & !is.null(omics_files())){
            samples <- lapply(omics_files(), colnames) %>% unlist() %>% unique()
            anno <- annotation()$annotation
            id_col_check <- sapply(colnames(anno), FUN = function(x) sum(samples %in% anno[[x]]))
            id <- colnames(anno)[which.max(id_col_check)]

            not_annotated_samples <- samples[which(!(samples %in% anno[[id]]))]
            extra_samples_in_anno <- anno[[id]][which(!(anno[[id]] %in% samples))]
            if (length(not_annotated_samples) == 0 & length(extra_samples_in_anno) == 0){
                return(NULL)
            }else {
                message1 <- NULL
                message2 <- NULL
                if (length(extra_samples_in_anno) != 0){
                    message1 <- shiny::HTML(paste('&ensp;', paste("<b>", length(extra_samples_in_anno), "</b>"), " samples are not present in the provided omics data collection", "<br/>"))
                }else if (length(not_annotated_samples) != 0){
                    message2 <- shiny::HTML(paste('&ensp;', paste("<b>", length(not_annotated_samples), "</b>"), " samples are not annotated in the provided annotation file", "<br/>"))
                }
                shinydashboard::box(
                    shiny::HTML(paste("<b>", p('WARNING!', style = "color:red"), "</b>")), message1, message2
                )
            }
        } else {
            shiny::HTML("<br/>")
        }
    })


    ############################################################################################
    #                               OMICS DATA : PRE-PROCESSING (OUTPUT)                       #
    ############################################################################################

    ############## omics matrices pre-processing #############
    processed_mat_list <- list(NULL)
    proc_matrices <- shiny::reactiveValues()

    output$process_omics_mat <- shiny::renderUI({
        if (is.null(omics_files()))
            return(NULL)
        else
            thetabs <- lapply(names(omics_files()), function(x) {

                processed_mat_list[[paste0('process_mat_', x)]] <- shiny::observeEvent(input[[paste0('process_mat_', x)]], {
                    processed_omics_mat <- omics_files()[[x]]
                    if (input[[paste0("remove0col_", x)]] == 'YES')
                        processed_omics_mat <- remove_zeros_cols (processed_omics_mat)
                    if (input[[paste0("norm_", x)]] == 'YES')
                        processed_omics_mat <- normalize_omics(processed_omics_mat)
                    proc_matrices[[paste0('p_', x)]] <- processed_omics_mat
                    print('Pre-processing...done!')
                    return(processed_omics_mat)
                })

                tab <- shiny::tabPanel(x,
                                       shiny::fluidRow(
                                           shiny::column(width = 12,
                                                         shiny::radioButtons(inputId = paste0("remove0col_", x), "Do you want to remove zeros columns?",
                                                                             choices = c('YES', 'NO'), selected = 'YES'
                                                         ),
                                                         shiny::radioButtons(inputId = paste0("norm_", x), "Do you want to normalize the omics matrix by column?",
                                                                             choices = c('YES', 'NO'), selected = 'YES'
                                                         )
                                           )
                                       ), shiny::hr(),
                                       shiny::fluidRow(
                                           shiny::column(width = 12, shiny::actionButton(paste0('process_mat_', x), paste0("Process matrix ", x))),
                                           #shiny::column(width = 6, shiny::uiOutput(paste0('download_omics_net_box_', x)))
                                       )
                )
                return(tab)
            })
        thetabs$width <- 6
        do.call(shinydashboard::tabBox, thetabs)
    })

    ############## Processed matrices info #############

    output$pro_matrices_sum_info <- shiny::renderUI({
        processed_mats <- shiny::reactiveValuesToList(proc_matrices)
        shiny::isolate({
            if (length(processed_mats) != 0){
                nsamples <- sapply(processed_mats, ncol)
                nfeatures <- sapply(processed_mats, nrow)
                shinydashboard::box(shiny::HTML(paste0(paste("<b>", names(processed_mats), "</b>"), " processed ",
                                                       paste("<b>", nfeatures, "</b>"), " features and ",
                                                       paste("<b>", nsamples, "</b>"), " samples <br/>")))
            } else {
                return(NULL)
            }
        })
    })

    ############## download processed omics matrices #############

    output$download_proc_mat_box <- shiny::renderUI({
        if (length(shiny::reactiveValuesToList(proc_matrices)) == 0) {
            return(NULL)
        } else {
            shinydashboard::box(shiny::h4(shiny::span('Download', style = "font-weight: bold")),  solidHeader = TRUE, collapsible = FALSE,
                                checkboxGroupInput('download_proc_mat',
                                                   'Select one or more processed matrices',
                                                   choices = names(shiny::reactiveValuesToList(proc_matrices))
                                ),

                                shiny::tags$hr(),
                                shiny::downloadButton(
                                    outputId = 'download_proc_mat_btn',
                                    label = "Download",
                                    icon = shiny::icon("download")
                                )
            )
        }
    })

    output$download_proc_mat_btn <- shiny::downloadHandler(
        filename = function() {
            paste("MoNETAshiny_processed_mat_", Sys.Date(), ".zip", sep = "")
        },
        content = function(file){

            temp_directory <- file.path(tempdir(), as.integer(Sys.time()))
            dir.create(temp_directory)

            to_download <- shiny::reactiveValuesToList(proc_matrices)[input$download_proc_mat]

            for (obj in names(to_download)) {
                file_name <- glue("{obj}.csv")
                readr::write_csv(as.data.frame(to_download[[obj]]), file.path(temp_directory, file_name))
            }

            zip::zip(
                zipfile = file,
                files = dir(temp_directory),
                root = temp_directory
            )
        }
    )




    ############## not mandatory intersection #############

    output$intersection <- shiny::renderUI({
        loaded_mats <- omics_files()
        processed_mats <- shiny::reactiveValuesToList(proc_matrices)
        if (length(loaded_mats) > 1 & length(processed_mats) == length(loaded_mats))
            shinydashboard::box(
                shiny::radioButtons(inputId = 'intersect_opt',
                                    label = 'Do you want to intersect omics matrices by samples?',
                                    choices = c('YES', 'NO'), selected = 'YES'
                ),
                shiny::hr(),
                shiny::actionButton(inputId = 'intersect_btn', label = 'Go')
            )
    })

    ############## intersection info #############

    output$intersection_info <- shiny::renderUI({
        int_mats <- shiny::reactiveValuesToList(intersected_omics_mat)$matrices
        shiny::isolate({
            if (length(int_mats) != 0){
                nsamples <- sapply(int_mats, ncol)
                nfeatures <- sapply(int_mats, nrow)
                shinydashboard::box(shiny::HTML(paste0(paste("<b>", names(int_mats), "</b>"), " processed ",
                                                       paste("<b>", nfeatures, "</b>"), " features and ",
                                                       paste("<b>", nsamples, "</b>"), " samples <br/>")))
            } else {
                return(NULL)
            }
        })
    })


    ############################################################################################
    #                        NETWORK : NETWORK INFERENCE (OUTPUT)                      #
    ############################################################################################

    ############## output parameters for OMICS network inference  #############

    omics_network_list <- list(NULL)
    gene_networks <- shiny::reactiveValues()

    output$omics_net_arguments <- shiny::renderUI({
        if (length(shiny::reactiveValuesToList(intersected_omics_mat)$matrices) == 0){
            if (length(omics_files()) == 1){
                req(proc_matrices)
                mats <- shiny::reactiveValuesToList(proc_matrices)
            }else{
                return(NULL)
            }
        }else {
            req(intersected_omics_mat)
            mats <- shiny::reactiveValuesToList(intersected_omics_mat)$matrices
        }

        thetabs <- lapply(names(mats), function(x) {

            omics_network_list[[paste0('generate_omics_net_', x)]] <- shiny::observeEvent(input[[paste0('generate_omics_net_', x)]], {
                progress <- shiny::Progress$new()
                on.exit(progress$close())

                progress$set(message = "MoNETA", detail = paste("Inferring the network..."), value = 0)
                net <- MoNETA::k_star_net(matrix = mats[[x]],
                                          sparsity = input[[paste0('sparsity_', x)]],
                                          distFun = input[[paste0('distFUN_', x)]],
                                          cores = input[[paste0('cores_', x)]],
                                          knn = input[[paste0('knn_', x)]],
                                          k_star = ifelse(input[[paste0('kstar_', x)]] == 'YES', TRUE, FALSE),
                                          MAX_ASSOC = Inf)

                progress$inc(1, detail = paste("Network constructed!"))
                gene_networks[[paste0(x, '_net')]] <- net
                print('Omics Network Inference...done!')

                output[[paste0(x, 'plot_net')]] <- visNetwork::renderVisNetwork({
                    ids <- data.frame(id = base::union(net$source, net$dest))
                    net_plot <- plot_net(edgeList = net, nodes_anno = ids, id_name = 'id', id_anno_color = 'id')
                    return(net_plot)
                })

                return(net)
            })

            omics_network_list[[paste0('update_omics_net_plot_', x)]] <-
                shiny::observeEvent(input[[paste0('update_omics_net_plot_', x)]], {

                    output[[paste0(x, 'plot_net')]] <- visNetwork::renderVisNetwork({
                        shiny::isolate({
                            net <- gene_networks[[paste0(x, '_net')]]
                            anno <- anno()
                            color <- if (input[[paste0(x, '_net_color')]] != '-') input[[paste0(x, '_net_color')]] else NA
                            shape <- if (input[[paste0(x, '_net_shape')]] != '-') input[[paste0(x, '_net_shape')]] else NA

                            samples <- base::union(net$source, net$dest)
                            id_col_check <- sapply(colnames(anno), FUN = function(x) sum(samples %in% anno[[x]]))
                            id <- colnames(anno)[which.max(id_col_check)]
                            f_anno <- anno[anno[[id]] %in% samples, ]
                            final_anno <- f_anno %>% tidyr::as_tibble(.) %>% dplyr::select(id, everything())

                            net_plot <- plot_net(edgeList = net, nodes_anno = final_anno,
                                                 id_name = id, id_anno_color = color, id_anno_shape = shape,
                                                 interactive = TRUE, wo_legend = FALSE)
                            return(net_plot)
                        })
                    })
                })

            output[[paste0(x, '_custom_net')]] <- shiny::renderUI({
                req(input[[paste0('generate_omics_net_', x)]])
                shiny::isolate({
                    if ( shiny::isTruthy(annotation())){
                        output <- tagList()
                        output[[1]] <-
                            shiny::selectInput(inputId = paste0(x, '_net_color'),
                                               label = 'Select the color of the nodes',
                                               #choices = c('-', colnames(annotation()$annotation)), selected = '-')
                                               choices = colnames(annotation()$annotation), selected = NULL)
                        output[[2]] <-
                            shiny::selectInput(inputId = paste0(x, '_net_shape'),
                                               label = 'Select the shape of the nodes',
                                               choices = c('-', colnames(annotation()$annotation)), selected = '-')
                        output
                    }else{
                        return(NULL)
                    }
                })
            })

            output[[paste0('update_omics_net_plot_', x)]] <- shiny::renderUI({
                req(input[[paste0('generate_omics_net_', x)]])
                shiny::isolate({
                    if ( shiny::isTruthy(annotation())){
                        shiny::actionButton(inputId = paste0('update_omics_net_plot_', x), label = 'Update plot')
                    }else{
                        return(NULL)
                    }
                })
            })

            tab <- shiny::tabPanel(x,
                                   shiny::fluidRow(
                                       shiny::column(width = 4,
                                                     shiny::selectInput(
                                                         inputId = paste0('distFUN_', x),
                                                         label = "Select a distance measure",
                                                         choices =  c("Euclidean", "Manhattan", "Cosine"),
                                                         selected = "Euclidean"
                                                     ),
                                                     shiny::numericInput(
                                                         inputId = paste0('sparsity_', x),
                                                         label = "Select the sparsity",
                                                         min = 0, max = 1, value = 0.7, step = 0.1
                                                     ),
                                                     shiny::numericInput(
                                                         inputId = paste0('knn_', x),
                                                         label = "Select the number of neighbors",
                                                         min = 1, max = 100, value = 25
                                                     ),
                                                     shiny::numericInput(
                                                         inputId = paste0('cores_', x),
                                                         label = "Select the number of cores",
                                                         min = 1, max = 200, value = 1
                                                     ),
                                                     shiny::radioButtons(inputId = paste0("kstar_", x), "Do you want to perform the kstar?",
                                                                         choices = c('YES', 'NO'), selected = 'YES'
                                                     )
                                       ),
                                       shiny::column(width = 8,
                                                     visNetwork::visNetworkOutput(paste0(x, 'plot_net')),
                                                     shiny::uiOutput(paste0(x, '_custom_net'))
                                       )
                                   ), shiny::hr(),
                                   shiny::fluidRow(
                                       shiny::column(width = 6, shiny::actionButton(paste0('generate_omics_net_', x), paste0("Generate network ", x))),
                                       shiny::column(width = 6, shiny::uiOutput(paste0('update_omics_net_plot_', x)))
                                   )
            )
            return(tab)
        })
        thetabs$width <- 8
        do.call(shinydashboard::tabBox, thetabs)
    })

    ############## download networks  #############

    output$download_net_box <- shiny::renderUI({
        if (length(shiny::reactiveValuesToList(gene_networks)) == 0) {
            return(NULL)
        } else {
            shinydashboard::box(shiny::h4(shiny::span('Download', style = "font-weight: bold")), solidHeader = TRUE, collapsible = FALSE, width = 4,
                                checkboxGroupInput('download_net',
                                                   'Select one or more inferred networks',
                                                   choices = names(shiny::reactiveValuesToList(gene_networks))
                                ),
                                shiny::tags$hr(),
                                shiny::downloadButton(
                                    outputId = 'download_net_btn',
                                    label = "Download",
                                    icon = shiny::icon("download")
                                )
            )
        }
    })

    output$download_net_btn <- shiny::downloadHandler(
        filename = function() {
            paste("MoNETAshiny_inferred_net_", Sys.Date(), ".zip", sep = "")
        },
        content = function(file){

            temp_directory <- file.path(tempdir(), as.integer(Sys.time()))
            dir.create(temp_directory)

            to_download <- shiny::reactiveValuesToList(gene_networks)[input$download_net]

            for (obj in names(to_download)) {
                file_name <- glue("{obj}.csv")
                readr::write_csv(as.data.frame(to_download[[obj]]), file.path(temp_directory, file_name))
            }

            zip::zip(
                zipfile = file,
                files = dir(temp_directory),
                root = temp_directory
            )
        }
    )


    ############################################################################################
    #                               NETWORK : NETWORK FILTERING (OUTPUT)                       #
    ############################################################################################

    ############## Filtering networks box #############

    filteredOmicsNetworks_list <- list(NULL)
    filteredOmicsNetworks <- shiny::reactiveValues()

    output$plot_net_box <- shiny::renderUI({
        if (length(shiny::reactiveValuesToList(gene_networks)) == 0){
            return(NULL)
        }else
            nets <- shiny::reactiveValuesToList(gene_networks)

        thetabs <- lapply(names(nets), function(x) {

            output[[paste0(x, '_info')]] <-
                shiny::renderUI({
                    input[[ paste0('button_', x)]]
                    shiny::isolate({
                        net <- nets[[x]]
                        fnet <- net[net$weight <= input[[paste0('omics_range', x)]],]
                        summary <- netSummary(tidyr::tibble(fnet))
                        res <- c()
                        for (i in 2:length(summary))
                            res[i-1] <- (paste0(paste('<b>', names(summary)[i], '</b>'), ': ', toString(summary[[i]]), ' <br/>'))
                        shiny::HTML(res)
                    })
                })

            filteredOmicsNetworks_list[[paste0('button_', x)]] <- shiny::observeEvent(input[[paste0('button_', x)]], {
                net <- nets[[x]]
                fnet <- net[net$weight <= input[[paste0('omics_range', x)]],]

                output[[paste0(x, "_weightDist")]] <- plotly::renderPlotly({
                    shiny::isolate({
                        plotly::ggplotly(ggplot2::ggplot(fnet, ggplot2::aes(x= weight)) + ggplot2::geom_histogram() +
                                     ggplot2::ggtitle('Weight Distribution') +
                                     ggplot2::xlab("Weight") + ggplot2::ylab("Counts") +
                                     ggplot2::theme(text = ggplot2::element_text(family="LM Roman 10"),
                                                    plot.title = ggplot2::element_text(size = 16, hjust = 0.5, face = 'bold'),
                                                    axis.title.x = ggplot2::element_text(size = 10),
                                                    axis.title.y = ggplot2::element_text(size = 10))
                        )
                    })
                })

                output[[paste0(x, "_nodeDegreeDist")]] <- plotly::renderPlotly({
                    shiny::isolate({
                        fnet1 <- fnet
                        fnet1$weight <- 1
                        graph <- igraph::graph_from_data_frame(fnet1, directed = FALSE)
                        degree_table <- data.frame(Degree = igraph::degree(graph))
                        p <- plotly::ggplotly(ggplot2::ggplot(degree_table, ggplot2::aes(x = Degree)) + ggplot2::geom_histogram() +
                                          ggplot2::ggtitle('Nodes Degree Distribution') +
                                          ggplot2::ylab("Number of nodes") + ggplot2::xlab("Node degree") +
                                          ggplot2::theme(text = ggplot2::element_text(family="LM Roman 10"),
                                                         plot.title = ggplot2::element_text(size = 16, hjust = 0.5, face = 'bold'),
                                                         axis.title.x = ggplot2::element_text(size = 10),
                                                         axis.title.y = ggplot2::element_text(size = 10))
                        )
                        return(p)
                    })
                })

                output[[paste0(x, '_plot_fnet')]] <- visNetwork::renderVisNetwork({
                    shiny::isolate({
                        ids <- data.frame(id = base::union(fnet$source, fnet$dest))
                        net_plot <- plot_net(edgeList = fnet, nodes_anno = ids, id_name = 'id', id_anno_color = 'id')
                        return(net_plot)
                    })
                })

                output[[paste0(x, '_custom_fnet')]] <- shiny::renderUI({
                    #shiny::isolate({
                    if ( shiny::isTruthy(annotation())){
                        output <- tagList()
                        output[[1]] <-
                            shiny::selectInput(inputId = paste0(x, '_fnet_color'),
                                               label = 'Select the color of the nodes',
                                               #choices = c('-', colnames(annotation()$annotation)), selected = '-')
                                               choices =colnames(annotation()$annotation), selected = NULL)
                        output[[2]] <-
                            shiny::selectInput(inputId = paste0(x, '_fnet_shape'),
                                               label = 'Select the shape of the nodes',
                                               choices = c('-', colnames(annotation()$annotation)), selected = '-')
                        output
                    }else{
                        return(NULL)
                    }
                    #})
                })
            })


            filteredOmicsNetworks_list[[paste0("update_", x)]] <- shiny::observeEvent(input[[paste0("update_", x)]], {
                #shiny::isolate({
                filt_net <- net[net$weight <= input[[paste0('omics_range', x)]],]
                filteredOmicsNetworks[[paste0('f_', x)]] <- filt_net
                print('Filering network...done!')

                return(filt_net)
                #})
            })

            filteredOmicsNetworks_list[[paste0('update_f_net_plot_', x)]] <-
                shiny::observeEvent(input[[paste0('update_f_net_plot_', x)]], {

                    output[[paste0(x, '_plot_fnet')]] <- visNetwork::renderVisNetwork({
                        shiny::isolate({
                            net <- nets[[x]]
                            filt_net <- net[net$weight <= input[[paste0('omics_range', x)]],]
                            anno <- anno()
                            color <- if (input[[paste0(x, '_fnet_color')]] != '-') input[[paste0(x, '_fnet_color')]] else NA
                            shape <- if (input[[paste0(x, '_fnet_shape')]] != '-') input[[paste0(x, '_fnet_shape')]] else NA

                            samples <- base::union(filt_net$source, filt_net$dest)
                            id_col_check <- sapply(colnames(anno), FUN = function(x) sum(samples %in% anno[[x]]))
                            id <- colnames(anno)[which.max(id_col_check)]
                            f_anno <- anno[anno[[id]] %in% samples, ]
                            final_anno <- f_anno %>% tidyr::as_tibble(.) %>%  dplyr::select(id, everything())

                            net_plot <- plot_net(edgeList = filt_net, nodes_anno = final_anno,
                                                 id_name = id, id_anno_color = color, id_anno_shape = shape,
                                                 interactive = TRUE, wo_legend = FALSE)
                            return(net_plot)
                        })
                    })
                })


            ### Update plot button
            output[[paste0('update_f_net_plot_', x)]] <- shiny::renderUI({
                #shiny::isolate({
                if ( shiny::isTruthy(annotation()) &  shiny::isTruthy(input[[paste0('button_', x)]])){
                    shiny::actionButton(inputId = paste0('update_f_net_plot_', x), label = 'Update plot')
                }else{
                    return(NULL)
                }
                #})
            })

            net <- nets[[x]]
            #output[[paste0('table_net', x)]] <- renderTable(net)
            min <-  round(min(net$weight), digits = 2)
            max <-  round(max(net$weight), digits = 2)

            tab <- shiny::tabPanel(x,
                                   shiny::fluidRow(
                                       shiny::column(width = 4,
                                                     shiny::sliderInput(inputId = paste0("omics_range", x), "Select a threshold:",
                                                                        min = min, max = max, value = max, step = 0.01),
                                                     shiny::hr(),
                                                     shiny::actionButton(inputId = paste0('button_', x), 'Show')
                                       ),
                                       shiny::column(width = 8,
                                                     shiny::uiOutput(paste0(x, '_info'))
                                       )
                                   ),
                                   shiny::hr(),
                                   shiny::fluidRow(
                                       shiny::column(width = 6,
                                                     plotly::plotlyOutput(paste0(x, "_weightDist"))
                                       ),
                                       shiny::column(width = 6,
                                                     plotly::plotlyOutput(paste0(x, "_nodeDegreeDist"))
                                       )
                                   ),
                                   shiny::fluidRow(
                                       shiny::column(width = 8,  visNetwork::visNetworkOutput(paste0(x, '_plot_fnet'))),
                                       shiny::column(width = 4, shiny::uiOutput(paste0(x, '_custom_fnet')))
                                   ),
                                   shiny::hr(),
                                   shiny::fluidRow(
                                       shiny::column(width = 6, shiny::actionButton(paste0("update_", x), paste('UPDATE ', x))),
                                       shiny::column(width = 6, shiny::uiOutput(paste0("update_f_net_plot_", x)))
                                   )
            )
            return(tab)
        })

        thetabs$width <- 8
        do.call(shinydashboard::tabBox, thetabs)
    })

    ############## Filtered networks info #############

    output$f_net_sum_info <- shiny::renderUI({
        f_nets <- shiny::reactiveValuesToList(filteredOmicsNetworks)
        shiny::isolate({
            if (length(f_nets) != 0){
                n_nodes <- sapply(f_nets, function(x){length(base::union(x[[1]], x[[2]]))})
                n_edges <- sapply(f_nets, nrow)
                shinydashboard::box(width = 4,
                                    shiny::HTML(paste0(paste("<b>", names(f_nets), "</b>"), " filtered: ",
                                                       paste("<b>", n_nodes, "</b>"), " nodes and ",
                                                       paste("<b>", n_edges, "</b>"), " edges <br/>")))
            } else {
                return(NULL)
            }
        })
    })

    ############## download filtered networks  #############

    output$download_fnet_box <- shiny::renderUI({
        if (length(shiny::reactiveValuesToList(filteredOmicsNetworks)) == 0) {
            return(NULL)
        } else {
            shinydashboard::box(shiny::h4(shiny::span('Download', style = "font-weight: bold")), solidHeader = TRUE, collapsible = FALSE, width = 4,
                                checkboxGroupInput('download_fnet',
                                                   'Select one or more filtered networks',
                                                   choices = names(shiny::reactiveValuesToList(filteredOmicsNetworks))
                                ),
                                shiny::tags$hr(),
                                shiny::downloadButton(
                                    outputId = 'download_fnet_btn',
                                    label = "Download",
                                    icon = shiny::icon("download")
                                )
            )
        }
    })

    output$download_fnet_btn <- shiny::downloadHandler(
        filename = function() {
            paste("MoNETAshiny_filtered_net_", Sys.Date(), ".zip", sep = "")
        },
        content = function(file){

            temp_directory <- file.path(tempdir(), as.integer(Sys.time()))
            dir.create(temp_directory)

            to_download <- shiny::reactiveValuesToList(filteredOmicsNetworks)[input$download_fnet]

            for (obj in names(to_download)) {
                file_name <- glue("{obj}.csv")
                readr::write_csv(as.data.frame(to_download[[obj]]), file.path(temp_directory, file_name))
            }

            zip::zip(
                zipfile = file,
                files = dir(temp_directory),
                root = temp_directory
            )
        }
    )


    ############################################################################################
    #                             RWR: LOADING NETWORKS (OUTPUT)                          #
    ############################################################################################

    ############## outputs loading info #############

    output$omics_net_sum_info <- shiny::renderUI({
        if (length(loaded_omics_net_list()) != 0){
            n_nodes <- sapply(loaded_omics_net_list(), function(x){length(unique(c(x[[1]], x[[2]])))})
            n_edges <- sapply(loaded_omics_net_list(), nrow)
            shinydashboard::box(shiny::HTML(paste0(paste("<b>", names(loaded_omics_net_list()), "</b>"), " loaded: ",  paste("<b>", n_nodes, "</b>"), " nodes and ", paste("<b>", n_edges, "</b>"), " edges <br/>")))
        } else {
            shiny::HTML("<br/>")
        }
    })

    ############## annotation1  info #############

    output$anno_info1 <- shiny::renderUI({
        input$load_anno_button1
        if ( shiny::isTruthy(annotation1()))
        {
            ncolumns <- ncol(annotation1()$annotation1)
            nrows <- nrow(annotation1()$annotation1)
            shinydashboard::box(
                shiny::HTML(paste0(paste("<b>", 'Annotation file', "</b>"), " loaded: ",  paste("<b>", ncolumns, "</b>"), " columns and ", paste("<b>", nrows, "</b>"), " rows <br/>")),
                shiny::HTML(paste0(paste("<b> Colnames </b>"), "<br/>")),
                shiny::HTML(paste('&ensp;', colnames(annotation1()$annotation1), "<br/>"))
            )
        } else {
            shiny::HTML("<br/>")
        }
    })

    output$anno_info_extra1 <- shiny::renderUI({
        if ( shiny::isTruthy(annotation1()) && length(omics_files()) != 0){
            samples <- lapply(omics_files(), colnames) %>% unlist() %>% unique()
            anno <- annotation1()$annotation
            id_col_check <- sapply(colnames(anno), FUN = function(x) sum(samples %in% anno[[x]]))
            id <- colnames(anno)[which.max(id_col_check)]
            not_annotated_samples <- samples[which(!(samples %in% anno[[id]]))]
            extra_samples_in_anno <- anno[[id]][which(!(anno[[id]] %in% samples))]
            if (length(not_annotated_samples) == 0 & length(extra_samples_in_anno) == 0){
                return(NULL)
            }else {
                message1 <- NULL
                message2 <- NULL
                if (length(extra_samples_in_anno) != 0){
                    message1 <- shiny::HTML(paste('&ensp;', paste("<b>", length(extra_samples_in_anno), "</b>"), " samples are not present in the omics data collection", "<br/>"))
                }else if (length(not_annotated_samples) != 0){
                    message2 <- shiny::HTML(paste('&ensp;', paste("<b>", length(not_annotated_samples), "</b>"), " samples are not annotated in the provided annotation file", "<br/>"))
                }
                shinydashboard::box(
                    shiny::HTML(paste("<b>", p('WARNING!', style = "color:red"), "</b>")), message1, message2
                )
            }
        } else {
            shiny::HTML("<br/>")
        }
    })

    output$anno_info_extra2 <- shiny::renderUI({
        if ( shiny::isTruthy(annotation1()) && length(loaded_omics_net_list()) != 0){
            samples <- lapply(loaded_omics_net_list(), function(x) {base::union(x[,1], x[,2])}) %>% unlist() %>% unique()

            anno <- annotation1()$annotation
            id_col_check <- sapply(colnames(anno), FUN = function(x) sum(samples %in% anno[[x]]))
            id <- colnames(anno)[which.max(id_col_check)]

            not_annotated_samples <- samples[which(!(samples %in% anno[[id]]))]
            extra_samples_in_anno <- anno[[id]][which(!(anno[[id]] %in% samples))]
            if (length(not_annotated_samples) == 0 & length(extra_samples_in_anno) == 0){
                return(NULL)
            }else {
                message1 <- NULL
                message2 <- NULL
                if (length(extra_samples_in_anno) != 0){
                    message1 <- shiny::HTML(paste('&ensp;', paste("<b>", length(extra_samples_in_anno), "</b>"), " samples are not present in the omics data collection", "<br/>"))
                }else if (length(not_annotated_samples) != 0){
                    message2 <- shiny::HTML(paste('&ensp;', paste("<b>", length(not_annotated_samples), "</b>"), " samples are not annotated in the provided annotation file", "<br/>"))
                }
                shinydashboard::box(
                    shiny::HTML(paste("<b>", p('WARNING!', style = "color:red"), "</b>")), message1, message2
                )
            }
        } else {
            shiny::HTML("<br/>")
        }
    })


    ############################################################################################
    #                               RWR: MULTIPLEX NETWORK (OUTPUT)                            #
    ############################################################################################

    output$multiplex_func <- shiny::renderUI({
        #if (is.null(g_net_list()) & is.null(l_net_list())){
        #  return(NULL)
        #}else {
        boxes <- list(
            shiny::fluidRow(
                shinydashboard::box(
                    shiny::radioButtons(inputId = 'weightMultiplex', label = 'Do you want to create a weighted multiplex network?',
                                        choices = c("YES", "NO"), selected = "NO"),
                    shiny::hr(),
                    shiny::actionButton(inputId = 'gen_multiplex_btn', label = 'GO')),
                shiny::htmlOutput(outputId = "multi_net_sum_info")
            ),
            shiny::fluidRow(
                shiny::conditionalPanel(
                    condition = 'input.gen_multiplex_btn >= 1 && input.weightMultiplex == "YES"',
                    shinydashboard::box(
                        shiny::radioButtons(inputId = 'pruneMultiplex', label = 'Do you want to prune the multiplex network?',
                                            choices = c('YES', 'NO'), selected = 'NO'),
                        shiny::conditionalPanel(
                            condition = 'input.pruneMultiplex == "YES"',
                            shiny::numericInput(inputId = 'pruneMultiplex_th', label = 'Select a threshold to prune the multiplex network',
                                                min = 0, max = 100, value = 50)
                        ),
                        shiny::hr(),
                        shiny::actionButton(inputId = 'prune_multiplex_btn', label = 'Go')
                    ),
                    shiny::htmlOutput(outputId = "pruned_multi_net_sum_info")
                )
            )
        )
        return(boxes)
        #}
    })

    ############## output multiplex net info #############

    output$multi_net_sum_info <- shiny::renderUI({
        input$gen_multiplex_btn
        shiny::isolate({
            req(multiplex_network)
            multiplex <- shiny::reactiveValuesToList(multiplex_network)$multiplex
            if (!is.null(multiplex)){
                n_layers <- multiplex %>%  dplyr::select(EdgeType) %>% unique(.) %>% nrow(.)
                n_nodes <- base::union(multiplex$source, multiplex$target) %>% length(.)
                n_edges <- multiplex %>% nrow(.)
                weighted <- ifelse(input$weightMultiplex == 'YES', 'Weighted', 'Unweighted')
                shinydashboard::box(shiny::HTML(paste0(paste("<b>", weighted ,"Multiplex network</b>"), " created: <br/>")),
                                    shiny::HTML(paste('&ensp;', "<b>", n_layers, "</b>", " layers, <br/>")),
                                    shiny::HTML(paste('&ensp;',"<b>", n_nodes, "</b>", " nodes <br/>")),
                                    shiny::HTML(paste('&ensp;',"<b>", n_edges, "</b>", "edges <br/>"))
                )
            } else {
                shiny::HTML("<br/>")
            }
        })
    })

    ############## output pruned multiplex net info #############

    output$pruned_multi_net_sum_info <- shiny::renderUI({
        input$prune_multiplex_btn
        shiny::isolate({
            pruned_multiplex <- shiny::reactiveValuesToList(multiplex_network)$pruned_multiplex
            if (!is.null(pruned_multiplex)){
                n_layers <- pruned_multiplex %>%  dplyr::select(EdgeType) %>% unique(.) %>% nrow(.)
                n_nodes <- base::union(pruned_multiplex$source, pruned_multiplex$target) %>% length(.)
                n_edges <- pruned_multiplex %>% nrow(.)
                shinydashboard::box(shiny::HTML(paste0(paste("<b>Weighted Multiplex network</b>"), " pruned: <br/>")),
                                    shiny::HTML(paste('&ensp;', "<b>", n_layers, "</b>", " layers, <br/>")),
                                    shiny::HTML(paste('&ensp;',"<b>", n_nodes, "</b>", " nodes <br/>")),
                                    shiny::HTML(paste('&ensp;',"<b>", n_edges, "</b>", "edges <br/>"))
                )
            } else {
                shiny::HTML("<br/>")
            }
        })
    })


    ############################################################################################
    #                                    RWR: PARAMETERS (OUTPUT)                              #
    ############################################################################################

    #### Taus
    output$tauBIO <- shiny::renderUI({
        numInput <- if (!is.null(g_net_list())) length(g_net_list()) else length(l_net_list())
        lapply(1:numInput, function(i) {
            shiny::numericInput(
                inputId = paste0('omics_tau', i),
                label = shiny::HTML("&tau;", i),
                min = 0,
                max = 1,
                value = as.double(1/numInput),
                step = 0.1)
        })
    })

    omics_tau_list <- shiny::reactive({
        numInput <- if (!is.null(g_net_list())) length(g_net_list()) else length(l_net_list())
        lapply(1:numInput, function(i) {
            input[[paste0("omics_tau", i)]]
        })
    })

    #### Layer Transition matrix
    observe({
        net_list <- if (!is.null(g_net_list())) g_net_list() else l_net_list()
        m <- shiny::reactive({
            L1 <- length(net_list)
            mat <- matrix(1, ncol = L1 , nrow = L1)
            mat <- mat/(L1-1)
            diag(mat) <- 0
            colnames(mat) <- rownames(mat) <- names(net_list)
            return(mat)
        })
        output$omics_trans_mat <- shiny::renderUI({
            div(
                shinyMatrix::matrixInput(inputId = "editable_trans_mat", value = m(), class = "numeric")
            )
        })
    })


    ############################################################################################
    #                               DIMENSIONALITY REDUCTION (OUTPUT)                          #
    ############################################################################################

    output$dr_opt <- shiny::renderUI({
        shinydashboard::box(width = 4,
                            shiny::radioButtons(inputId = 'emb_opt', label = 'Do you want to denoise the RWR-mat via embedding?',
                                                choices = c('YES', 'NO'), selected = 'NO'),
                            shiny::conditionalPanel(
                                condition = 'input.emb_opt== "YES"',
                                shiny::numericInput(inputId = 'dim_emb', label = 'Select the latent dimensions of the embedding',
                                                    min = 1, max = 100, value = 64),
                                shiny::numericInput(inputId = 'cores_emb', label = 'Select the number of CPU for the embedding function',
                                                    min = 1, max = 100, value = 1)
                            ),
                            shiny::hr(),
                            shiny::selectInput(inputId = 'dr_method', label = 'Select a dimensionality reduction method',
                                               choices = c('PCA', 'UMAP', 'tSNE'), selected = 'UMAP'),
                            shiny::conditionalPanel(
                                condition = 'input.dr_method == "PCA"',
                                shiny::numericInput(inputId = 'emb_size', 'Select the SIZE of the output embedding', value = 64),
                            ),
                            shiny::conditionalPanel(
                                condition = 'input.dr_method == "tSNE"',
                                shiny::numericInput(inputId = 'emb_size_tsne', 'Select the SIZE of the output embedding', min = 1, max = 3, value = 2),
                                shiny::numericInput(inputId = 'max_iter', label = 'Select the number of iterations', value = 20000),
                                shiny::numericInput(inputId = 'perplexity', label = 'Select the PERPLEXITY parameter', min = 1, max = 100, value = 5)
                            ),
                            shiny::conditionalPanel(
                                condition = 'input.dr_method != "PCA"',
                                shiny::numericInput(inputId = 'threads', label = 'Select the number of CPU for umap function', value = 1)
                            ),
                            shiny::hr(),
                            shiny::actionButton("dr_btn", "RUN")
        )
    })

    output$dr_plot_box <- shiny::renderUI({
        input$dr_btn
        shiny::isolate({
            if (length(shiny::reactiveValuesToList(dr_output)) == 0){
                return(NULL)
            }else{
                method <- ifelse(input$dr_method == 'UMAP', 'UMAP of',
                                 ifelse(input$dr_method == 'PCA', 'PCA of', 'tSNE of'))
                emb <- ifelse(input$emb_opt == 'YES', 'the Embedded RWR matrix', 'the RWR matrix')
                dr_mat <- shiny::reactiveValuesToList(dr_output)$dr_mat
                dr_plot <- MoNETA::plot_2D_matrix(coord = dr_mat[1:2,], nodes_anno = data.frame(id = colnames(dr_mat)),
                                                  id_name = 'id',
                                                  interactive = TRUE, wo_legend = FALSE) %>%
                    plotly::layout(font = list(family = "LMRoman10", color = "black"),
                                   title = list(text = paste(method, emb), font = list(family = "LMRoman10", color = "black", size = 20)),
                                   xaxis = list(title = list(text = paste0(input$dr_method, '1'))),
                                   yaxis = list(title = list(text =  paste0(input$dr_method, '2')))
                    )%>%
                    plotly::config(
                        toImageButtonOptions = list(filename = paste('MoNETAshiny', gsub(x = method, replacement = '_', pattern = ' '),
                                                                     gsub(x = emb, replacement = '_', pattern = ' '), 'plot', sep='_'),
                                                    format = "png", width = 1800, height = 800)
                    )
                return(shinydashboard::box(width = 8, dr_plot, shiny::hr(), shiny::downloadButton(outputId = 'download_dr_plot', label = 'Download plot')))
            }
        })
    })

    ############## download filtered networks  #############

    output$download_dr_box <- shiny::renderUI({
        if (length(shiny::reactiveValuesToList(dr_output)) == 0) {
            return(NULL)
        } else {
            shinydashboard::box(shiny::h4(shiny::span('Download', style = "font-weight: bold")), solidHeader = TRUE, collapsible = FALSE, width = 4,
                                checkboxGroupInput('download_dr',
                                                   'Select one or more low-dimensional RWR matrices',
                                                   choices = names(shiny::reactiveValuesToList(dr_output))
                                ),
                                shiny::tags$hr(),
                                shiny::downloadButton(
                                    outputId = 'download_dr_btn',
                                    label = "Download",
                                    icon = shiny::icon("download")
                                )
            )
        }
    })

    output$download_dr_btn <- shiny::downloadHandler(
        filename = function() {
            paste("MoNETAshiny_dimensionalityReduction_", Sys.Date(), ".zip", sep = "")
        },
        content = function(file){

            temp_directory <- file.path(tempdir(), as.integer(Sys.time()))
            dir.create(temp_directory)

            to_download <- shiny::reactiveValuesToList(dr_output)[input$download_dr]

            for (obj in names(to_download)) {
                file_name <- if (obj == 'dr_mat') glue("{input$dr_method}_{obj}.csv") else glue("{obj}.csv")
                readr::write_csv(as.data.frame(to_download[[obj]]), file.path(temp_directory, file_name))
            }

            zip::zip(
                zipfile = file,
                files = dir(temp_directory),
                root = temp_directory
            )
        }
    )

    ############## download plot cluster in 2D #############

    output$download_dr_plot <- shiny::downloadHandler(
        filename = function() { paste0('MoNETAshiny_', input$dr_method, '_plot.png', sep='') },
        content = function(file) {
            ggplot2::ggsave(filename = file, plot = dr_plot(), device = "png", width = 10, height = 5, limitsize = FALSE)
        }
    )


    ############################################################################################
    #                                     CLUSTERING  (OUTPUT)                                 #
    ############################################################################################

    ############## customized Cluster selector #############

    output$custom_anno_selector <- shiny::renderUI({
        if (input$cluster_method == 'annotation'){
            anno_tab <- anno()
            output <- tagList()
            output[[1]] <-
                shiny::selectInput(inputId = 'custom_anno_selector', label = 'Select a column to color data points',
                                   choices = colnames(anno_tab), selected = NULL)
            output[[2]] <-
                shiny::selectInput(inputId = 'custom_anno_shape_selector', label = 'Select a column to shape data points',
                                   choices = colnames(anno_tab), selected = NULL)
            output
        } else if (input$cluster_method == 'annotation1'){
            anno_tab <- anno1()
            output <- tagList()
            output[[1]] <-
                shiny::selectInput(inputId = 'custom_anno_selector', label = 'Select a column to color data points',
                                   choices = colnames(anno_tab), selected = NULL)
            output[[2]] <-
                shiny::selectInput(inputId = 'custom_anno_shape_selector', label = 'Select a column to shape data points',
                                   choices = colnames(anno_tab), selected = NULL)
            output
        }
    })

    ############## output cluster box #############

    output$cl_opt <- shiny::renderUI({
        anno_name <- if (! shiny::isTruthy(annotation())) NULL else 'annotation'
        anno1_name <-  if (! shiny::isTruthy(input$anno_file1)) NULL else 'annotation1'
        anno_list <- list(anno_name, anno1_name)

        output$shape_opt <- shiny::renderUI({
            if ( ( shiny::isTruthy(annotation()) |  shiny::isTruthy(input$anno_file1) ) & input$cluster_method %in% c('kmeans', 'dbscan', 'hclust') ){
                anno_table <- if ( shiny::isTruthy(annotation())) anno() else anno1()
                output <-  tagList()
                output[[1]] <-
                    shiny::radioButtons(inputId = 'shape_bool', label = 'Do you want to use particular shape for points?', c('YES', 'NO'), selected = 'NO')
                output[[2]] <-
                    shiny::conditionalPanel(
                        condition = 'input.shape_bool == "YES"',
                        shiny::selectInput(inputId = 'select_shape', label = 'Select a column to shape data points',  choices = colnames(anno_table))
                    )
                output
            }else if (input$cluster_method %in% c('kmeans', 'dbscan', 'hclust') & (! shiny::isTruthy(annotation()) & ! shiny::isTruthy(input$anno_file1)) ){
                return(NULL)
            }
        })

        shinydashboard::box(width = 4,
                            shiny::conditionalPanel(
                                condition = 'input.emb_opt == "YES"',
                                shiny::radioButtons(inputId = 'cl_on_emb', label = 'Do you want to cluster the embedded matrix?',
                                                    choices = c('YES', 'NO'), selected = 'NO'),
                            ),
                            shiny::selectInput(inputId = 'cluster_method', label = 'Select a clustering method',
                                               choices = list(
                                                   'algorithms' = c("kmeans", "dbscan" , "hclust"),
                                                   'customize' = anno_list[!sapply(anno_list, is.null)]
                                               ),
                                               selected = "kmeans"
                            ),
                            shiny::conditionalPanel(
                                condition = "input.cluster_method == 'annotation' | input.cluster_method == 'annotation1' ",
                                shiny::uiOutput('custom_anno_selector'),
                                shiny::uiOutput('custom_anno_shape_elector')
                            ),
                            shiny::conditionalPanel(
                                condition = "input.cluster_method == 'kmeans'",
                                shiny::numericInput(inputId = "k",
                                                    label = "Number of clusters for kmeans",
                                                    min = 2, max= 500, value = 50
                                )
                            ),
                            shiny::conditionalPanel(
                                condition = "input.cluster_method == 'dbscan'",
                                shiny::sliderInput(inputId = "eps", label = "epsilon",
                                                   min = 0, max = 1,value = 0.1, step = 0.01
                                )
                            ),
                            shiny::conditionalPanel(
                                condition = "input.cluster_method == 'hclust'",
                                shiny::sliderInput("h", "Select h for cutree function:",
                                                   min = 0, max = 100, value = 10, step = 0.01
                                )
                            ),
                            shiny::uiOutput('shape_opt'),
                            shiny::hr(),
                            shiny::actionButton("cl_btn", "RUN")
        )
    })

    ############## plot cluster in 2D #############

    output$cl_plot_box <- shiny::renderUI({
        input$cl_btn
        shiny::isolate({
            if (is.null(shiny::reactiveValuesToList(clusters)$cluster)){
                return(NULL)
            }else{
                clusters <- shiny::reactiveValuesToList(clusters)
                clu_anno <- clusters$cluster %>% dplyr::arrange(clust)
                if (is.null(clu_anno$shape)){
                    cl_plot <- MoNETA::plot_2D_matrix(coord = clusters$dr_mat[1:2,], nodes_anno = clu_anno,
                                                      id_name = 'id', id_anno_color = 'clust',
                                                      interactive = TRUE, wo_legend = FALSE) %>%
                        plotly::layout(font = list(family = "LMRoman10", color = "black"),
                                       title = list(text = title(), font = list(family = "LMRoman10", color = "black", size = 20)),
                                       xaxis = list(title = list(text = paste0(input$dr_method, '1'))),
                                       yaxis = list(title = list(text =  paste0(input$dr_method, '2')))
                        ) %>%
                        plotly::config(
                            toImageButtonOptions = list(filename = paste('MoNETAshiny', input$dr_method, input$cluster_method, 'plot', sep = '_'),
                                                        format = "png", width = 1800, height = 800)
                        )
                } else{
                    cl_plot <- MoNETA::plot_2D_matrix(coord = clusters$dr_mat[1:2,], nodes_anno = clu_anno,
                                                      id_name = 'id', id_anno_color = 'clust', id_anno_shape = 'shape',
                                                      interactive = TRUE, wo_legend = FALSE) %>%
                        plotly::layout(font = list(family = "LMRoman10", color = "black"),
                                       title = list(text = title(), font = list(family = "LMRoman10", color = "black", size = 20)),
                                       xaxis = list(title = list(text = paste0(input$dr_method, '1'))),
                                       yaxis = list(title = list(text =  paste0(input$dr_method, '2')))
                        ) %>%
                        plotly::config(
                            toImageButtonOptions = list(filename = paste('MoNETAshiny', input$dr_method, input$cluster_method, 'plot', sep = '_'),
                                                        format = "png", width = 1800, height = 800)
                        )
                }
                return(shinydashboard::box(width = 8,
                                           cl_plot, shiny::hr(), shiny::downloadButton(outputId = 'download_cl_plot', label = 'Download plot')))
            }
        })
    })

    ### dendrogram
    output$plotHclust_box <- shiny::renderUI({
        if (!is.null(shiny::reactiveValuesToList(clusters)$tree)){
            output$plotHclust <- shiny::renderPlot({
                input$cl_btn
                shiny::isolate({
                    if (is.null(shiny::reactiveValuesToList(clusters)$num_clust)) {
                        tree <- shiny::reactiveValuesToList(clusters)$tree
                        tree %>%
                            dendextend::set("labels", '') %>%
                            #dendextend::set("branches_k_color", k = round(base::attr(tree, "height"), digits = 2)) %>%
                            plot(horiz=FALSE, axes=TRUE)
                    } else{
                        tree <- shiny::reactiveValuesToList(clusters)$tree
                        num_clust <- shiny::reactiveValuesToList(clusters)$num_clust
                        the_bars <- shiny::reactiveValuesToList(clusters)$cluster %>%
                            tibble::column_to_rownames(var = "id") %>%
                            dplyr::mutate_if(is.character, as.factor) %>%
                            dplyr::mutate_if(is.factor, as.numeric)
                        par(mar = c(10, 3, 3, 4) + 0.1,
                            xpd = NA) # allow content to go into outer margin
                        tree_plot <- tree %>%
                            #dendextend::set("labels_col", k = num_clust) %>%
                            dendextend::set("labels", '') %>%
                            dendextend::set("branches_k_color", k = num_clust) %>%
                            plot(horiz=FALSE, axes=TRUE)
                        final_tree_plot <- tree_plot +
                            abline(h = input$h, lty = 'dotdash') +
                            dendextend::colored_bars(colors = the_bars, dend = tree,
                                         rowLabels = colnames(the_bars),
                                         #y_shift = -10
                            )
                        return(final_tree_plot)
                    }
                })
            })
            shinydashboard::box(width = 8, shiny::plotOutput(outputId = "plotHclust"))
        }else {
            return(NULL)
        }
    })

    ### table

    output$cluster_table_box <- shiny::renderUI({
        if (!is.null(shiny::reactiveValuesToList(clusters)$cluster)){
            output$cluster_table <- DT::renderDT(server = FALSE, {
                input$cl_btn
                shiny::isolate({
                    DT::datatable(shiny::reactiveValuesToList(clusters)$cluster,
                                  extensions = 'Buttons',
                                  options = list(dom = "Blfrtip",
                                                 buttons = list("copy",
                                                                list(extend = 'csv',   filename =  paste("MoNETA", input$cluster_method, 'table',sep = "_")),
                                                                list(extend = 'excel', filename =  paste("MoNETA", input$cluster_method, 'table',sep = "_"))
                                                 ),
                                                 paging = TRUE,
                                                 pageLength = 10,
                                                 scrollX = TRUE,
                                                 scrollY = TRUE,
                                                 autoWidth = TRUE
                                  ),
                                  selection = 'multiple',
                                  filter = 'top',
                                  rownames = FALSE
                    )
                })
            })
            shinydashboard::box(width = 4, DT::DTOutput('cluster_table'))
        } else{
            return(NULL)
        }
    })

    ############## download plot cluster in 2D #############

    output$download_cl_plot <- shiny::downloadHandler(
        filename = function() { paste('MoNETAshiny', input$dr_method, input$cluster_method, 'plot.png', sep='_') },
        content = function(file) {
            ggplot2::ggsave(filename = file, plot = cl_plot(), device = "png", width = 10, height = 5, limitsize = FALSE)
        }
    )

}

