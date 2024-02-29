#' MoNETA Shiny app
#'
#' @import shiny
#' @import shinydashboard
#' @import shinyalert
#' @import shinyFiles
#' @import shinyMatrix
#' @import shinyjs
#' @import shinycssloaders
#' @import readr
#' @import glue
#' @import conflicted
#' @import fpc
#' @import parallel
#' @import plotly
#' @import ggplot2
#' @import visNetwork
#' @import tidyverse
#' @import dendextend
#' @import magrittr
#' @import dplyr
#' @import igraph
#' @import vroom
#' @import zip
#' @importFrom utils data read.csv
#' @importFrom stats as.dendrogram dist kmeans hclust
#' @importFrom graphics abline par
#' @param MAXreq shiny max request size
#' @return Shiny app
#' @export


MoNETAshiny = function(MAXreq = 10000) {
    options(shiny.maxRequestSize = MAXreq * 1024^2)
    shiny::shinyApp(ui, server)

}


netSummary <- function(network){

    nodes <- c(network$source, network$target)
    degree <- table(nodes) %>%
        dplyr::as_tibble(.) %>%
        dplyr::arrange(desc(n))
    top10 <- c(degree[c(1:10), 'nodes'])
    netSummary_list <- list( 'dim'= dim(network),
                             'edges'= dim(network)[1],
                             'nodes'= length(unique(nodes)),
                             'maxConnections'= degree[[1, 'n']]#,
                             #'maxDegreeNodes'= top10$nodes
    )
    return(netSummary_list)
}


ui <- shinydashboard::dashboardPage(

    shinydashboard::dashboardHeader(title = shiny::span("MoNETA ", style = "color: white; font-size: 28px")),

    shinydashboard::dashboardSidebar(
        shinyjs::useShinyjs(),
        shinydashboard::sidebarMenu(id = 'tabs',
                                    shiny::tags$head(shiny::tags$link(rel = "stylesheet", type = "text/css", href = "custom.css")),
                                    shiny::tags$head(shiny::tags$style(".inactiveLink {
                            pointer-events: none;
                           cursor: default;
                           }")),
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
        shinydashboard::tabItems(shinydashboard::tabItem(tabName = "mat_sub_1",
                                                         shiny::fluidRow(
                                                             shiny::column(width = 6,
                                                                           shiny::uiOutput('jump2P1.1')),
                                                             shiny::column(width = 6,
                                                                           shiny::uiOutput('jump2P2'), align = 'right'
                                                             )
                                                         ),

                                                         shiny::fluidRow(

                                                             shiny::column(width = 6,
                                                                           shiny::fluidRow(
                                                                               shiny::column(width = 12,
                                                                                             shiny::tags$hr(),
                                                                                             shinydashboard::box(
                                                                                                 title = shiny::h2(shiny::span("Upload omics matrices", style = "font-weight: bold")),
                                                                                                 width = 12, status = "primary", solidHeader = FALSE,
                                                                                                 shiny::radioButtons(inputId = 'omics_example_opt',
                                                                                                                     label = shiny::h4(shiny::span('Load the example dataset', style = "font-weight: bold")),
                                                                                                                     choices = c('Yes', 'No'), selected = 'No'),
                                                                                                 shiny::conditionalPanel(
                                                                                                     condition = 'input.omics_example_opt == "Yes"',
                                                                                                     shiny::checkboxGroupInput(inputId = 'omics_example_files',
                                                                                                                               label = shiny::h4(shiny::span('Select one or more omics matrices', style = "font-weight: bold")),
                                                                                                                               choices = c("GliomaCNV_norm","GliomaMethylation_norm", "GliomaExpression_norm" ),
                                                                                                                               selected = c("GliomaCNV_norm","GliomaMethylation_norm", "GliomaExpression_norm" )),
                                                                                                     shiny::tags$hr()
                                                                                                 ),
                                                                                                 shiny::conditionalPanel(
                                                                                                     condition = 'input.omics_example_opt == "No"',
                                                                                                     shiny::fileInput(inputId = "omics_files", label = shiny::h4(shiny::span('Select one or more files', style = "font-weight: bold")),
                                                                                                                      multiple = TRUE),
                                                                                                     shiny::uiOutput('omics_names')
                                                                                                 ),
                                                                                                 shiny::actionButton('load_mat_button', 'Load')
                                                                                             )
                                                                               ),

                                                                               shiny::column(width = 12,
                                                                                             shiny::htmlOutput(outputId = "omics_sum_info")
                                                                               ),

                                                                               shiny::column(width = 12,
                                                                                             shinydashboard::box(
                                                                                                 width = 12, status = "primary", solidHeader = FALSE,
                                                                                                 title = shiny::h2(shiny::span("Upload annotation file", style = "font-weight: bold")),
                                                                                                 shiny::radioButtons(inputId = 'anno_example_opt',
                                                                                                                     label = shiny::h4(shiny::span('Load the example annotation file', style = "font-weight: bold")),
                                                                                                                     choices = c('Yes', 'No'), selected = 'No'),
                                                                                                 shiny::conditionalPanel(
                                                                                                     condition = 'input.anno_example_opt == "No"',
                                                                                                     shiny::fileInput(inputId = "anno_file", label = shiny::h4(shiny::span('Select a file', style = "font-weight: bold")),
                                                                                                                      multiple = FALSE)
                                                                                                 ),
                                                                                                 shiny::tags$hr(),
                                                                                                 shiny::actionButton('load_anno_button', 'Load')
                                                                                             )
                                                                               ),

                                                                               shiny::column(width = 12,
                                                                                             shiny::htmlOutput(outputId = "anno_info")
                                                                               ),

                                                                               shiny::column(width = 12,
                                                                                             shiny::htmlOutput(outputId = "anno_info_extra")
                                                                               )
                                                                           )
                                                             ),
                                                             shiny::column(width = 6,
                                                                           shiny::tags$hr(),
                                                                           shinydashboard::infoBox(
                                                                               title = shiny::h3(shiny::span('Information', style = "font-weight: bold")),
                                                                               value = shiny::uiOutput('info_box_1'),
                                                                               subtitle = NULL,
                                                                               icon = shiny::icon("info"), color = "aqua", width = NULL,
                                                                               fill = FALSE)
                                                             )
                                                         )
        ),
        shinydashboard::tabItem(tabName = "mat_sub_2",
                                shiny::fluidRow(
                                    shiny::column(width = 6,
                                                  shiny::uiOutput('back2P1')),
                                    shiny::column(width = 6,
                                                  shiny::uiOutput('jump2P3'),  align = 'right'
                                    )
                                ),
                                shiny::fluidRow(
                                    shiny::column(width = 6,
                                                  shiny::fluidRow(
                                                      shiny::tags$hr(),
                                                      shiny::column(width = 12,
                                                                    #shiny::h2(shiny::span("Omics matrices Pre-processing", style = "font-weight: bold")),
                                                                    shiny::uiOutput("process_omics_mat"),
                                                      ),
                                                      shiny::column(width = 12,
                                                                    shiny::uiOutput('pro_matrices_sum_info')
                                                      ),
                                                      shiny::column(width = 12,
                                                                    shiny::uiOutput('intersection')
                                                      ),
                                                      shiny::column(width = 12,
                                                                    shiny::uiOutput('intersection_info')
                                                      )
                                                  )
                                    ),
                                    shiny::column(width = 6,
                                                  shiny::tags$hr(),
                                                  shiny::fluidRow(
                                                            shinydashboard::infoBox(
                                                                title = shiny::h3(shiny::span('Information', style = "font-weight: bold")),
                                                                value = shiny::uiOutput('info_box_2'),
                                                                subtitle = NULL,
                                                                icon = shiny::icon("info"), color = "aqua", width = NULL,
                                                                fill = FALSE),
                                                            shiny::uiOutput('download_proc_mat_box')
                                                  )
                                    )
                                )
        ),

        shinydashboard::tabItem(tabName = 'net_sub_1',
                                shiny::fluidRow(
                                    shiny::column(width = 6,
                                                  shiny::uiOutput('back2P2')),
                                    shiny::column(width = 6,
                                                  shiny::uiOutput('jump2P4'),  align = 'right'
                                    )
                                ),
                                shiny::fluidRow(
                                    shiny::column(width = 6,
                                                  shiny::tags$hr(),
                                                  shiny::h2(shiny::span("Omics Network Inference", style = "font-weight: bold")),
                                                  shiny::fluidRow(
                                                      shiny::uiOutput("omics_net_arguments"),
                                                      shiny::uiOutput("download_net_box")
                                                  )
                                    ),
                                    shiny::column(width = 6,
                                                  shiny::tags$hr(),
                                                  shiny::column(width = 12,
                                                                shinydashboard::infoBox(
                                                                    title = shiny::h3(shiny::span('Information', style = "font-weight: bold")),
                                                                    value = shiny::uiOutput('info_box_3'),
                                                                    subtitle = NULL,
                                                                    icon = shiny::icon("info"), color = "aqua", width = NULL,
                                                                    fill = FALSE)

                                                  )
                                    )
                                )
        ),
        shinydashboard::tabItem(tabName = 'net_sub_2',
                                shiny::fluidRow(
                                    shiny::column(width = 6,
                                                  shiny::uiOutput('back2P3')),
                                    shiny::column(width = 6,
                                                  shiny::uiOutput('jump2P5'),  align = 'right'
                                    )
                                ),
                                shiny::fluidRow(
                                    shiny::column(width = 6,
                                                  shiny::tags$hr(),
                                                  shiny::h2(shiny::span("Omics Network Filtering", style = "font-weight: bold")),
                                                  shiny::fluidRow(
                                                      shiny::uiOutput("plot_net_box")
                                                  )
                                    ),
                                    shiny::column(width = 6,
                                                  shiny::tags$hr(),
                                                  shiny::column(width = 12,
                                                                shinydashboard::infoBox(
                                                                    title = shiny::h3(shiny::span('Information', style = "font-weight: bold")),
                                                                    value = shiny::uiOutput('info_box_4'),
                                                                    subtitle = NULL,
                                                                    icon = shiny::icon("info"), color = "aqua", width = NULL,
                                                                    fill = FALSE),
                                                                shiny::uiOutput("download_fnet_box")
                                                  )
                                    )
                                )
        ),

        shinydashboard::tabItem(tabName = "rwr_tab_1",

                                shiny::fluidRow(
                                    shiny::column(width = 6,
                                                  shiny::uiOutput('back2P1_from1.1')
                                    ),
                                    shiny::column(width = 6,
                                                  shiny::uiOutput('jump2P5.1'),  align = 'right'
                                    )
                                ),

                                shiny::fluidRow(
                                    shiny::column(width = 6,
                                                  shiny::fluidRow(
                                                      shiny::column(width = 12,
                                                                    shiny::tags$hr(),
                                                                    shinydashboard::box(
                                                                        title = shiny::h2(shiny::span("Upload omics networks", style = "font-weight: bold")),
                                                                        width = 12, status = "primary", solidHeader = FALSE,
                                                                        shiny::fileInput("omics_net_files", label = shiny::h4(shiny::span( 'Select one or more files', style = "font-weight: bold")),
                                                                                         multiple = TRUE,
                                                                                         accept = c("text/csv", '.RDS', "text/comma-separated-values,text/plain", ".csv")
                                                                        ),
                                                                        shiny::uiOutput('omics_net_names'),
                                                                        shiny::actionButton('load_omics_net_button', 'Load')
                                                                    )
                                                      ),
                                                      shiny::column(width = 12,
                                                                    shiny::htmlOutput(outputId = "omics_net_sum_info")
                                                      ),
                                                      shiny::column(width = 12,
                                                                    shinydashboard::box(
                                                                        title = shiny::h2(shiny::span("Upload annotation file", style = "font-weight: bold")),
                                                                        width = 12, status = "primary", solidHeader = FALSE,
                                                                        shiny::fileInput(inputId = "anno_file1",
                                                                                         label = shiny::h4(shiny::span( 'Select a file', style = "font-weight: bold")),
                                                                                         multiple = FALSE),
                                                                        shiny::actionButton('load_anno_button1', 'Load')
                                                                    ),
                                                                    shiny::htmlOutput(outputId = "anno_info1")
                                                      )

                                                  )
                                    ),
                                    shiny::column(width = 6,
                                                  shiny::tags$hr(),
                                                  shinydashboard::infoBox(
                                                      title = shiny::h3(shiny::span('Information', style = "font-weight: bold")),
                                                      value = shiny::uiOutput('info_box_5'),
                                                      subtitle = NULL,
                                                      icon = shiny::icon("info"), color = "aqua", width = NULL,
                                                      fill = FALSE),
                                                  shiny::htmlOutput(outputId = "anno_info_extra1"),
                                                  shiny::htmlOutput(outputId = "anno_info_extra2")
                                    )

                                )
        ),

        shinydashboard::tabItem(tabName = "rwr_tab_2",
                                shiny::fluidRow(
                                    shiny::column(width = 6,
                                                  shiny::uiOutput('back2P4'),
                                                  shiny::uiOutput('back2P1.1')),
                                    shiny::column(width = 6,
                                                  shiny::uiOutput('jump2P6'),  align = 'right'
                                    )
                                ),
                                shiny::fluidRow(
                                    shiny::column(width = 6,
                                                  shiny::fluidRow(
                                                      shiny::column(width=12,
                                                                    shiny::tags$hr(),
                                                                    shinydashboard::box(width = 12, status = "primary", solidHeader = FALSE,
                                                                                        title =  shiny::h2(shiny::span("Construction of Multiplex network", style = "font-weight: bold")),
                                                                                        shiny::radioButtons(inputId = 'weightMultiplex',
                                                                                                            label = shiny::h4(shiny::span('Do you want to create a weighted multiplex network?', style = "font-weight: bold")),
                                                                                                            choices = c("YES", "NO"), selected = "NO"),
                                                                                        shiny::tags$hr(),
                                                                                        shiny::conditionalPanel(
                                                                                            condition = 'input.weightMultiplex == "YES"',
                                                                                            shiny::radioButtons(inputId = 'pruneMultiplex',
                                                                                                                label = shiny::h4(shiny::span('Do you want to prune the multiplex network?', style = "font-weight: bold")),
                                                                                                                choices = c('YES', 'NO'), selected = 'NO'),
                                                                                            shiny::conditionalPanel(
                                                                                                condition = 'input.pruneMultiplex == "YES"',
                                                                                                shiny::numericInput(inputId = 'pruneMultiplex_th',
                                                                                                                    label = shiny::h4(shiny::span('Select a threshold to prune the multiplex network', style = "font-weight: bold")),
                                                                                                                    min = 0, max = 100, value = 50)
                                                                                            ),
                                                                                            shiny::tags$hr()
                                                                                        ),
                                                                                        shiny::actionButton(inputId = 'gen_multiplex_btn', label = 'GO'))
                                                      ),
                                                      shiny::column(width=12,
                                                                    shiny::htmlOutput(outputId = "multi_net_sum_info")
                                                      )
                                                  )
                                    ),
                                    shiny::column(width = 6,
                                                  shiny::tags$hr(),
                                                  shinydashboard::infoBox(
                                                      title = shiny::h3(shiny::span('Information', style = "font-weight: bold")),
                                                      value = shiny::uiOutput('info_box_6'),
                                                      subtitle = NULL,
                                                      icon = shiny::icon("info"), color = "aqua", width = NULL,
                                                      fill = FALSE)
                                    )
                                )
        ),
        shinydashboard::tabItem(tabName = "rwr_tab_3",
                                shiny::fluidRow(
                                    shiny::column(width = 6,
                                                  shiny::uiOutput('back2P5')),
                                    shiny::column(width = 6,
                                                  shiny::uiOutput('jump2P7'),  align = 'right'
                                    )
                                ),

                                shiny::fluidRow(

                                    shiny::column(width = 6,
                                                  shiny::tags$hr(),
                                                  shinydashboard::box(width=12,
                                                                      status = "primary", solidHeader = FALSE,
                                                                      title = shiny::h2(shiny::span("Random Walk with Restart", style = "font-weight: bold")),
                                                                      shinydashboard::tabBox(width=12,
                                                                                             title = shiny::column(width = 6, shiny::uiOutput('download_rwr_mat')),
                                                                                             shiny::tabPanel(title = shiny::h3(shiny::span('Restart', style = "font-weight: bold")),
                                                                                                             shiny::fluidRow(shiny::column(width = 12,
                                                                                                                                           shiny::sliderInput('restart', shiny::h4(shiny::span('Select the restarting parameter (r)', style = "font-weight: bold")),
                                                                                                                                                              min = 0, max = 1, value = 0.7, step = 0.01)
                                                                                                             )
                                                                                                             ),
                                                                                                             shiny::tags$hr(),
                                                                                                             shiny::fluidRow(shiny::column(width = 12,
                                                                                                                                           shiny::radioButtons('tao_opt',
                                                                                                                                                               label = shiny::h4(shiny::HTML('<span style:"font-family = LM Roman 10"><b>Restarting Probabilities (&tau;) </b></span>')),
                                                                                                                                                               choices = c('Custumize restarting probabilities per layer', 'Use default restarting probabilities'),
                                                                                                                                                               selected = 'Use default restarting probabilities'),
                                                                                                                                           shiny::uiOutput(outputId = "tauBIO")
                                                                                                             )
                                                                                                             )
                                                                                             ),
                                                                                             shiny::tabPanel(title = shiny::h3(shiny::span('Transition', style = "font-weight: bold")),
                                                                                                             shiny::fluidRow(
                                                                                                                 shiny::column(width = 12,
                                                                                                                               shiny::sliderInput('delta',
                                                                                                                                                  shiny::h4(shiny::HTML('<span style:"font-family = LM Roman 10"><b>Select the transition parameter (&delta;) </b></span>')),
                                                                                                                                                  min = 0, max = 1, value = 0.5, step = 0.01)
                                                                                                                 )
                                                                                                             ),
                                                                                                             shiny::tags$hr(),
                                                                                                             shiny::fluidRow(
                                                                                                                 shiny::column(width = 12,
                                                                                                                               shiny::radioButtons('omics_jump_neigh',
                                                                                                                                                   shiny::h4(shiny::span('RWR-NF option', style = "font-weight: bold")),
                                                                                                                                                   c('YES', 'NO'), selected = 'NO'),
                                                                                                                               shiny::conditionalPanel(
                                                                                                                                   condition = 'input.omics_jump_neigh == "YES"',
                                                                                                                                   shiny::tags$hr(),
                                                                                                                                   shiny::radioButtons('omics_weight_jump_neigh',
                                                                                                                                                       shiny::h4(shiny::span('Weight inter-layers edges', style = "font-weight: bold")),
                                                                                                                                                       c('YES', 'NO'), selected = 'YES')
                                                                                                                               )
                                                                                                                 )
                                                                                                             ),
                                                                                                             shiny::tags$hr(),
                                                                                                             shiny::fluidRow(
                                                                                                                 shiny::column(width = 12,
                                                                                                                               shiny::h4(shiny::span('Transition layer matrices', style = "font-weight: bold")),
                                                                                                                               shiny::radioButtons(inputId = 'bioInf_transition',
                                                                                                                                                   label = 'Do you want to create a biologically informed layer transition matrix?',
                                                                                                                                                   choices = c('YES', 'NO'), selected = 'NO'),
                                                                                                                               shiny::conditionalPanel(
                                                                                                                                   condition = 'input.bioInf_transition == "NO"',
                                                                                                                                   shiny::uiOutput('omics_trans_mat')
                                                                                                                               )
                                                                                                                 )
                                                                                                             )
                                                                                             )
                                                                      ),

                                                                      shiny::fluidRow(
                                                                          shiny::column(width = 12,
                                                                                        shinydashboard::box(width = 12,
                                                                                                            shiny::numericInput('cores_rwr',
                                                                                                                                shiny::h4(shiny::span('Select the number of cores to run the RWR', style = "font-weight: bold")),
                                                                                                                                min = 1, max = 190, value = 1)
                                                                                        )
                                                                          )
                                                                      ),
                                                                      shiny::fluidRow(
                                                                          shiny::tags$hr(),
                                                                          shiny::column(width = 6, shiny::actionButton('rwr_button', 'Run RWR')),
                                                                          #shiny::column(width = 6, shiny::uiOutput('download_rwr_mat')),
                                                                      )
                                                  )
                                    ),
                                    shiny::column(width = 6,
                                                  shiny::tags$hr(),
                                                  shinydashboard::infoBox(
                                                      title = shiny::h3(shiny::span('Information', style = "font-weight: bold")),
                                                      value = shiny::uiOutput('info_box_7'),
                                                      subtitle = NULL,
                                                      icon = shiny::icon("info"), color = "aqua", width = NULL,
                                                      fill = FALSE)
                                    )
                                )
        ),

        shinydashboard::tabItem(tabName = "dr_tab",
                                shiny::fluidRow(
                                    shiny::column(width = 6,
                                                  shiny::uiOutput('back2P6')),
                                    shiny::column(width = 6,
                                                  shiny::uiOutput('jump2P8'),  align = 'right'
                                    )
                                ),
                                shiny::h2(shiny::span("Dimensionality Reduction Methods", style = "font-weight: bold")),
                                shiny::fluidRow(
                                    shiny::uiOutput("dr_opt"),
                                    shiny::uiOutput('dr_plot_box'),
                                    shiny::uiOutput('download_dr_box')
                                )
        ),

        shinydashboard::tabItem(tabName = "cl_tab",
                                shiny::fluidRow(
                                    shiny::column(width = 6,
                                                  shiny::uiOutput('back2P7'))
                                ),
                                shiny::h2(shiny::span("Clustering Analysis", style = "font-weight: bold")),
                                shiny::fluidRow(
                                    shiny::uiOutput("cl_opt"),
                                    shinycssloaders::withSpinner(shiny::uiOutput(outputId = "cl_plot_box"))
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

    shiny::tags$head(shiny::tags$style(shiny::HTML('* {font-family: "LM Roman 10"}; font-size: 1em;')))

)


server <- function(input, output, session) {

    ############################################################################################
    #                            OMICS DATA: LOADING DATA (INPUT)                       #
    ############################################################################################

    output$info_box_1 <- renderText({
        HTML("<br/> <span style='font-weight:normal;'>
            <p align='justify'>
            Welcome to the MoNETA Shiny app, an R package for network-based multi-omics data integration.<br/>
            You can start the pipeline by uploading omics matrices on the following page,
            or click the <b> Skip </b> button in the top left-hand corner to directly load patient networks. </span> <br/>
            </p>
            <hr style='border-top: 1px solid white;'>

            <b> Upload Omic Matrices </b>: <br/> <span style='font-weight:normal;'>
            <p align='justify'>
            These matrices should have features on rows and the same or non-overlapping sets of samples on columns.
            Accepted formats include .csv, .txt, and .RDS. </span> <br/>
            </p>
            <hr style='border-top: 1px solid white;'>

            <b> Use Example Dataset </b>: <br/>  <span style='font-weight:normal;'>
             <p align='justify'>
            Alternatively, you can upload the example dataset consisting of three omic matrices related to the analysis of 788 glioma patients,
            generated by <span style= 'color: blue;'> <u>  <a href='https://www.ceccarellilab.org/'  target='_blank'> Ceccarelli Lab</a></u>. </span> </span>
            </p>
            <hr style='border-top: 1px solid white;'>

            <b> Upload Annotation File </b>: <br/> <span style='font-weight:normal;'>
            <p align='justify'>
            This step is not mandatory but useful to enrich subsequent plots.
            The file should be a dataframe with a variable number of columns,
            where one column must contain the sample IDs, the same ones present as column names in the uploaded files.
            Accepted formats include .csv, .txt, and .RDS. </span> <br/>
            </p>
            <hr style='border-top: 1px solid white;'>

            <b> Use Example Annotation File </b>: <br/>  <span style='font-weight:normal;'>
            <p align='justify'>
            If you have selected the Example Dataset, you can use the related annotation file.
            The 788 samples are classified based on three different criteria: <b>histology</b>, genetic background (<b>IDH-status</b>)
            and DNA methylation profiles (<b>Supervised.DNA.Methylation.Cluster</b>).
            The sampleIDs are stored in the <b>Case</b> column.
            For further details, please refer to the following paper
            <span style= 'color: blue;'> <u>  <a href='https://www.ceccarellilab.org/'  target='_blank'> (link)</a></u></span>. </span>
            </p>
            <hr style='border-top: 1px solid white;'>")
    })


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
                    shinyalert::shinyalert("Type Error", "Uploaded Data is not a numeric matrix", closeOnClickOutside = TRUE, type = "error")
                    returnValue()
                }
            }

        } else {
            data("GBM_mtx", package = 'MoNETA')
            req(input$omics_example_files)
            example_data <- GBM_mtx[input$omics_example_files]
            omicsFiles <- example_data
        }

        omicsFiles
    })


    observeEvent(input$load_mat_button, {
        shinyalert::shinyalert(
            title = "Wait",
            text = "Waiting for data loading",
            size = "xs",
            closeOnEsc = TRUE,
            closeOnClickOutside = TRUE,
            html = TRUE,
            type = "info",
            showConfirmButton = TRUE,
            confirmButtonText = "OK",
            confirmButtonCol = "#004192",
            showCancelButton = FALSE,
            imageUrl = "",
            animation = TRUE
        )

    },
    ignoreNULL = FALSE,
    ignoreInit = TRUE
    )

    ### Load annotation file
    annotation <- reactiveVal(NULL)
    shiny::observeEvent(input$load_anno_button, {
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
                    annotation(listFiles)
                }else {
                    shinyalert::shinyalert("Type Error", "Uploaded Data is not a dataframe",closeOnClickOutside = TRUE, type = "error")
                    returnValue()
                }
            }
        }else if (input$anno_example_opt == 'Yes'){
            data("GBM_pdata", package = 'MoNETA')
            listFiles[['annotation']] <- GBM_pdata
            annotation(listFiles)
        }else{
            return(NULL)
        }
    })

    observeEvent(input$load_anno_button, {
        shinyalert::shinyalert(
            title = "Wait",
            text = "Waiting for data loading",
            size = "xs",
            closeOnEsc = TRUE,
            closeOnClickOutside = TRUE,
            html = TRUE,
            type = "info",
            showConfirmButton = TRUE,
            confirmButtonText = "OK",
            confirmButtonCol = "#004192",
            showCancelButton = FALSE,
            imageUrl = "",
            animation = TRUE
        )

    },
    ignoreNULL = FALSE,
    ignoreInit = TRUE
    )

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



    ###

    output$jump2P2 <- renderUI({
        if (length(omics_files()) != 0){
            actionButton('jump2P2', label = 'Next', icon("paper-plane"),
                         style="color: #fff; background-color: #337ab7; border-color: #2e6da4; padding:20px")
        }
    })

    observeEvent(input$jump2P2, {
        shinydashboard::updateTabItems(session, inputId = "tabs", selected = "mat_sub_2")
        shinyjs::runjs('$(".sidebar-menu .treeview").removeClass("active"); $("#mat_tab").closest(".treeview").addClass("active");')
    })

    ###
    output$jump2P1.1 <- renderUI({
        actionButton('jump2P1.1', label = 'Skip', icon("paper-plane"),
                     style="color: #fff; background-color: #337ab7; border-color: #2e6da4; padding:20px")
    })

    observeEvent(input$jump2P1.1, {
        shinydashboard::updateTabItems(session, inputId = "tabs", selected = "rwr_tab_1")
    })


    ############################################################################################
    #                        OMICS DATA: PRE-PROCESSING (INPUT)                      #
    ############################################################################################

    shinyjs::addCssClass(selector = "a[data-value='mat_sub_2']", class = "inactiveLink")
    observe({
        if (length(omics_files()) != 0){
            shinyjs::removeCssClass(selector = "a[data-value='mat_sub_2']", class = "inactiveLink")
        }
    })

    ############
    output$info_box_2 <- renderText({
        HTML("<br/> <span style='font-weight:normal;'>
        <p align='justify'>
            Welcome to the MoNETA Shiny app, an R package for network-based multi-omics data integration.<br/>
            In this section you can manage the omics matrices before inferring the networks. </span> <br/>
            </p> <hr style='border-top: 1px solid white;'>

            <b> Omics Matrices Pre-processing </b>: <br/> <span style='font-weight:normal;'>
            <p align='justify'>
            Here, you can decide for each omics matrix whether to remove samples with no detected features  and/or normalize them by column.
            When the <b> Submit </b> button is pressed, the omics matrices will be processed <u>all at once</u> </span> <br/>
            </p> <hr style='border-top: 1px solid white;'>

            <b> Omics Matrices Intersection </b>: <br/>  <span style='font-weight:normal;'>
            <p align='justify'>
            This box will be shown only if you have loaded more than one omics matrix and the <b> Submit </b> button have already been  pressed.
            If you select <b> Yes</b>, only the overlapping sets of samples will be kept.</span>
            </p> <hr style='border-top: 1px solid white;'>

            <span style='font-weight:normal;'>
            <p align='justify'>
            The processed matrices can be downloaded from the <b>download box</b> displayed below. </span>
            </p> <hr style='border-top: 1px solid white;'>
            ")
    })

    ############

    intersected_omics_mat <- shiny::reactiveValues()
    shiny::observeEvent(input$intersect_btn, {
        processed_mats <- shiny::reactiveValuesToList(proc_matrices)
        #shiny::isolate({
        if (input$intersect_opt == 'YES'){
            intersected_omics_mat$matrices <- MoNETA::get_intersection_matrices(processed_mats)
        }else {
            intersected_omics_mat$matrices <- processed_mats
        }

        shinyalert::shinyalert(
            title = "Success",
            text = 'Intersection step done! <br/>  Now press <b> "Next" </b> in the top right-hand corner to continue.',
            closeOnEsc = TRUE,
            closeOnClickOutside = TRUE,
            html = TRUE,
            type = "success",
            showConfirmButton = TRUE,
            confirmButtonText = "OK",
            confirmButtonCol = "#004192",
            showCancelButton = FALSE,
            imageUrl = "",
            animation = TRUE
        )

        #})
    })



    ############################################################################################
    #                        NETWORK : NETWORK INFERENCE (INPUT)                      #
    ############################################################################################

    shinyjs::addCssClass(selector = "a[data-value='net_sub_1']", class = "inactiveLink")
    observe({
        if (length(omics_files()) == 1 && length(reactiveValuesToList(proc_matrices)) != 0){
            shinyjs::removeCssClass(selector = "a[data-value='net_sub_1']", class = "inactiveLink")
        } else if (length(omics_files()) >= 1 && length(reactiveValuesToList(proc_matrices)) != 0 && shiny::isTruthy(input$intersect_btn)){
            shinyjs::removeCssClass(selector = "a[data-value='net_sub_1']", class = "inactiveLink")
        }else{
            return(NULL)
        }

    })

    ############
    output$info_box_3 <- renderText({
        HTML("<br/> <span style='font-weight:normal;'>
        <p align='justify'>
            Welcome to the MoNETA Shiny app, an R package for network-based multi-omics data integration.<br/>
            In this section you can create the networks starting from each omics matrix. </span> <br/>
            </p> <hr style='border-top: 1px solid white;'>

            <b> Omics Network Inference </b>: <br/> <span style='font-weight:normal;'>
            <p align='justify'>
            Here, you can set the parameters of MoNETA <b>k_star_net</b> function: <br/>
            i) <b>distance metrics</b>; <br/>
            ii) the <b>sparsity</b> of the data in the features space; <br/>
            iii) the maximum number of neighbors per node (<b>knn</b>);<br/>
            iv) the number of <b>cores</b> for parallelization; <br/>
            v) dynamic or fixed choice of number of neighbors per node (<b>k_star</b>). <br/>
            *Euclidean and Cosine distances are recommended for continuous data,
            while Manhattan distance is suitable for discrete data.
            </p> <hr style='border-top: 1px solid white;'>

            <span style='color: red;'>
            <p align='justify'>
            After clicking the button to infer a network, please wait for the process to finish before clicking the next one.
            If you click them one after the other, the second process will be queued, and the <b>Wait</b> message will be displayed accordingly. </span>
            </p> <hr style='border-top: 1px solid white;'> </span>

            <b>Network Visualization:</b> <br/>  <span style='font-weight:normal;'>
            <p align='justify'>
            Once the network is created, the plot is shown. The nodes are colored according to their nodeID by default,
            but, if an annotation file is loaded, the color and the shape can be updated by clicking the <b<Update plot</b> button.</span>
            </p> <hr style='border-top: 1px solid white;'>

            <span style='font-weight:normal;'>
            <p align='justify'>
            The networks can be downloaded from the <b>download box</b> located on the left-hand side. </span>
            </p> <hr style='border-top: 1px solid white;'>
            ")
    })


    ############################################################################################
    #                        NETWORK : NETWORK FILTERING (INPUT)                      #
    ############################################################################################

    shinyjs::addCssClass(selector = "a[data-value='net_sub_2']", class = "inactiveLink")
    observe({
        if (length(shiny::reactiveValuesToList(gene_networks)) == length(omics_files())){
            shinyjs::removeCssClass(selector = "a[data-value='net_sub_2']", class = "inactiveLink")
        }else{
            return(NULL)
        }
    })

    output$info_box_4 <- renderText({
        HTML("<br/> <span style='font-weight:normal;'>
            <p align='justify'>
            Welcome to the MoNETA Shiny app, an R package for network-based multi-omics data integration.<br/>
            Here you can prune the edges of inferred networks by filtering the <b>weight</b> column. </span> <br/>
            </p> <hr style='border-top: 1px solid white;'>

            <b> Show button </b>: <br/> <span style='font-weight:normal;'>
            <p align='justify'>
            When the <b>Show</b> button is pressed, the network and the related
            edge weight and node degree distributions (based on the chosen filtering threshold) will be displayed.
            This action does not trigger the update of the networks. </span> <br/>
            </p> <hr style='border-top: 1px solid white;'>

            <b> Submit button </b>: <br/> <span style='font-weight:normal;'>
            <p align='justify'>
            When the <b>Submit</b> button is pressed, the networks will be updated
            <u>all at once</u>, allowing you to proceed. </span> <br/>
            </p> <hr style='border-top: 1px solid white;'>

            <span style='font-weight:normal;'>
            <p align='justify'>
            The filtered networks can be downloaded from the <b>download box</b> displayed below. </span>
            </p> <hr style='border-top: 1px solid white;'>
            ")
    })


    ############################################################################################
    #                                 RWR: LOADING NETWORKS (INPUT)                            #
    ############################################################################################

    output$info_box_5 <- renderText({
        HTML("<br/> <span style='font-weight:normal;'>
            <p align='justify'>
            Welcome to the MoNETA Shiny app, an R package for network-based multi-omics data integration.<br/>
            You can start the pipeline here by uploading omics networks,
            or click the <b> Back </b> button in the top left-hand corner to return back to the 'Omics Matrices Uploading' page. </span> <br/>
            </p> <hr style='border-top: 1px solid white;'>

            <b> Upload Omics Networks </b>: <br/> <span style='font-weight:normal;'>
            <p align='justify'>
            The networks are undirected and are loaded as tables with three columns: source, target, and weight.
            The first two columns contain the samples/nodes, while the last one contains the weights of the edges.
            The higher the weight, the stronger the interaction. A 'weight' column containing only '1' represents an unweighted network.</span> <br/>
            </p> <hr style='border-top: 1px solid white;'>

            <b> Upload Annotation File </b>: <br/> <span style='font-weight:normal;'>
            <p align='justify'>
            This step is not mandatory but useful to enrich subsequent plots.
            The file should be a dataframe with a variable number of columns,
            where one column must contain the sample IDs representing the nodes in the uploaded files. </span> <br/>
            </p> <hr style='border-top: 1px solid white;'>

            <span style='font-weight:normal;'>
            <p align='justify'>
            Accepted formats include .csv, .txt, and .RDS. </span>
             </p> <hr style='border-top: 1px solid white;'>
            ")
    })


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
                    shinyalert::shinyalert("Type Error", "Uploaded Data is not a dataframe",closeOnClickOutside = TRUE, type = "error")
                    returnValue()
                }else if (ncol(file) != 3){
                    shinyalert::shinyalert("Column Error", "Uploaded Data has not 3 columns",closeOnClickOutside = TRUE, type = "error")
                    returnValue()
                }
            }
        }
        return(listFiles)
    })


    observeEvent(input$load_omics_net_button, {
        shinyalert::shinyalert(
            title = "Wait",
            text = "Waiting for data loading",
            size = "xs",
            closeOnEsc = TRUE,
            closeOnClickOutside = TRUE,
            html = TRUE,
            type = "info",
            showConfirmButton = TRUE,
            confirmButtonText = "OK",
            confirmButtonCol = "#004192",
            showCancelButton = FALSE,
            imageUrl = "",
            animation = TRUE
        )

    },
    ignoreNULL = FALSE,
    ignoreInit = TRUE
    )


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

    annotation1 <- reactiveVal(NULL)
    shiny::observeEvent(input$load_anno_button1, {
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
                shinyalert::shinyalert("Type Error", "Uploaded Data is not a dataframe",closeOnClickOutside = TRUE, type = "error")
                returnValue()
            }
        }
        annotation1(listFiles)
    })

    observeEvent(input$load_anno_button1, {
        shinyalert::shinyalert(
            title = "Wait",
            text = "Waiting for data loading",
            size = "xs",
            closeOnEsc = TRUE,
            closeOnClickOutside = TRUE,
            html = TRUE,
            type = "info",
            showConfirmButton = TRUE,
            confirmButtonText = "OK",
            confirmButtonCol = "#004192",
            showCancelButton = FALSE,
            imageUrl = "",
            animation = TRUE
        )

    },
    ignoreNULL = FALSE,
    ignoreInit = TRUE
    )


    ############################################################################################
    #                                RWR: MULTIPLEX NETWORK (INPUT)                            #
    ############################################################################################

    shinyjs::addCssClass(selector = "a[data-value='rwr_tab_2']", class = "inactiveLink")
    observe({
        if (length(shiny::reactiveValuesToList(filteredOmicsNetworks)) == length(omics_files())){
            shinyjs::removeCssClass(selector = "a[data-value='rwr_tab_2']", class = "inactiveLink")
        }else{
            return(NULL)
        }
    })
    observe({
        if (length(loaded_omics_net_list()) != 0 ){
            shinyjs::removeCssClass(selector = "a[data-value='rwr_tab_2']", class = "inactiveLink")
        }else{
            return(NULL)
        }
    })

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

    multiplex_network <- shiny::eventReactive(input$gen_multiplex_btn, {
        print('Generating the multiplex network...')
        net_list <- if (!is.null(g_net_list())) g_net_list() else l_net_list()

        if(input$weightMultiplex == 'YES'){
            multiplex <- MoNETA::create_multiplex(net_list, weighted = TRUE)
            multiplex_type <- 'Weighted'
            if(input$pruneMultiplex == 'YES'){
                if (!is.numeric(input$pruneMultiplex_th)) {
                    shinyalert::shinyalert(title = 'Error', text = 'The threshold must be a real number.',closeOnClickOutside = TRUE,type = 'error')
                    return(NULL)
                }
                multiplex <- MoNETA::prune_multiplex_network(multiplex, input$pruneMultiplex_th)
            }else{
                multiplex <- multiplex
            }
        }else{
            multiplex <- MoNETA::create_multiplex(net_list, weighted = FALSE)
            multiplex_type <- 'Unweighted'
        }

        shinyalert::shinyalert(title = 'Success', type = 'success',closeOnClickOutside = TRUE,
                               text = paste(multiplex_type,' Multiplex Network created. <br/> Now press <b> "Next" </b> in the top left-hand corner to continue'), html = TRUE)
        return(multiplex)
    })


    ############################################################################################
    #                                   RWR: PARAMETERS (INPUT)                                #
    ############################################################################################

    shinyjs::addCssClass(selector = "a[data-value='rwr_tab_3']", class = "inactiveLink")
    observe({
        if (isTruthy(input$gen_multiplex_btn)){
            shinyjs::removeCssClass(selector = "a[data-value='rwr_tab_3']", class = "inactiveLink")
        }else{
            return(NULL)
        }
    })

    ######################## RWR-SIM ###########################

    RWR_output <- shiny::reactiveValues()
    shiny::observeEvent(input$rwr_button, {
        if (is.null(multiplex_network())){
            return(NULL)
        } else{
            print(omics_tau_list())
            print(unlist(omics_tau_list()))
            print(sum(unlist(omics_tau_list())))
            if (!is.na(omics_tau_list())) {
                if (sum(unlist(omics_tau_list())) != 1) {
                    shinyalert::shinyalert("Error", "The sum of restarting probabilities per layer (taus) must be equal to 1",closeOnClickOutside = TRUE, type = "error")
                    return(NULL)
                }
            }
            print( 'skip')
            max_cores <- parallel::detectCores()
            if (!is.numeric(input$cores_rwr) | input$cores_rwr < 1 | input$cores_rwr > max_cores){
                shinyalert::shinyalert(title = 'Error', text = paste('Select a number of cores between 1 and', max_cores, '.'),closeOnClickOutside = TRUE, type = 'error')
                return(NULL)
            }

            req(multiplex_network)
            omics_multiplex <- multiplex_network()

            if (input$bioInf_transition == 'YES') {
                net_list <- if (!is.null(g_net_list())) g_net_list() else l_net_list()
                trans_mat <- MoNETA::create_layer_transition_matrix(net_list)
            }else {
                trans_mat <- input$editable_trans_mat
            }


            shinyalert::shinyalert(
                title = "Wait",
                text = "Waiting for RWR similarity matrix generation. <br/> Please check the <b> progressbar </b>  in the bottom right-hand corner.",
                closeOnEsc = TRUE,
                closeOnClickOutside = TRUE,
                html = TRUE,
                type = "info",
                showConfirmButton = TRUE,
                confirmButtonText = "OK",
                confirmButtonCol = "#004192",
                showCancelButton = FALSE,
                imageUrl = "",
                animation = TRUE
            )

            progress <- shiny::Progress$new()
            on.exit(progress$close())

            progress$set(message = "MoNETA", detail = paste("Doing RWR"), value = 0)

            rwr_simMat <- MoNETA::gen_sim_mat_M(network = omics_multiplex,
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
                    paste('- Weighted jump_neigh edges:            ', input$omics_weight_jump_neigh),''
                ),
                'Transition layers matrix' = trans_mat
            )

            RWR_output$rwr_simMat <- rwr_simMat
            RWR_output$Summary <- Summary

            shinyalert::shinyalert(title = 'Success', type = 'success',closeOnClickOutside = TRUE,
                                   text = 'Random Walk with Restart done! <br/> Now press <b> "Next" </b> in the top left-hand corner to continue', html = TRUE)
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
                    sink(file = glue::glue(temp_directory, "/{obj}.txt"), append =  TRUE)
                    print(to_download[[obj]])
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

    shinyjs::addCssClass(selector = "a[data-value='dr_tab']", class = "inactiveLink")
    observe({
        if (length(shiny::reactiveValuesToList(RWR_output)) != 0){
            shinyjs::removeCssClass(selector = "a[data-value='dr_tab']", class = "inactiveLink")
        }else{
            return(NULL)
        }
    })


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

            shinyalert::shinyalert(
                title = "Wait",
                text = "It might take a little time at the end. <br/> Please check the <b> progressbar </b> in the bottom right-hand corner.",
                closeOnEsc = TRUE,
                closeOnClickOutside = TRUE,
                html = TRUE,
                type = "info",
                showConfirmButton = TRUE,
                confirmButtonText = "OK",
                confirmButtonCol = "#004192",
                showCancelButton = FALSE,
                imageUrl = "",
                animation = TRUE
            )

            progress <- shiny::Progress$new()
            on.exit(progress$close())
            progress$set(message = "MoNETA", detail = paste("Doing Embedding"), value = 0)

            if (input$emb_opt == 'YES'){
                emb_mat <- MoNETA::get_embedding(matrix = rwr_mat, embedding_size = input$dim_emb, cores = input$cores_emb)
                dr_output$emb_mat <- emb_mat
            }else{
                emb_mat <- rwr_mat
            }

            if (input$dr_method == 'UMAP'){
                progress$inc(0.8, detail = paste("Doing UMAP"))
                dr_mat <- MoNETA::get_parallel_umap_embedding(matrix = emb_mat,
                                                              embedding_size = 2,
                                                              n_threads = input$threads)
            }else if (input$dr_method == 'PCA'){
                progress$inc(0.8, detail = paste("Doing PCA"))
                dr_mat <- MoNETA::get_pca_embedding(matrix = emb_mat, embedding_size = input$emb_size)
            }else if (input$dr_method == 'tSNE'){
                progress$inc(0.8, detail = paste("Doing tSNE"))
                dr_mat <- get_tsne_embedding(matrix = emb_mat,
                                             embedding_size = input$emb_size_tsne,
                                             perplexity = input$perplexity, max_iter = input$max_iter,
                                             num_threads = input$threads)
            }

            progress$inc(1, detail = paste("Done!"))
            dr_output$dr_mat <- dr_mat

            shinyalert::shinyalert(title = 'Success', type = 'success', closeOnClickOutside = TRUE,
                                   text = 'Dimensionality Reduction step done! <br/> Now press <b> "Next" </b> in the top left-hand corner to continue', html = TRUE)
        }
    })


    dr_plot <- shiny::reactive({
        method <- ifelse(input$dr_method == 'UMAP', 'UMAP of',
                         ifelse(input$dr_method == 'PCA', 'PCA of', 'tSNE of'))
        emb <- ifelse(input$emb_opt == 'YES', 'the Embedded RWR matrix', 'the RWR matrix')
        dr_mat <- shiny::reactiveValuesToList(dr_output)$dr_mat
        dr_plot <- MoNETA::plot_2D_matrix(coord = dr_mat[1:2,], nodes_anno = data.frame(id = colnames(dr_mat)),
                                          id_name = 'id', interactive = FALSE, wo_legend = FALSE) +
            ggplot2::ggtitle(paste(method, emb)) +
            ggplot2::xlab(paste0(input$dr_method, '1')) +  ggplot2::ylab( paste0(input$dr_method, '2')) +
            ggplot2::theme(text = ggplot2::element_text(family="Times New Roman"),
                           plot.title = ggplot2::element_text(size = 16, hjust = 0.5, face = 'bold'),
                           axis.title.x = ggplot2::element_text(size = 10),
                           axis.title.y = ggplot2::element_text(size = 10))
        return(dr_plot)
    })


    ############################################################################################
    #                                         CLUSTERING (INPUT)                               #
    ############################################################################################

    shinyjs::addCssClass(selector = "a[data-value='cl_tab']", class = "inactiveLink")
    observe({
        if (length(shiny::reactiveValuesToList(dr_output)) != 0){
            shinyjs::removeCssClass(selector = "a[data-value='cl_tab']", class = "inactiveLink")
        }else{
            return(NULL)
        }
    })

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
            gPlot_mat <- dplyr::tibble("id" = rownames(t_mat))

            shinyalert::shinyalert(
                title = "Wait",
                text = paste0("It might take a little time at the end of the clustering analysis. <br/>",
                              "Please check the <b> progressbar </b>  in the bottom right-hand corner."
                ),
                closeOnEsc = TRUE,
                closeOnClickOutside = TRUE,
                html = TRUE,
                type = "info",
                showConfirmButton = TRUE,
                confirmButtonText = "OK",
                confirmButtonCol = "#004192",
                showCancelButton = FALSE,
                imageUrl = "",
                animation = TRUE
            )

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

                if ((  !is.null(anno()) |  shiny::isTruthy(input$anno_file1)) && ifelse(input$shape_bool == 'YES', TRUE, FALSE)){
                    anno <- if ( length(annotation()) != 0) anno() else anno1()
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
                cluster$num_clust <- length(unique(cluster$clust))
                progress$inc(1, detail = paste("Done!"))

                if ((  !is.null(anno())  |  shiny::isTruthy(input$anno_file1)) && ifelse(input$shape_bool == 'YES', TRUE, FALSE)){
                    anno <- if ( length(annotation()) != 0) anno() else anno1()
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
                    if ((  !is.null(anno())  |  shiny::isTruthy(input$anno_file1)) && ifelse(input$shape_bool == 'YES', TRUE, FALSE)){
                        anno <- if ( !is.null(anno()) ) anno() else anno1()
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

            shinyalert::shinyalert(
                title = "Success",
                text = paste0("Clustering Analysis done! <br/> You have completed MoNETA pipeline!"
                ),
                closeOnEsc = TRUE,
                closeOnClickOutside = TRUE,
                html = TRUE,
                type = "success",
                showConfirmButton = TRUE,
                confirmButtonText = "OK",
                confirmButtonCol = "#004192",
                showCancelButton = FALSE,
                imageUrl = "",
                animation = TRUE
            )
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
            return(paste0('DBSCAN clustering (epsilon: ', input$eps, ', #cluster: ', shiny::reactiveValuesToList(clusters)$num_clust, ')' ))
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
                ggplot2::ggtitle(title()) +
                ggplot2::xlab(paste0(input$dr_method, '1')) +  ggplot2::ylab( paste0(input$dr_method, '2')) +
                ggplot2::theme(text = ggplot2::element_text(family="Times New Roman"),
                               plot.title = ggplot2::element_text(size = 16, hjust = 0.5, face = 'bold'),
                               axis.title.x = ggplot2::element_text(size = 10),
                               axis.title.y = ggplot2::element_text(size = 10))
        } else{
            cl_plot <- MoNETA::plot_2D_matrix(coord = clusters$dr_mat[1:2,], nodes_anno = clu_anno,
                                              id_name = 'id', id_anno_color = 'clust', id_anno_shape = 'shape',
                                              interactive = FALSE, wo_legend = FALSE) +
                ggplot2::ggtitle(title()) +
                ggplot2::xlab(paste0(input$dr_method, '1')) +  ggplot2::ylab( paste0(input$dr_method, '2')) +
                ggplot2::theme(text = ggplot2::element_text(family="Times New Roman"),
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
            shinydashboard::box(width = 12,
                                shiny::HTML(" "),
                                shiny::HTML(paste0(
                                    "<span style='font-family: LM Roman 10;' > <h4>",
                                    paste("<b>", names(omics_files()), "</b> </span> "),
                                    "<span style='font-family: LM Roman 10;' > loaded: ", paste("<b>", nfeatures, "</b>"), " features and", paste("<b>", nsamples, "</b>"), " samples <br/>"),
                                    '</span> </h4>')
            )
        } else {
            shiny::HTML("<br/>")
        }
    })

    output$anno_info <- shiny::renderUI({
        input$load_anno_button
        if ( length(annotation()) != 0 )
        {
            ncolumns <- ncol(annotation()$annotation)
            nrows <- nrow(annotation()$annotation)
            shinydashboard::box(width = 12,
                                shiny::HTML("<span style='font-family: LM Roman 10;' >"),
                                shiny::HTML(paste0(paste("<h4> <b>", 'Annotation file', "</b>"),
                                                   "<span style='font-family: LM Roman 10;' > loaded: ",
                                                   paste("<b>", ncolumns, "</b>"), " columns and ",
                                                   paste("<b>", nrows, "</b>"), " rows <br/> </span> ")),
                                shiny::HTML(paste0(paste("<b> Column names: </b>"), "<br/>")),
                                shiny::HTML(paste("&ensp; <span style='font-family: LM Roman 10;' >", colnames(annotation()$annotation), "<br/>")),
                                shiny::HTML('</h4> </span> ')
            )
        } else {
            shiny::HTML("<br/>")
        }
    })


    output$anno_info_extra <- shiny::renderUI({
        #input$load_anno_button
        if ( length(annotation()) != 0 & !is.null(omics_files())){
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
                shinydashboard::box(width = 12,
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
    #processed_mat_list <- list(NULL)
    proc_matrices <- shiny::reactiveValues()
    shiny::observeEvent(
        eventExpr = {
            buttons <- paste0("process_mat_", names(omics_files()))
            list_of_buttons = NULL
            for(var in buttons) {
                list_of_buttons <- append(list_of_buttons, input[[var]])
            }
            req(list_of_buttons)
        },
        handlerExpr = {
            for (x in names(omics_files())){
                processed_omics_mat <- omics_files()[[x]]
                if (input[[paste0("remove0col_", x)]] == 'YES')
                    processed_omics_mat <- MoNETA::remove_zeros_cols(processed_omics_mat)
                if (input[[paste0("norm_", x)]] == 'YES')
                    processed_omics_mat <- MoNETA::normalize_omics(processed_omics_mat)
                proc_matrices[[paste0('p_', x)]] <- processed_omics_mat
            }
            message('All matrices processed!')

            if (length(omics_files()) == 1){
                shinyalert::shinyalert(
                    title = "Success",
                    text = 'Processig step done! <br/> Now press <b> "Next" </b> in the top right-hand corner to continue.',
                    closeOnEsc = TRUE,
                    closeOnClickOutside = TRUE,
                    html = TRUE,
                    type = "success",
                    showConfirmButton = TRUE,
                    confirmButtonText = "OK",
                    confirmButtonCol = "#004192",
                    showCancelButton = FALSE,
                    imageUrl = "",
                    animation = TRUE
                )
            }else if (length(omics_files()) > 1){
                shinyalert::shinyalert(
                    title = "Success",
                    text = 'Processig step done! <br/>  Now move on to the <b> Intersection </b> step.',
                    closeOnEsc = TRUE,
                    closeOnClickOutside = TRUE,
                    html = TRUE,
                    type = "success",
                    showConfirmButton = TRUE,
                    confirmButtonText = "OK",
                    confirmButtonCol = "#004192",
                    showCancelButton = FALSE,
                    imageUrl = "",
                    animation = TRUE
                )
            }
        },
        ignoreInit = T
    )

    output$process_omics_mat <- shiny::renderUI({
        if (is.null(omics_files()))
            return(NULL)
        else{
            thetabs <- lapply(names(omics_files()), function(x) {

                tab <- shiny::tabPanel(title = x,
                                       shiny::fluidRow(
                                           shiny::column(width = 12,
                                                         shiny::radioButtons(inputId = paste0("remove0col_", x),
                                                                             label = shiny::h4(shiny::span("Do you want to remove zeros columns?", style = "font-weight: bold")),
                                                                             choices = c('YES', 'NO'), selected = 'NO'
                                                         ),
                                                         shiny::radioButtons(inputId = paste0("norm_", x),
                                                                             label =  shiny::h4(shiny::span("Do you want to normalize the omics matrix by column?", style = "font-weight: bold")),
                                                                             choices = c('YES', 'NO'), selected = 'NO'
                                                         )
                                           )
                                       ), shiny::hr(),
                                       shiny::fluidRow(
                                           shiny::column(width = 12, shiny::actionButton(paste0('process_mat_', x), "Submit"))
                                       )
                )
                return(tab)
            })
            thetabs$width <- 12
            do.call(shinydashboard::tabBox, thetabs)
        }
    })

    ############## Processed matrices info #############

    output$pro_matrices_sum_info <- shiny::renderUI({
        processed_mats <- shiny::reactiveValuesToList(proc_matrices)
        shiny::isolate({
            if (length(processed_mats) != 0){
                nsamples <- sapply(processed_mats, ncol)
                nfeatures <- sapply(processed_mats, nrow)
                shinydashboard::box(width = 12,
                                    shiny::HTML(paste0(paste("<b>", names(processed_mats), "</b>"), " processed: ",
                                                       paste("<b>", nfeatures, "</b>"), " features and ", paste("<b>", nsamples, "</b>"), " samples <br/>")))
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
            shinydashboard::box(width = 12, solidHeader = TRUE, collapsible = FALSE,
                                title = shiny::h3(shiny::span('Download', style = "font-weight: bold")),
                                checkboxGroupInput('download_proc_mat',
                                                   label = shiny::h4(shiny::span('Select one or more processed matrices', style = "font-weight: bold")),
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
                file_name <- glue::glue("{obj}.csv")
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
            shinydashboard::box(width = 12,  status = "primary", solidHeader = FALSE,
                                shiny::h2(shiny::span("Omics matrices Intersection", style = "font-weight: bold")),
                                shiny::radioButtons(inputId = 'intersect_opt',
                                                    label = shiny::h4(shiny::span('Do you want to intersect omics matrices by samples?', style = "font-weight: bold")),
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
                shinydashboard::box(width = 12,
                                    shiny::HTML(paste0(paste("<b>", names(int_mats), "</b>"), " processed: ",
                                                       paste("<b>", nfeatures, "</b>"), " features and ", paste("<b>", nsamples, "</b>"), " samples <br/>")))
            } else {
                return(NULL)
            }
        })
    })

    ###
    output$back2P1 <- renderUI({
        if (length(omics_files()) >= 1 ){
            actionButton('back2P1', label = 'Back', icon("paper-plane"),
                         style="color: #fff; background-color: #337ab7; border-color: #2e6da4; padding:20px")
        }else{
            return(NULL)
        }
    })

    observeEvent(input$back2P1, {
        shinydashboard::updateTabItems(session, inputId = "tabs",
                                       selected = "mat_sub_1")
    })

    ###
    output$jump2P3 <- renderUI({
        if (length(omics_files()) == 1 && length(reactiveValuesToList(proc_matrices)) != 0){
            actionButton('jump2P3', label = 'Next', icon("paper-plane"),
                         style="color: #fff; background-color: #337ab7; border-color: #2e6da4; padding:20px")
        } else if (length(omics_files()) >= 1 && length(reactiveValuesToList(proc_matrices)) != 0 && isTruthy(input$intersect_btn)){
            actionButton('jump2P3', label = 'Next', icon("paper-plane"),
                         style="color: #fff; background-color: #337ab7; border-color: #2e6da4; padding:20px")
        }else{
            return(NULL)
        }
    })

    observeEvent(input$jump2P3, {
        shinydashboard::updateTabItems(session, inputId = "tabs",
                                       selected = "net_sub_1")
    })


    ############################################################################################
    #                        NETWORK : NETWORK INFERENCE (OUTPUT)                      #
    ############################################################################################

    ############## output parameters for OMICS network inference  #############

    omics_network_list <- list(NULL)
    gene_networks <- shiny::reactiveValues()

    count_net <- shiny::reactiveValues()
    observe({
        if (length(shiny::reactiveValuesToList(intersected_omics_mat)$matrices) == 0){
            if (length(omics_files()) == 1){
                input_names <- names(shiny::reactiveValuesToList(proc_matrices))
            }else{
                return(NULL)
            }
        }else {
            input_names <- names(shiny::reactiveValuesToList(intersected_omics_mat)$matrices)
        }

        for (x in input_names){
            count_net[[x]] <- 0
        }
    })

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

                max_cores <- parallel::detectCores()

                if (!is.numeric(input[[paste0('sparsity_', x)]]) | input[[paste0('sparsity_', x)]] < 0){
                    shinyalert::shinyalert(title = 'Error', closeOnClickOutside = TRUE, text = 'Sparsity must be a positive real number.',type = 'error')
                    return(NULL)
                }
                if (!is.integer(input[[paste0('knn_', x)]]) | input[[paste0('knn_', x)]] <= 0){
                    shinyalert::shinyalert(title = 'Error', closeOnClickOutside = TRUE, text = 'Knn must be a positive integer.', type = 'error')
                    return(NULL)
                }
                if (!is.numeric(input[[paste0('cores_', x)]]) | input[[paste0('cores_', x)]] < 1 | input[[paste0('cores_', x)]] > max_cores){
                    shinyalert::shinyalert(title = 'Error', closeOnClickOutside = TRUE, text = paste('Select a number of cores between 1 and', max_cores, '.'), type = 'error')
                    return(NULL)
                }
                if (input[[paste0('cores_', x)]] == max_cores){

                    shinyalert::shinyalert(title = 'Error', closeOnClickOutside = TRUE,
                                           text = tagList('You have selected the maximum number of cores available.',
                                                          #shiny::radioButtons(inputId = 'warn_continue',
                                                          #             label = 'Do you want to proceed anyway?',
                                                          #             choices = c('Yes', 'No'), selected = 'No')
                                           ),
                                           type = 'error', html = TRUE
                    )


                    #shinyalert::shinyalert(title = 'Warning',
                    #                       text = tagList('You have selected the maximum number of cores available.',
                    #                                      shiny::radioButtons(inputId = 'warn_continue',
                    #                                                   label = 'Do you want to proceed anyway?',
                    #                                                   choices = c('Yes', 'No'), selected = 'No')
                    #                                      #shiny::textInput(inputId = 'warn_continue',
                    #                                      #                 label = 'Do you want to proceed anyway? (type "yes" or "no")')
                    #                                      ),
                    #                       type = 'warning', html = TRUE
                    #                       )
                    #shinyalert(html = TRUE, text = tagList(
                    #    textInput("warn_continue", "Do you want to proceed anyway? (y/n)"),
                    #    )
                    #)
                    #print(input$shinyalert)
                    #print(input$warn_continue)
                    return(NULL)
                }

                shinyalert::shinyalert(
                    title = "Wait",
                    text = paste0("Waiting for the inference of", paste("<b>", x, "</b>"), "network. <br/>",
                                  "Please check the <b> progressbar </b> in the bottom right-hand corner."
                    ),
                    closeOnEsc = TRUE,
                    closeOnClickOutside = TRUE,
                    html = TRUE,
                    type = "info",
                    showConfirmButton = TRUE,
                    confirmButtonText = "OK",
                    confirmButtonCol = "#004192",
                    showCancelButton = FALSE,
                    imageUrl = "",
                    animation = TRUE
                )

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

                count_net[[x]] <- 1
                print(shiny::reactiveValuesToList(count_net))
                net_created <- sum(unlist(shiny::reactiveValuesToList(count_net)))


                if (net_created == length(omics_files())){
                    shinyalert::shinyalert(
                        title = "Success",
                        text = paste0(paste("<b>", x, "</b>"), " network created.<br/>",
                                      paste("<b> Network inferring status: </b>"), net_created, '/', length(omics_files()),"<br/>",
                                      'Now press <b> "Next" </b> in the top right-hand corner to continue.'
                        ),
                        closeOnEsc = TRUE,
                        closeOnClickOutside = TRUE,
                        html = TRUE,
                        type = "success",
                        showConfirmButton = TRUE,
                        confirmButtonText = "OK",
                        confirmButtonCol = "#004192",
                        showCancelButton = FALSE,
                        imageUrl = "",
                        animation = TRUE
                    )
                } else {
                    shinyalert::shinyalert(
                        title = "Success",
                        text = paste0(paste("<b>", x, "</b>"), " network created. <br/>",
                                      paste("<b> Network inferring status: </b>"), net_created, '/', length(omics_files())
                        ),
                        closeOnEsc = TRUE,
                        closeOnClickOutside = TRUE,
                        html = TRUE,
                        type = "success",
                        showConfirmButton = TRUE,
                        confirmButtonText = "OK",
                        confirmButtonCol = "#004192",
                        showCancelButton = FALSE,
                        imageUrl = "",
                        animation = TRUE
                    )
                }

                output[[paste0(x, "_spinner")]] <- renderUI({
                    shinycssloaders::withSpinner(visNetwork::visNetworkOutput(paste0(x, 'plot_net')))
                })

                output[[paste0(x, 'plot_net')]] <- visNetwork::renderVisNetwork({
                    ids <- data.frame(id = base::union(net$source, net$dest))
                    net_plot <- MoNETA::plot_net(edgeList = net, nodes_anno = ids, id_name = 'id', id_anno_color = 'id')
                    return(net_plot)
                })

                return(net)
            })


            omics_network_list[[paste0('update_omics_net_plot_', x)]] <-
                shiny::observeEvent(input[[paste0('update_omics_net_plot_', x)]], {

                    output[[paste0(x, "_spinner")]] <- renderUI({
                        shinycssloaders::withSpinner(visNetwork::visNetworkOutput(paste0(x, 'plot_net')))
                    })

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
                            final_anno <- f_anno %>% dplyr::as_tibble(.) %>% dplyr::select(id, everything())

                            net_plot <- MoNETA::plot_net(edgeList = net, nodes_anno = final_anno,
                                                         id_name = id, id_anno_color = color, id_anno_shape = shape,
                                                         interactive = TRUE, wo_legend = FALSE)
                            return(net_plot)
                        })
                    })
                })

            output[[paste0(x, '_custom_net')]] <- shiny::renderUI({
                req(input[[paste0('generate_omics_net_', x)]])
                shiny::isolate({
                    if (length(annotation()) != 0){
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
                    if ( length(annotation()) != 0){
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
                                                         label = shiny::h4(shiny::span("Select a distance measure", style = "font-weight: bold")),
                                                         choices =  c("Euclidean", "Manhattan", "Cosine"),
                                                         selected = "Euclidean"
                                                     ),
                                                     shiny::numericInput(
                                                         inputId = paste0('sparsity_', x),
                                                         label = shiny::h4(shiny::span("Select the sparsity", style = "font-weight: bold")),
                                                         min = 0, max = 1, value = 0.7, step = 0.1
                                                     ),
                                                     shiny::numericInput(
                                                         inputId = paste0('knn_', x),
                                                         label = shiny::h4(shiny::span("Select the number of neighbors (knn)", style = "font-weight: bold")),
                                                         min = 1, max = 100, value = 25
                                                     ),
                                                     shiny::numericInput(
                                                         inputId = paste0('cores_', x),
                                                         label = shiny::h4(shiny::span("Select the number of cores", style = "font-weight: bold")),
                                                         min = 1, max = 200, value = 1
                                                     ),
                                                     shiny::radioButtons(inputId = paste0("kstar_", x),
                                                                         shiny::h4(shiny::span("Do you want to perform the kstar?", style = "font-weight: bold")),
                                                                         choices = c('YES', 'NO'), selected = 'YES'
                                                     )
                                       ),
                                       shiny::column(width = 8,
                                                     #shinycssloaders::withSpinner(visNetwork::visNetworkOutput(paste0(x, 'plot_net'))),
                                                     shiny::uiOutput(paste0(x, "_spinner")),
                                                     shiny::uiOutput(paste0(x, '_custom_net'))
                                       )
                                   ), shiny::hr(),
                                   shiny::fluidRow(
                                       shiny::column(width = 6, shiny::actionButton(paste0('generate_omics_net_', x), paste0("Generate network ", x))),
                                       shiny::column(width = 6, shiny::uiOutput(paste0('update_omics_net_plot_', x)), align = 'right')
                                   )
            )
            return(tab)
        })
        thetabs$width <- 12
        do.call(shinydashboard::tabBox, thetabs)
    })

    ############## download networks  #############

    output$download_net_box <- shiny::renderUI({
        if (length(shiny::reactiveValuesToList(gene_networks)) == 0) {
            return(NULL)
        } else {
            shinydashboard::box(shiny::h3(shiny::span('Download', style = "font-weight: bold")), solidHeader = TRUE, collapsible = FALSE, width = 12,
                                checkboxGroupInput('download_net',
                                                   label = shiny::h4(shiny::span('Select one or more networks', style = "font-weight: bold")),
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
                file_name <- glue::glue("{obj}.csv")
                readr::write_csv(as.data.frame(to_download[[obj]]), file.path(temp_directory, file_name))
            }

            zip::zip(
                zipfile = file,
                files = dir(temp_directory),
                root = temp_directory
            )
        }
    )


    ###
    output$back2P2 <- renderUI({
        if (length(omics_files()) == 1 && length(reactiveValuesToList(proc_matrices)) != 0){
            actionButton('back2P2', label = 'Back', icon("paper-plane"),
                         style="color: #fff; background-color: #337ab7; border-color: #2e6da4; padding:20px")
        } else if (length(omics_files()) >= 1 && length(reactiveValuesToList(proc_matrices)) != 0 && isTruthy(input$intersect_btn)){
            actionButton('back2P2', label = 'Back', icon("paper-plane"),
                         style="color: #fff; background-color: #337ab7; border-color: #2e6da4; padding:20px")
        }else{
            return(NULL)
        }
    })

    observeEvent(input$back2P2, {
        shinydashboard::updateTabItems(session, inputId = "tabs",
                                       selected = "mat_sub_2")
    })

    ###
    output$jump2P4 <- renderUI({
        if (length(shiny::reactiveValuesToList(gene_networks)) == length(omics_files())){
            actionButton('jump2P4', label = 'Next', icon("paper-plane"),
                         style="color: #fff; background-color: #337ab7; border-color: #2e6da4; padding:20px")
        }else{
            return(NULL)
        }
    })

    observeEvent(input$jump2P4, {
        shinydashboard::updateTabItems(session, inputId = "tabs",
                                       selected = "net_sub_2")
    })


    ############################################################################################
    #                               NETWORK : NETWORK FILTERING (OUTPUT)                       #
    ############################################################################################

    ############## Filtering networks box #############

    filteredOmicsNetworks_list <- list(NULL)
    filteredOmicsNetworks <- shiny::reactiveValues()

    count_unet <- shiny::reactiveValues()
    observe({
        input_names <- names(shiny::reactiveValuesToList(gene_networks))
        for (x in input_names){
            count_unet[[x]] <- 0
        }
    })


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
                        summary <- netSummary(dplyr::tibble(fnet))
                        res <- c()
                        for (i in 2:length(summary))
                            res[i-1] <- (paste0(paste('<b>', names(summary)[i], '</b>'),
                                                ': <span style="font-family: LM Roman 10;"> ', toString(summary[[i]]), '</span> <br/>'))
                        shiny::HTML(paste('<h4>', res, '</h4> '))
                    })
                })

            filteredOmicsNetworks_list[[paste0('button_', x)]] <- shiny::observeEvent(input[[paste0('button_', x)]], {
                net <- nets[[x]]
                fnet <- net[net$weight <= input[[paste0('omics_range', x)]],]

                output[[paste0(x, "_spinner_weightDist")]] <- renderUI({
                    shinycssloaders::withSpinner(plotly::plotlyOutput(paste0(x, "_weightDist")))
                })

                output[[paste0(x, "_weightDist")]] <- plotly::renderPlotly({
                    #shiny::isolate({
                    plotly::ggplotly(ggplot2::ggplot(fnet, ggplot2::aes(x= weight)) + ggplot2::geom_histogram() +
                                         ggplot2::ggtitle('Weight Distribution') +
                                         ggplot2::xlab("Weight") + ggplot2::ylab("Counts") +
                                         ggplot2::theme(text = ggplot2::element_text(family="LM Roman 10"),
                                                        plot.title = ggplot2::element_text(size = 16, hjust = 0.5, face = 'bold'),
                                                        axis.title.x = ggplot2::element_text(size = 10),
                                                        axis.title.y = ggplot2::element_text(size = 10))
                    )
                    #})
                })


                output[[paste0(x, "_spinner_nodeDegreeDist")]] <- renderUI({
                    shinycssloaders::withSpinner(plotly::plotlyOutput(paste0(x, "_nodeDegreeDist")))
                })

                output[[paste0(x, "_nodeDegreeDist")]] <- plotly::renderPlotly({
                    #shiny::isolate({
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
                    #})
                })


                output[[paste0(x, "_spinner")]] <- renderUI({
                    shinycssloaders::withSpinner(visNetwork::visNetworkOutput(paste0(x, "_plot_fnet")))
                })

                output[[paste0(x, '_plot_fnet')]] <- visNetwork::renderVisNetwork({
                    shiny::isolate({
                        ids <- data.frame(id = base::union(fnet$source, fnet$dest))
                        net_plot <- MoNETA::plot_net(edgeList = fnet, nodes_anno = ids, id_name = 'id', id_anno_color = 'id')
                        return(net_plot)
                    })
                })

                output[[paste0(x, '_custom_fnet')]] <- shiny::renderUI({
                    #shiny::isolate({
                    if (length(annotation()) != 0){
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


            #filteredOmicsNetworks_list[[paste0("update_", x)]] <- shiny::observeEvent(input[[paste0("update_", x)]], {
            filteredOmicsNetworks_list$update_all_nets <- shiny::observeEvent(input$update_all_nets, {
                #shiny::isolate({
                filt_net <- net[net$weight <= input[[paste0('omics_range', x)]],]
                filteredOmicsNetworks[[paste0('f_', x)]] <- filt_net
                print('Filering network...done!')

                count_unet[[x]] <- 1
                print(shiny::reactiveValuesToList(count_unet))
                net_update <- sum(unlist(shiny::reactiveValuesToList(count_unet)))


                if (net_update == length(omics_files())){
                    shinyalert::shinyalert(
                        title = "Success",
                        text = paste0(paste("<b>", x, "</b>"), " network updated. <br/>",
                                      paste("<b> Network update status: </b>"), net_update, '/', length(omics_files()), "<br/>",
                                      'Now press <b> "Next" </b> in the top right-hand corner to continue.'
                        ),
                        closeOnEsc = TRUE,
                        closeOnClickOutside = TRUE,
                        html = TRUE,
                        type = "success",
                        showConfirmButton = TRUE,
                        confirmButtonText = "OK",
                        confirmButtonCol = "#004192",
                        showCancelButton = FALSE,
                        imageUrl = "",
                        animation = TRUE
                    )
                } else{
                    shinyalert::shinyalert(
                        title = "Success",
                        text = paste0(paste("<b>", x, "</b>"), " network updated. <br/>",
                                      paste("<b> Network update status: </b>"), net_update, '/', length(omics_files())
                        ),
                        closeOnEsc = TRUE,
                        closeOnClickOutside = TRUE,
                        html = TRUE,
                        type = "success",
                        showConfirmButton = TRUE,
                        confirmButtonText = "OK",
                        confirmButtonCol = "#004192",
                        showCancelButton = FALSE,
                        imageUrl = "",
                        animation = TRUE
                    )
                }
                return(filt_net)
                #})
            })

            filteredOmicsNetworks_list[[paste0('update_f_net_plot_', x)]] <-
                shiny::observeEvent(input[[paste0('update_f_net_plot_', x)]], {

                    output[[paste0(x, "_spinner")]] <- renderUI({
                        shinycssloaders::withSpinner(visNetwork::visNetworkOutput(paste0(x, "_plot_fnet")))
                    })

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
                            final_anno <- f_anno %>% dplyr::as_tibble(.) %>%  dplyr::select(id, everything())

                            net_plot <- MoNETA::plot_net(edgeList = filt_net, nodes_anno = final_anno,
                                                         id_name = id, id_anno_color = color, id_anno_shape = shape,
                                                         interactive = TRUE, wo_legend = FALSE)
                            return(net_plot)
                        })
                    })
                })


            ### Update plot button
            output[[paste0('update_f_net_plot_', x)]] <- shiny::renderUI({
                #shiny::isolate({
                if ( length(annotation()) != 0 &  shiny::isTruthy(input[[paste0('button_', x)]])){
                    shiny::actionButton(inputId = paste0('update_f_net_plot_', x), label = 'Update plot')
                }else{
                    return(NULL)
                }
                #})
            })

            net <- nets[[x]]
            min <-  round(min(net$weight), digits = 2)
            max <-  round(max(net$weight), digits = 2)

            tab <- shiny::tabPanel(x,
                                   shiny::fluidRow(
                                       shiny::column(width = 8,
                                                     shiny::sliderInput(inputId = paste0("omics_range", x),
                                                                        label = shiny::h4(shiny::span( 'Select a threshold:', style = "font-weight: bold")),
                                                                        min = min, max = max, value = max, step = 0.01)
                                       ),
                                       shiny::column(width = 4,
                                                     fluidRow(
                                                         shiny::tags$hr(style='border-top: 1px solid white;'),
                                                         shiny::tags$hr(style='border-top: 1px solid white;'),
                                                         shiny::actionButton(inputId = paste0('button_', x), 'Show'), align = 'center')
                                       )
                                   ),

                                   shiny::fluidRow(
                                       shinydashboard::tabBox( width = 12,
                                                               shiny::tabPanel(title = 'Distributions',
                                                                               shiny::fluidRow(
                                                                                   shiny::column(width = 6,
                                                                                                 shiny::uiOutput(paste0(x, "_spinner_weightDist"))
                                                                                   ),
                                                                                   shiny::column(width = 6,
                                                                                                 shiny::uiOutput(paste0(x, "_spinner_nodeDegreeDist"))
                                                                                   )
                                                                               )
                                                               ),
                                                               shiny::tabPanel(title = 'Network',
                                                                               shiny::fluidRow(
                                                                                   shiny::column(width = 8,
                                                                                                 shiny::uiOutput(paste0(x, "_spinner"))
                                                                                   ),
                                                                                   shiny::column(width = 4,
                                                                                                 shiny::uiOutput(paste0(x, '_info')),
                                                                                                 shiny::tags$hr(style='border-top: 1px solid white;'),
                                                                                                 shiny::uiOutput(paste0(x, '_custom_fnet')),
                                                                                   )
                                                                               ),
                                                                               shiny::fluidRow(
                                                                                   shiny::column(width = 12, shiny::uiOutput(paste0("update_f_net_plot_", x)), align = 'right')
                                                                               )
                                                               )
                                       )
                                   )
            )
            return(tab)
        })

        thetabs$width <- 12
        thetabs$title <- actionButton(inputId = 'update_all_nets', label = 'Submit')
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
            shinydashboard::box(shiny::h3(shiny::span('Download', style = "font-weight: bold")), solidHeader = TRUE, collapsible = FALSE, width = 12,
                                checkboxGroupInput('download_fnet',
                                                   label = shiny::h4(shiny::span('Select one or more filtered networks', style = "font-weight: bold")),
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
                file_name <- glue::glue("{obj}.csv")
                readr::write_csv(as.data.frame(to_download[[obj]]), file.path(temp_directory, file_name))
            }

            zip::zip(
                zipfile = file,
                files = dir(temp_directory),
                root = temp_directory
            )
        }
    )

    ###
    output$back2P3 <- renderUI({
        if (length(shiny::reactiveValuesToList(gene_networks)) == length(omics_files())){
            actionButton('back2P3', label = 'Back', icon("paper-plane"),
                         style="color: #fff; background-color: #337ab7; border-color: #2e6da4; padding:20px")
        }else{
            return(NULL)
        }
    })

    observeEvent(input$back2P3, {
        shinydashboard::updateTabItems(session, inputId = "tabs",
                                       selected = "net_sub_1")
    })

    ###
    output$jump2P5 <- renderUI({
        if (length(shiny::reactiveValuesToList(filteredOmicsNetworks)) == length(omics_files())){
            actionButton('jump2P5', label = 'Next', icon("paper-plane"),
                         style="color: #fff; background-color: #337ab7; border-color: #2e6da4; padding:20px")
        }else{
            return(NULL)
        }
    })

    observeEvent(input$jump2P5, {
        shinydashboard::updateTabItems(session, inputId = "tabs",
                                       selected = "rwr_tab_2")
    })


    ############################################################################################
    #                             RWR: LOADING NETWORKS (OUTPUT)                          #
    ############################################################################################

    ############## outputs loading info #############

    output$omics_net_sum_info <- shiny::renderUI({
        if (length(loaded_omics_net_list()) != 0){
            n_nodes <- sapply(loaded_omics_net_list(), function(x){length(unique(c(x[[1]], x[[2]])))})
            n_edges <- sapply(loaded_omics_net_list(), nrow)
            shinydashboard::box(width = 12,
                                shiny::HTML(paste0(paste("<b>", names(loaded_omics_net_list()), "</b>"), " loaded: ",  paste("<b>", n_nodes, "</b>"), " nodes and ", paste("<b>", n_edges, "</b>"), " edges <br/>")))
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
            shinydashboard::box(width = 12,
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
                shinydashboard::box(width = 12,
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
                shinydashboard::box(width = 12,
                                    shiny::HTML(paste("<b>", p('WARNING!', style = "color:red"), "</b>")), message1, message2
                )
            }
        } else {
            shiny::HTML("<br/>")
        }
    })

    ###
    output$back2P1_from1.1 <- renderUI({
        actionButton('back2P1_from1.1', label = 'Back', icon("paper-plane"),
                     style="color: #fff; background-color: #337ab7; border-color: #2e6da4; padding:20px")
    })

    observeEvent(input$back2P1_from1.1, {
        shinydashboard::updateTabItems(session, inputId = "tabs",
                                       selected = "mat_sub_1")
    })

    ###
    output$jump2P5.1 <- renderUI({
        if (length(loaded_omics_net_list()) != 0){
            actionButton('jump2P5.1', label = 'Next', icon("paper-plane"),
                         style="color: #fff; background-color: #337ab7; border-color: #2e6da4; padding:20px")
        }else{
            return(NULL)
        }
    })

    observeEvent(input$jump2P5.1, {
        shinydashboard::updateTabItems(session, inputId = "tabs",
                                       selected = "rwr_tab_2")
    })


    ############################################################################################
    #                               RWR: MULTIPLEX NETWORK (OUTPUT)                            #
    ############################################################################################

    output$info_box_6 <- renderText({
        HTML("<br/> <span style='font-weight:normal;'>
            <p align='justify'>
            Welcome to the MoNETA Shiny app, an R package for network-based multi-omics data integration.<br/>
            In this section of the app, you can create a multiplex network,
            which is a multilayered network with as many layers as the number of networks inferred in the previous steps. </span> <br/>
            </p> <hr style='border-top: 1px solid white;'>

            <b> Weighted Multiplex </b>: <br/> <span style='font-weight:normal;'>
            <p align='justify'>
             If you choose to create a <u>weighted</u> multiplex network, the edge weights of each
             network will be retained and contribute to the node visit probability
             in the next phase of the pipeline. Subsequently,
             you can filter the multiplex network by selecting a threshold.  </span> <br/>
            </p> <hr style='border-top: 1px solid white;'>

            <b> Unweighted Multiplex </b>: <br/> <span style='font-weight:normal;'>
            <p align='justify'>
             If you decide not to weight the multiplex network, the weight column will be filled with <b>1</b>,
             and consequently, each <span style='color:red;'>intra-layer </span> edge will have an equal probability of being traversed by the random walker. </span> <br/>
            </p> <hr style='border-top: 1px solid white;'>
            ")
    })

    ############## output multiplex net info #############

    output$multi_net_sum_info <- shiny::renderUI({
        input$gen_multiplex_btn
        shiny::isolate({
            multiplex <- multiplex_network()
            if (!is.null(multiplex)){
                multiplex_type <- ifelse(input$weightMultiplex == 'YES', 'Weighted', 'Unweighted')
                n_layers <- multiplex %>%  dplyr::select(EdgeType) %>% unique(.) %>% nrow(.)
                n_nodes <- base::union(multiplex$source, multiplex$target) %>% length(.)
                n_edges <- multiplex %>% nrow(.)
                shinydashboard::box(width = 12,
                                    shiny::HTML(paste0(paste("<b>",multiplex_type, "Multiplex network:</b>"), "<br/>")),
                                    shiny::HTML(paste('&ensp;', "<b>", n_layers, "</b>", " layers, <br/>")),
                                    shiny::HTML(paste('&ensp;',"<b>", n_nodes, "</b>", " nodes <br/>")),
                                    shiny::HTML(paste('&ensp;',"<b>", n_edges, "</b>", "edges <br/>"))
                )
            } else {
                shiny::HTML("<br/>")
            }
        })
    })

    ###
    output$back2P4 <- renderUI({
        if (length(shiny::reactiveValuesToList(filteredOmicsNetworks)) == length(omics_files())){
            actionButton('back2P4', label = 'Back', icon("paper-plane"),
                         style="color: #fff; background-color: #337ab7; border-color: #2e6da4; padding:20px")
        }else{
            return(NULL)
        }
    })
    output$back2P1.1 <- renderUI({
        if (!is.null(loaded_omics_net_list())){
            actionButton('back2P1.1', label = 'Back', icon("paper-plane"),
                         style="color: #fff; background-color: #337ab7; border-color: #2e6da4; padding:20px")
        }else{
            return(NULL)
        }
    })

    observeEvent(input$back2P4, {
        shinydashboard::updateTabItems(session, inputId = "tabs",
                                       selected = "net_sub_2")
    })

    observeEvent(input$back2P1.1, {
        shinydashboard::updateTabItems(session, inputId = "tabs",
                                       selected = "rwr_tab_1")
    })

    ###
    output$jump2P6 <- renderUI({
        if (isTruthy(input$gen_multiplex_btn)){
            actionButton('jump2P6', label = 'Next', icon("paper-plane"),
                         style="color: #fff; background-color: #337ab7; border-color: #2e6da4; padding:20px")}else{
                             return(NULL)
                         }
    })

    observeEvent(input$jump2P6, {
        shinydashboard::updateTabItems(session, inputId = "tabs",
                                       selected = "rwr_tab_3")
    })

    ############################################################################################
    #                                    RWR: PARAMETERS (OUTPUT)                              #
    ############################################################################################

    output$info_box_7 <- renderText({
        HTML("<br/> <span style='font-weight:normal;'>
            <p align='justify'>
            Welcome to the MoNETA Shiny app, an R package for network-based multi-omics data integration.<br/>
            Here, you can integrate multi-omics network through the Random Walk with Restart (RWR) algorithm,
            tath simulates the traversal of an imaginary particle. <br/>
            <hr style='border-top: 1px solid white;'>

            The random walker can:<br/>
            i) walk within a layer of the multiplex through <b><span style='color:red;'>intra-layer</span></b> edges; <br/>
            ii) return back to the <span style='color:red;'>seed</span>  (<b>Restart</b> panel); <br/>
            iii) jump to another node in a different layer (<b>Transition</b> panel). <br/>
            </p> </span> <hr style='border-top: 1px solid white;'>

            <b> Restart panel</b>: <br/> <span style='font-weight:normal;'>
            <p align='justify'>
                 The parameter <b>r</b> represents the probability for the random walker
                 to return back to the seed, while the vector parameter <b>&tau;</b> measures
                 the probability of restarting in the seed of each layer. If you <u>use
                 default restarting probabilities</u>, all layers will have the same probability of being visited during the restart. <br/>
            </p> </span>
             <hr style='border-top: 1px solid white;'>

            <b>Transition panel</b>: <br/> <span style='font-weight:normal;'>
            <p align='justify'>
                The parameter <b>&delta;</b> represents the general probability for the random walker
                to jump between layers of the multiplex network. <br/>
                &ensp; Actually, it is possible to specified the jumping probability
                between all possible pairs of layers by modifing the <b>Transition Layer Matrix</b>.
                If the <u>biologically informed</u> strategy is selected, the random walker can transition between
                omics layers with probability proportional to the number of shared relationships. <br/>
                &ensp; By default, the <b><span style='color:red;'>inter-layers</span></b> edges connect nodes representing the same entity across layers, but
                if you select the <b>RWR-NF</b> option, the edge set will also include their neighborhood across layers.
                Only if RWR-NF is selected, it is possible to weight these edges.
                For further details, please refer to the following paper
                <span style= 'color: blue;'> <u>  <a href='https://bmcbioinformatics.biomedcentral.com/articles/10.1186/s12859-021-04029-3' target='_blank'> (link)</a></u></span>.<br/>
             </p> </span>
             <hr style='border-top: 1px solid white;'>
            ")
    })

    #### Taus
    output$tauBIO <- shiny::renderUI({
        if (input$tao_opt == 'Custumize restarting probabilities per layer'){
            net_list <- if (!is.null(g_net_list())) g_net_list() else l_net_list()
            numInput <- if (!is.null(g_net_list())) length(g_net_list()) else length(l_net_list())
            lapply(1:numInput, function(i) {
                shiny::numericInput(
                    inputId = paste0('omics_tau', i),
                    label = shiny::HTML("&tau;", i, "(", names(net_list)[i] ,")"),
                    min = 0,
                    max = 1,
                    value = 0,
                    step = 0.1)
            })
        }else{
            return(NULL)
        }
    })

    omics_tau_list <- shiny::reactive({
        if (input$tao_opt != 'Custumize restarting probabilities per layer'){
            return(NA)
        }else{
            numInput <- if (!is.null(g_net_list())) length(g_net_list()) else length(l_net_list())
            lapply(1:numInput, function(i) {
                input[[paste0("omics_tau", i)]]
            })
        }
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

    ###
    output$back2P5 <- renderUI({

        if (isTruthy(input$gen_multiplex_btn)){
            actionButton('back2P5', label = 'Back', icon("paper-plane"),
                         style="color: #fff; background-color: #337ab7; border-color: #2e6da4; padding:20px")
        }else{
            return(NULL)
        }
    })

    observeEvent(input$back2P5, {
        shinydashboard::updateTabItems(session, inputId = "tabs",
                                       selected = "rwr_tab_2")
    })

    ###
    output$jump2P7 <- renderUI({
        if (length(shiny::reactiveValuesToList(RWR_output)) != 0){
            actionButton('jump2P7', label = 'Next', icon("paper-plane"),
                         style="color: #fff; background-color: #337ab7; border-color: #2e6da4; padding:20px")
        }else{
            return(NULL)
        }
    })

    observeEvent(input$jump2P7, {
        shinydashboard::updateTabItems(session, inputId = "tabs",
                                       selected = "dr_tab")
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
                                shiny::numericInput(inputId = 'threads', label = 'Select the number of CORES', value = 1)
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

                output$dr_plot <- plotly::renderPlotly({
                    isolate({
                        MoNETA::plot_2D_matrix(coord = dr_mat[1:2,], nodes_anno = data.frame(id = colnames(dr_mat)),
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
                    })

                })

                return(shinydashboard::box(width = 8,
                                           shinycssloaders::withSpinner(plotly::plotlyOutput('dr_plot')),
                                           shiny::hr(),
                                           shiny::downloadButton(outputId = 'download_dr_plot', label = 'Download plot')))
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
                file_name <- if (obj == 'dr_mat') glue::glue("{input$dr_method}_{obj}.csv") else glue::glue("{obj}.csv")
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

    ###
    output$back2P6 <- renderUI({
        if (length(shiny::reactiveValuesToList(RWR_output)) != 0){
            actionButton('back2P6', label = 'Back', icon("paper-plane"),
                         style="color: #fff; background-color: #337ab7; border-color: #2e6da4; padding:20px")
        }else{
            return(NULL)
        }
    })

    observeEvent(input$back2P6, {
        shinydashboard::updateTabItems(session, inputId = "tabs",
                                       selected = "rwr_tab_3")
    })

    ###
    output$jump2P8 <- renderUI({
        if (length(shiny::reactiveValuesToList(dr_output)) != 0){
            actionButton('jump2P8', label = 'Next', icon("paper-plane"),
                         style="color: #fff; background-color: #337ab7; border-color: #2e6da4; padding:20px")
        }else{
            return(NULL)
        }
    })

    observeEvent(input$jump2P8, {
        shinydashboard::updateTabItems(session, inputId = "tabs",
                                       selected = "cl_tab")
    })



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
        }else{
            return(NULL)
        }
    })

    ############## output cluster box #############

    output$cl_opt <- shiny::renderUI({
        anno_name <- if (length(annotation()) == 0) NULL else 'annotation'
        anno1_name <-  if (length(annotation1()) == 0) NULL else 'annotation1'
        anno_list <- list(anno_name, anno1_name)

        output$shape_opt <- shiny::renderUI({
            if ( ( !is.null(anno()) |  !is.null(anno1()) ) & input$cluster_method %in% c('kmeans', 'dbscan', 'hclust') ){
                anno_table <- if ( !is.null(anno()) ) anno() else anno1()
                output <-  tagList()
                output[[1]] <-
                    shiny::radioButtons(inputId = 'shape_bool', label = 'Do you want to use particular shape for points?', c('YES', 'NO'), selected = 'NO')
                output[[2]] <-
                    shiny::conditionalPanel(
                        condition = 'input.shape_bool == "YES"',
                        shiny::selectInput(inputId = 'select_shape', label = 'Select a column to shape data points',  choices = colnames(anno_table))
                    )
                output
            }else if (input$cluster_method %in% c('kmeans', 'dbscan', 'hclust') & (  is.null(anno()) & is.null(anno1()) ) ){
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
                    output$cl_plot <- plotly::renderPlotly({
                        isolate({
                            MoNETA::plot_2D_matrix(coord = clusters$dr_mat[1:2,], nodes_anno = clu_anno,
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
                        })
                    })
                } else{
                    output$cl_plot <- plotly::renderPlotly({
                        isolate({
                            MoNETA::plot_2D_matrix(coord = clusters$dr_mat[1:2,], nodes_anno = clu_anno,
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
                        })
                    })
                }
                return(shinydashboard::box(width = 8,
                                           shinycssloaders::withSpinner(plotly::plotlyOutput('cl_plot')),
                                           shiny::hr(),
                                           shiny::downloadButton(outputId = 'download_cl_plot', label = 'Download plot')))
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

                        shinyalert::shinyalert(
                            title = "Next step",
                            text = paste0("i) Look at the <b> dendrogram </b> that was just displayed; <br/>
                                          ii) <b> move the slider </b> to select a threshold <b> 'h' </b> to cut the dendrogram; <br/>
                                          iii) click the <b> 'Run' </b> button again."
                            ),
                            closeOnEsc = TRUE,
                            closeOnClickOutside = TRUE,
                            html = TRUE,
                            type = "info",
                            showConfirmButton = TRUE,
                            confirmButtonText = "OK",
                            confirmButtonCol = "#004192",
                            showCancelButton = FALSE,
                            imageUrl = "",
                            animation = TRUE
                        )

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
            shinydashboard::box(width = 8,
                                shinycssloaders::withSpinner(shiny::plotOutput(outputId = "plotHclust")))
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

    ###
    output$back2P7 <- renderUI({
        if (length(shiny::reactiveValuesToList(dr_output)) != 0){
            actionButton('back2P7', label = 'Back', icon("paper-plane"),
                         style="color: #fff; background-color: #337ab7; border-color: #2e6da4; padding:20px")
        }else{
            return(NULL)
        }
    })

    observeEvent(input$back2P7, {
        shinydashboard::updateTabItems(session, inputId = "tabs",
                                       selected = "dr_tab")
    })

}

