#' Plot umap
#'
#' @import igraph
#' @import visNetwork
#' @import tidyverse
#' @import network
#' @import RColorBrewer
#' @import FNN
#' @import ggrepel
#' @import fpc
#' @import ggfortify
#' @import survival
#' @import survminer
#' @import plotly
#' @import umap
#' @import dplyr
#' @import tibble
#' @param matrix A Matrix with samples on rows that has to be plotted
#' @param nodes_anno Annotation DataFrame with all the information for each sample
#' @param id_name String for identification of sample in nodes_anno
#' @param id_anno_color String for the column necessary for distinguish cases in the nodes_anno, it will be used for giving a specific color for each case
#' @param id_anno_shape String for the column necessary for distinguish cases in the nodes_anno, it will be used for giving a specific shape for each case
#' @param title Title of the plot
#' @return A plot
#' @export

plot_umap <- function(matrix, nodes_anno, id_name, id_anno_color = NA, id_anno_shape = NA, title = "") {

    matrix <- t(matrix)
    nodes_anno <- nodes_anno
    nodes_anno <- nodes_anno %>% dplyr::filter(nodes_anno[[id_name]] %in%
                                            base::rownames(matrix))
    #base::colnames(nodes_anno)[1:3] <- c("id", "group", "group2")
    if (is.na(id_anno_color) & is.na(id_anno_shape)) {
        nodes_anno <- nodes_anno %>% select(id_name)
        colnames(nodes_anno) <- "id"
        umap_coord <- umap::umap(d = matrix)
        gPlot_umap_data <- umap_coord$layout %>% tibble::as_tibble(rownames = "id") %>%
            dplyr::left_join(nodes_anno, by = "id")
        gPlot_umap <- gPlot_umap_data %>% plotly::plot_ly(x = ~V1, y = ~V2, type = "scatter",
                                                          mode = "markers",
                                                          marker = list(size = 5), text = ~id) %>%
            plotly::layout(xaxis = list(zeroline = F), yaxis = list(zeroline = F),
                           title = title)

    } else if (!is.na(id_anno_color) & !is.na(id_anno_shape)) {
        nodes_anno <- nodes_anno %>% select(id_name, id_anno_color, id_anno_shape)
        colnames(nodes_anno) <- c("id", "group1", "group2")
        umap_coord <- umap::umap(d = matrix)
        gPlot_umap_data <- umap_coord$layout %>% tibble::as_tibble(rownames = "id") %>%
            dplyr::left_join(nodes_anno, by = "id")
        gPlot_umap <- gPlot_umap_data %>% plotly::plot_ly(x = ~V1, y = ~V2, type = "scatter",
                                                          color = ~group1, mode = "markers", symbol = ~group2,
                                                          marker = list(size = 5), text = ~id) %>%
            plotly::layout(xaxis = list(zeroline = F), yaxis = list(zeroline = F),
                           title = title)

    } else if (!is.na(id_anno_color)) {
        nodes_anno <- nodes_anno %>% select(id_name, id_anno_color)
        colnames(nodes_anno) <- c("id", "group1")
        umap_coord <- umap::umap(d = matrix)
        gPlot_umap_data <- umap_coord$layout %>% tibble::as_tibble(rownames = "id") %>%
            dplyr::left_join(nodes_anno, by = "id")
        gPlot_umap <- gPlot_umap_data %>% plotly::plot_ly(x = ~V1, y = ~V2, type = "scatter",
                                                          color = ~group1, mode = "markers",
                                                          marker = list(size = 5), text = ~id) %>%
            plotly::layout(xaxis = list(zeroline = F), yaxis = list(zeroline = F),
                           title = title)

    } else if (!is.na(id_anno_shape)) {
        nodes_anno <- nodes_anno %>% select(id_name, id_anno_shape)
        colnames(nodes_anno) <- c("id", "group2")
        umap_coord <- umap::umap(d = matrix)
        gPlot_umap_data <- umap_coord$layout %>% tibble::as_tibble(rownames = "id") %>%
            dplyr::left_join(nodes_anno, by = "id")
        gPlot_umap <- gPlot_umap_data %>% plotly::plot_ly(x = ~V1, y = ~V2, type = "scatter",
                                                          symbol = ~group2, mode = "markers",
                                                          marker = list(size = 5), text = ~id) %>%
            plotly::layout(xaxis = list(zeroline = F), yaxis = list(zeroline = F),
                           title = title)

    }

    gPlot_umap



}
