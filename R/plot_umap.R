#' Plot umap
#'
#' @importFrom igraph V
#' @import visNetwork
#' @import tidyverse
#' @import RColorBrewer
#' @import FNN
#' @import ggrepel
#' @import fpc
#' @import ggfortify
#' @import survival
#' @import survminer
#' @importFrom plotly orca
#' @import dplyr
#' @import tibble
#' @importFrom ggplot2 aes ggplot geom_point theme_bw labs theme
#' @importFrom umap umap.knn
#' @param matrix A Matrix with samples on columns that has to be plotted
#' @param nodes_anno Annotation DataFrame with all the information for each sample
#' @param id_name String for identification of sample in nodes_anno
#' @param id_anno_color String for the column necessary for distinguish cases in the nodes_anno, it will be used for giving a specific color for each case
#' @param id_anno_shape String for the column necessary for distinguish cases in the nodes_anno, it will be used for giving a specific shape for each case, it must be associated to a discrete variable
#' @param interactive A boolean flag, if TRUE returns an interactive plot
#' @param wo_legend A boolean flag, if TRUE returns the plot without legend
#' @param title Title of the plot
#' @param umap_coord It represents previously computed umap coordinates, if not present, it will be computed
#' @param shapes It is possible to give the shapes list for plotting, id_anno_shape is not NA
#' @return A plot
#' @export

plot_umap <- function(matrix, nodes_anno, id_name, id_anno_color = NA, id_anno_shape = NA,
                      interactive = TRUE, wo_legend = FALSE, title = "",
                      umap_coord = NULL, shapes = NULL) {

    matrix <- t(matrix)

    nodes_anno <- nodes_anno %>% dplyr::filter(nodes_anno[[id_name]] %in%
                                            base::rownames(matrix))



    if (is.na(id_anno_color) & is.na(id_anno_shape)) {
        nodes_anno <- nodes_anno %>% dplyr::select(dplyr::all_of(id_name))
        colnames(nodes_anno) <- "id"
    } else if (!is.na(id_anno_color) & !is.na(id_anno_shape)) {
        nodes_anno <- nodes_anno %>% dplyr::select(dplyr::all_of( c(id_name, id_anno_color, id_anno_shape) ))
        colnames(nodes_anno) <- c("id", "group1", "group2")
    } else if (!is.na(id_anno_color)) {
        nodes_anno <- nodes_anno %>% dplyr::select(dplyr::all_of( c(id_name, id_anno_color) ))
        colnames(nodes_anno) <- c("id", "group1")
    } else if (!is.na(id_anno_shape)) {
        nodes_anno <- nodes_anno %>% dplyr::select(dplyr::all_of( c(id_name, id_anno_shape) ))
        colnames(nodes_anno) <- c("id", "group2")
    }


    if(is.null(umap_coord)) {
        umap_coord <- umap::umap(d = matrix)
    } else {
        umap_coord <- t(umap_coord)
    }


    gPlot_umap_data <- umap_coord$layout
    colnames(gPlot_umap_data) = c("V1", "V2")
    gPlot_umap_data <- gPlot_umap_data %>% tibble::as_tibble(rownames = "id", .name_repair = 'unique')
    gPlot_umap_data = gPlot_umap_data %>% dplyr::left_join(nodes_anno, by = "id")

    if (!is.na(id_anno_shape)) {
        if (is.null(shapes)) {
            vals <- plotly::schema(F)$traces$scatter$attributes$marker$symbol$values
            vals <- grep("[a-z]", vals, value = T)

            split_pos <- grep("-open-dot", vals)
            opendot = vals[split_pos]
            vals = vals[-split_pos]

            split_pos <- grep("-dot", vals)
            dot = vals[split_pos]
            vals = vals[-split_pos]

            split_pos <- grep("-open", vals)
            open = vals[split_pos]
            vals = vals[-split_pos]

            px = which(vals == "x")
            psq = which(vals == "square")

            vals[px] = "square"
            vals[psq] = "x"

            vals = c(vals, opendot, dot, open)

            vals = vals[!(stringr::str_ends(vals, "down") | stringr::str_ends(vals, "up") |
                              stringr::str_ends(vals, "left") | stringr::str_ends(vals, "right") |
                              stringr::str_ends(vals, "open") | stringr::str_ends(vals, "up") |
                              stringr::str_ends(vals, "ne") | stringr::str_ends(vals, "se") |
                              stringr::str_ends(vals, "sw") | stringr::str_ends(vals, "nw") )]
        } else {
            vals = shapes
        }


    }




    if (is.na(id_anno_color) & is.na(id_anno_shape)) {
        if (interactive) {
            gPlot_umap <- gPlot_umap_data %>% plotly::plot_ly(x = ~V1, y = ~V2, type = "scatter",
                                                              mode = "markers",
                                                              marker = list(size = 5), text = ~id,
                                                              colors = "Set2") %>%
                plotly::layout(xaxis = list(zeroline = F), yaxis = list(zeroline = F),
                               title = title)
        } else {
            gPlot_umap <- ggplot(gPlot_umap_data, aes(x = V1, y = V2)) +
                geom_point(color = "blue")  +
                theme_bw() +
                labs(title = title)
        }
    } else if (!is.na(id_anno_color) & !is.na(id_anno_shape)) {
        if (interactive) {
            gPlot_umap <- gPlot_umap_data %>%
                plotly::plot_ly(x = ~V1, y = ~V2, type = "scatter",
                                color = ~group1, mode = "markers", symbol = ~group2,
                                symbols = vals,
                                marker = list(size = 5), text = ~id,
                                name = paste0(gPlot_umap_data$group1, "\n", gPlot_umap_data$group2),
                                colors = "Set2") %>%
                plotly::layout(xaxis = list(zeroline = F), yaxis = list(zeroline = F),
                               title = title)
        } else {
            gPlot_umap <- ggplot(gPlot_umap_data, aes(x = V1, y = V2)) +
                geom_point(aes(color = group1, shape = group2))  +
                theme_bw() +
                labs(title = title)
        }
    } else if (!is.na(id_anno_color)) {
        if (interactive) {
            gPlot_umap <- gPlot_umap_data %>% plotly::plot_ly(x = ~V1, y = ~V2, type = "scatter",
                                                              color = ~group1, mode = "markers",
                                                              marker = list(size = 5), text = ~id,
                                                              colors = "Set2") %>%
                plotly::layout(xaxis = list(zeroline = F), yaxis = list(zeroline = F),
                               title = title)
        } else {
            gPlot_umap <- ggplot(gPlot_umap_data, aes(x = V1, y = V2)) +
                geom_point(aes(color = group1))  +
                theme_bw()  +
                labs(title = title)
        }
    } else if (!is.na(id_anno_shape)) {
        if (interactive) {
            gPlot_umap <- gPlot_umap_data %>%
                plotly::plot_ly(x = ~V1, y = ~V2, type = "scatter",
                                mode = "markers", symbol = ~group2,
                                symbols = vals,
                                marker = list(size = 5), text = ~id,
                                name = paste0(gPlot_umap_data$group2),
                                colors = "Set2") %>%
                plotly::layout(xaxis = list(zeroline = F), yaxis = list(zeroline = F),
                               title = title)
        } else {
            gPlot_umap <- ggplot(gPlot_umap_data, aes(x = V1, y = V2)) +
                geom_point(color = "blue", aes(shape = group2))  +
                theme_bw() +
                labs(title = title)
        }
    }


    if (!interactive & wo_legend) {
        gPlot_umap <- gPlot_umap + theme(legend.position = "none")
    }

    gPlot_umap


}
