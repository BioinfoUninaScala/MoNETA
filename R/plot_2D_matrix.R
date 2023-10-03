#' Plot bimedinational reduction of a matrix
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
#' @importFrom ggplot2 aes scale_color_brewer ggplot geom_point theme_bw labs theme
#' @importFrom stringr str_ends
#' @param coord a data frame with x and y coordinates on rows and observations on columns
#' @param nodes_anno annotation data frame with observations on rows and annotations on columns
#' @param id_name name of the column in \emph{nodes_anno} containing unique sample ids. There has to be correspondence on \emph{coord} columns names
#' @param id_anno_color  name of the column in \emph{nodes_anno} whose values are mapped to colors to observations. Notice that colors can be recycled it eh number of distinct values exceed available shapes
#' @param id_anno_shape name of the column in \emph{nodes_anno} whose values are mapped to shapes to observations. Notice that shapes can be recycled it eh number of distinct values exceed available shapes
#' @param interactive logical, if TRUE returns a plotly interactive plot, otherwise a standard ggplot is returned
#' @param wo_legend logical, if TRUE returns the plot without legend
#' @param title title of the plot
#' @param shapes  shapes list to be used for plotting, id_anno_shape is not NA
#' @return a ggplot or a plotly object
#' @export

plot_2D_matrix <- function(coord, nodes_anno, id_name, id_anno_color = NA, id_anno_shape = NA,
                               interactive = TRUE, wo_legend = FALSE, title = "", shapes = NULL) {



    coord <- t(coord)
    colnames(coord) <- c("V1","V2")

    nodes_anno <- nodes_anno %>% dplyr::filter(nodes_anno[[id_name]] %in%
                                                   base::rownames(coord))

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




    gPlot_data <- coord %>% tibble::as_tibble(rownames = "id", .name_repair = 'unique') %>%
        dplyr::left_join(nodes_anno, by = "id")

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


    gPlot_data[is.na(gPlot_data)] = "NA"

    if (is.na(id_anno_color) & is.na(id_anno_shape)) {
        if (interactive) {
            gPlot <- gPlot_data %>% plotly::plot_ly(x = ~V1, y = ~V2, type = "scatter",
                                                              mode = "markers",
                                                              marker = list(size = 5), text = ~id,
                                                              colors = "Set2") %>%
                plotly::layout(xaxis = list(zeroline = F), yaxis = list(zeroline = F),
                               title = title)
        } else {
            gPlot <- ggplot(gPlot_data, aes(x = V1, y = V2)) +
                geom_point(color = "blue")  +
                theme_bw() +
                labs(title = title)
        }
    } else if (!is.na(id_anno_color) & !is.na(id_anno_shape)) {
        if (interactive) {
            gPlot <- gPlot_data %>%
                plotly::plot_ly(x = ~V1, y = ~V2, type = "scatter",
                                color = ~group1, mode = "markers", symbol = ~group2,
                                symbols = vals,
                                marker = list(size = 5), text = ~id,
                                name = paste0(gPlot_data$group1, "\n", gPlot_data$group2),
                                colors = "Set2") %>%
                plotly::layout(xaxis = list(zeroline = F), yaxis = list(zeroline = F),
                               title = title)
        } else {
            gPlot <- ggplot(gPlot_data, aes(x = V1, y = V2)) +
                geom_point(aes(color = group1, shape = group2))  +
                theme_bw() +
                labs(title = title, color = id_anno_color, shape = id_anno_shape) +
                scale_color_brewer(palette="Set2")
        }
    } else if (!is.na(id_anno_color)) {
        if (interactive) {
            gPlot <- gPlot_data %>% plotly::plot_ly(x = ~V1, y = ~V2, type = "scatter",
                                                              color = ~group1, mode = "markers",
                                                              marker = list(size = 5), text = ~id,
                                                              colors = "Set2") %>%
                plotly::layout(xaxis = list(zeroline = F), yaxis = list(zeroline = F),
                               title = title)
        } else {
            gPlot <- ggplot(gPlot_data, aes(x = V1, y = V2)) +
                geom_point(aes(color = group1))  +
                theme_bw()  +
                labs(title = title, color = id_anno_color) +
                scale_color_brewer(palette="Set2")
        }
    } else if (!is.na(id_anno_shape)) {
        if (interactive) {
            gPlot <- gPlot_data %>%
                plotly::plot_ly(x = ~V1, y = ~V2, type = "scatter",
                                mode = "markers", symbol = ~group2,
                                symbols = vals,
                                marker = list(size = 5), text = ~id,
                                name = paste0(gPlot_data$group2),
                                colors = "Set2") %>%
                plotly::layout(xaxis = list(zeroline = F), yaxis = list(zeroline = F),
                               title = title)
        } else {
            gPlot <- ggplot(gPlot_data, aes(x = V1, y = V2)) +
                geom_point(color = "blue", aes(shape = group2))  +
                theme_bw() +
                labs(title = title, shape = id_anno_shape)
        }
    }


    if (!interactive & wo_legend) {
        gPlot <- gPlot + ggplot2::theme(legend.position = "none")
    }

    gPlot

}
