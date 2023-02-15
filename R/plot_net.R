#' Plot network
#'
#' @importFrom igraph V
#' @importFrom network as.color
#' @import visNetwork
#' @import tidyverse
#' @import RColorBrewer
#' @import dplyr
#' @param edgeList DataFrame representing a graph
#' @param nodes_anno Annotation DataFrame with all the information for each sample
#' @param id_name String for identification of sample in nodes_anno
#' @param id_anno_color String for the column necessary for distinguish cases in the nodes_anno, it will be used for giving a specific color for each case
#' @param id_anno_shape String for the column necessary for distinguish cases in the nodes_anno, it will be used for giving a specific shape for each case, it must be associated to a discrete variable
#' @param title Title of plot to be displayed
#' @param html Boolean, if you want to save the plot in a html file
#' @param wo_legend A boolean flag, if TRUE returns the plot without legend
#' @export


plot_net <- function (edgeList, nodes_anno, id_name, id_anno_color = NA, id_anno_shape = NA, title = "", html = FALSE, wo_legend = TRUE)
{
    base::colnames(edgeList)[1:2] <- c("from", "to")
    graph <- igraph::graph_from_data_frame(edgeList, directed = FALSE)

    nodes_anno <- nodes_anno %>% dplyr::filter(nodes_anno[[id_name]] %in%
                                            base::union(edgeList[[1]], edgeList[[2]]))

    nodes <- data.frame(id = nodes_anno[[id_name]])

    if (!is.na(id_anno_color)) {
        nodes$id_color = nodes_anno[[id_anno_color]]
        colors_disp = colorRampPalette(RColorBrewer::brewer.pal(n = 8,
                                                                name = "Dark2"))(length(unique(network::as.color(nodes$id_color))))
        nodes$color <- colors_disp[network::as.color(nodes$id_color) %% length(colors_disp) +1]

        lnodes <- unique(nodes[, c("id_color", "color"), drop = FALSE])
        lnodes$label = lnodes$id_color
    }
    if (!is.na(id_anno_shape)) {
        nodes$id_shape = nodes_anno[[id_anno_shape]]
        shapes = c("square", "triangle", "box", "circle", "dot", "star",
                  "ellipse", "database", "text", "diamond")

        unique_gshape_len = length(unique(nodes$id_shape))

        shapes = rep(shapes, ceiling(unique_gshape_len / length(shapes)))

        nodes$shape = shapes[match(nodes$id_shape, unique(nodes$id_shape))]

    }







    p = visNetwork::visNetwork(nodes, edgeList, main = title) %>%
        visNetwork::visIgraphLayout() %>%
        visNetwork::visNodes(labelHighlightBold = TRUE, font = list(size=20))


    if(!is.na(id_anno_color)) {
        p = p %>% visNetwork::visOptions(
                                        highlightNearest = list(enabled = TRUE),
                                         selectedBy = "id_color")
    }
    if(!wo_legend) {
        p = p %>% visNetwork::visLegend(useGroups = FALSE, addNodes = lnodes)
    }

    if(html) {
        p = p %>% visNetwork::visSave(file = paste0(title, ".html"))
    }

    p
}
