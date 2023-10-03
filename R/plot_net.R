#' Plot omics network with annotations
#'
#' @importFrom igraph V
#' @importFrom network as.color
#' @import visNetwork
#' @import tidyverse
#' @import RColorBrewer
#' @import dplyr
#' @param edgeList a data frame containing the edge list of a network in the first two columns containing node ids.
#' @param nodes_anno annotation data frame with nodes on rows and annotations on columns
#' @param id_name name of the column in \emph{nodes_anno} containing unique sample ids
#' @param id_anno_color name of the column in \emph{nodes_anno} whose values are mapped to colors to observations. Notice that colors can be recycled it eh number of distinct values exceed available shapes
#' @param id_anno_shape name of the column in \emph{nodes_anno} whose values are mapped to shapes to observations. Notice that shapes can be recycled it eh number of distinct values exceed available shapes
#' @param title title of plot
#' @param html logical, if TRUE the plot is produced in html format
#' @param wo_legend logical, if TRUE returns the plot without legend
#' @param interactive logical, if TRUE returns an interactive plot
#' @export


plot_net <- function (edgeList, nodes_anno, id_name, id_anno_color = NA, id_anno_shape = NA, title = "", html = FALSE, wo_legend = TRUE, interactive = TRUE)
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

    if (!interactive) {
        p = p %>%
            visInteraction(dragNodes = FALSE, dragView = FALSE, zoomView = FALSE)
    }

    if(html) {
        p = p %>% visNetwork::visSave(file = paste0(title, ".html"))
    }

    p
}
