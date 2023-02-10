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
#' @param title Title of plot to be displayed
#' @param html Boolean, if you want to save the plot in a html file
#' @export


plot_net <- function (edgeList, nodes_anno, title, html = FALSE)
{
    base::colnames(edgeList)[1:2] <- c("from", "to")
    graph <- igraph::graph_from_data_frame(edgeList, directed = FALSE)
    nodes_anno <- nodes_anno %>% dplyr::filter(nodes_anno[[1]] %in%
                                            base::union(edgeList[[1]], edgeList[[2]]))
    nodes <- data.frame(id = nodes_anno[[1]], anno = nodes_anno[[2]])
    colors_disp = colorRampPalette(RColorBrewer::brewer.pal(n = 8,
                                           name = "Dark2"))(length(unique(network::as.color(nodes$anno))))
    nodes$color <- colors_disp[network::as.color(nodes$anno) %% length(colors_disp) +1]
    if (!(html)) {
        visNetwork::visNetwork(nodes, edgeList, main = title) %>% visNetwork::visIgraphLayout() %>%
            visNetwork::visOptions(highlightNearest = list(enabled = TRUE,
                                               degree = 3), selectedBy = "anno")
    } else {
        visNetwork::visNetwork(nodes, edgeList, main = title) %>% visNetwork::visLegend(main = title) %>%
            visNetwork::visIgraphLayout() %>% visNetwork::visOptions(highlightNearest = list(enabled = TRUE,
                                                                     degree = 3), selectedBy = "anno") %>% visNetwork::visSave(file = paste0(title,
                                                                                                                                 ".html"))
    }
}
