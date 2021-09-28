#' Plot network
#'
#' @import igraph
#' @import visNetwork
#' @import tidyverse
#' @import network
#' @import RColorBrewer
#' @import dplyr
#' @param edgeList DataFrame representing a graph
#' @param nodes_anno Annotation DataFrame with all the information for each sample
#' @param name Title of plot to be displayed
#' @param html Boolean, if you want to save the plot in a html file
#' @export


plot_net <- function (edgeList, nodes_anno, name, html = FALSE)
{
    base::colnames(edgeList) <- c("from", "to")
    graph <- igraph::graph_from_data_frame(edgeList, directed = FALSE)
    cluster <- igraph::cluster_louvain(graph)
    cluster_df <- igraph::membership(cluster) %>% tibble::enframe()
    base::colnames(cluster_df) <- c("id", "group")
    nodes_anno <- nodes_anno %>% dplyr::filter(nodes_anno[[1]] %in%
                                            base::union(edgeList[[1]], edgeList[[2]]))
    nodes <- data.frame(id = nodes_anno[[1]], anno = nodes_anno[[2]])
    nodes <- dplyr::left_join(nodes, cluster_df)
    nodes$color <- RColorBrewer::brewer.pal(n = length(unique(network::as.color(nodes$anno))),
                              name = "Dark2")[network::as.color(nodes$anno)]
    if (!(html)) {
        visNetwork::visNetwork(nodes, edgeList, main = name) %>% visNetwork::visIgraphLayout() %>%
            visNetwork::visOptions(highlightNearest = list(enabled = TRUE,
                                               degree = 3), selectedBy = "anno")
    }
    else {
        visNetwork::visNetwork(nodes, edgeList, main = name) %>% visNetwork::visLegend(main = name) %>%
            visNetwork::visIgraphLayout() %>% visNetwork::visOptions(highlightNearest = list(enabled = TRUE,
                                                                     degree = 3), selectedBy = "anno") %>% visNetwork::visSave(file = paste0(name,
                                                                                                                                 ".html"))
    }
}
