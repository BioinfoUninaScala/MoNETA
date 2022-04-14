#' Plot parallel umap db
#'
#' @importFrom uwot umap
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
#' @import ggplot2
#' @import plotly
#' @import umap
#' @import dplyr
#' @param matrix A DataFrame with samples on columns that has to be plotted
#' @param k Reachability distance, dbscan
#' @param interactive A boolean flag, it TRUE returns an interactive plot
#' @param n_neighbors The size of local neighborhood. Larger values result in more global views of the manifold, while smaller values result in more local data being preserved. In general values should be in the range 2 to 100
#' @param n_threads Number of thread, NULL if single core
#' @param n_sgd_threads Number of thread for stochastic gradient descent, if > 1 it will not be reproducible anyway
#' @param grain_size Minimum amount of work to do on each thread
#' @export



plot_parallel_umap_db <- function(matrix, k=1, interactive = TRUE,
                                  n_neighbors = 15, n_threads = NULL, n_sgd_threads = 0,
                                  grain_size = 1){
    matrix = as.data.frame(matrix)
    cols = colnames(matrix)
    matrix = t(matrix)
    umap_coord <- uwot::umap(matrix, n_components = 2, n_neighbors = n_neighbors,
                             n_threads = n_threads, n_sgd_threads = n_sgd_threads,
                             grain_size = grain_size)
    rownames(umap_coord) = cols
    gPlot_umap_data <- umap_coord %>% tibble::as_tibble(rownames = "id")

    ds.norm = fpc::dbscan(umap_coord, k)
    gPlot_umap_data$clust = factor(ds.norm$cluster)
    hc.norm.cent = gPlot_umap_data %>% dplyr::group_by(clust) %>% dplyr::select(V1,
                                                                                V2, clust) %>% dplyr::summarize_all(mean)
    if (interactive) {
        gPlot_umap_data_key <- plotly::highlight_key(gPlot_umap_data, ~id)
    } else {
        gPlot_umap_data_key <- gPlot_umap_data
    }

    gPlot_umap_cl <- (ggplot(gPlot_umap_data_key, aes(x = V1,
                                                      y = V2, colour = clust)) + geom_point(alpha = 0.3) +
                          theme_bw() + ggrepel::geom_label_repel(aes(label = clust), data = hc.norm.cent) +
                          ggtitle(paste("DBSCAN epslion:", k)))

    gPlot_umap_cl

}
