#' Plot umap db
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
#' @param matrix A Matrix with samples on rows that has to be plotted
#' @param type Type of matrix, used as input for umap, see umap documentation
#' @param k Reachability distance, dbscan
#' @export



plot_umap_db <- function(matrix, type="data", k=1){

    matrix = t(matrix)
    umap_coord <- umap::umap(d = matrix, input = type)
    gPlot_umap_data <- umap_coord$layout %>% tibble::as_tibble(rownames = "id")

    ds.norm = fpc::dbscan(umap_coord$layout, k)
    gPlot_umap_data$clust = factor(ds.norm$cluster)
    hc.norm.cent = gPlot_umap_data %>% dplyr::group_by(clust) %>% dplyr::select(V1,
                                                                  V2) %>% dplyr::summarize_all(mean)
    gPlot_umap_data_key <- ifelse(flag, plotly::highlight_key(gPlot_umap_data, ~id), gPlot_umap_data)

    gPlot_umap_cl <- (ggplot(gPlot_umap_data_key, aes(x = V1,
                                                      y = V2, colour = clust)) + geom_point(alpha = 0.3) +
                          theme_bw() + geom_label_repel(aes(label = clust), data = hc.norm.cent) +
                          ggtitle(paste(id, "DBSCAN epslion:", k)))

    gPlot_umap_cl

}
