#' Prune low weight edges from network
#'
#' @param networka data frame containing the edge list of a network with two columns containing nodes ids and one column named weight containing edge weights
#' @param weight_remv a threshold value on the edge weight. Edges with weight less then \emph{weight_remv} will be removed from the network.
#' @return data frame contianing pruned network as an edge list
#' @export

prune_multiplex_network <- function(network, weight_remv = NULL) {

    if (is.null(weight_remv)) {
        weight_remv <- max(network[, "weight"])
    }

    network[network[, "weight"] < weight_remv,]

}
