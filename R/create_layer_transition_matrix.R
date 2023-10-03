#' Creates matrix containing transition probabilities between omics layers of multiples network.
#'
#' @importFrom dplyr inner_join
#' @param omics_list
#' a named list of data frames containing an omics multiplex network.
#' Each data frame in the list is required to have at least two columns containing nodes ids respectively named source and destination
#' and one column named weight containing edge weights.
#' @return A squared n x n matrix (where n is the number of omics layers) containing  transition probabilities between omics layers.
#' @export

create_layer_transition_matrix <- function(omics_list) {

    # @param multiplex A multiplex, an object that comprises different omics.

    #layers <- unique(multiplex$EdgeType)
    layers <- names(omics_list)
    jump_mat <- matrix(0, nrow = length(layers), ncol = length(layers), dimnames = list(layers,layers))

    comb <- utils::combn(layers,2)
    for (i in 1:ncol(comb))
        jump_mat[ comb[1,i],comb[2,i] ] <- max(nrow(dplyr::inner_join(omics_list[[ comb[1,i] ]], omics_list[[ comb[2,i] ]], by = c("source","dest"))), 1) /
        (nrow(omics_list[[ comb[1,i] ]]) + nrow(omics_list[[ comb[2,i] ]]))

    jump_mat[lower.tri(jump_mat)] <- t(jump_mat)[lower.tri(jump_mat)]

    jump_mat <- jump_mat/base::rowSums(jump_mat)
    jump_mat
}



