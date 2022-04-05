#' Create jump matrix containing the probabilities to switch from one omics to another
#'
#' @importFrom dplyr inner_join
#' @param omics_list A list of DataFrames containings omics data, each is a DataFrame with three columns: source, destination, weight
#' @return A Matrix containing the probabilities to switch from one omics to another
#' @export

create_jump_matrix <- function(omics_list) {

    # @param multiplex A multiplex, an object that comprises different omics. It is a DataFrame containing these information: EdgeType, source, target, weight
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



