#' Create jump matrix
#'
#' @importFrom dplyr inner_join
#' @param multiplex a multiplex
#' @return jump_mat a jump matrix
#' @export

create_jump_matrix <- function(multiplex) {

    layers <- unique(multiplex$EdgeType)
    jump_mat <- matrix(0, nrow = length(layers), ncol = length(layers), dimnames = list(layers,layers))

    comb <- utils::combn(layers,2)
    for (i in 1:ncol(comb))
        jump_mat[ comb[1,i],comb[2,i] ] <- nrow(dplyr::inner_join(net_list[[ comb[1,i] ]], net_list[[ comb[2,i] ]], by = c("source","dest"))) /
        (nrow(net_list[[ comb[1,i] ]]) + nrow(net_list[[ comb[2,i] ]]))

    jump_mat[lower.tri(jump_mat)] <- t(jump_mat)[lower.tri(jump_mat)]

    # base:: o Matrix::
    jump_mat <- jump_mat/base::rowSums(jump_mat)
    jump_mat
}
