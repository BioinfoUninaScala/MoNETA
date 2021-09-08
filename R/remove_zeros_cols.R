#' Create jump matrix
#'
#' @param mtx a matrix
#' @return the same matrix as input without the zeros columns
#' @export


remove_zeros_cols <- function(mtx) {
    mtx[, colSums(mtx) > 0]
}
