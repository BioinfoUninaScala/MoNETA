#' This function clean a Matrix removing the columns with all zeros
#'
#' @param matrix A Matrix
#' @return A Matrix of the same type of mtx without columns with all zeros
#' @export


remove_zeros_cols <- function(matrix) {
    matrix[, colSums(matrix) > 0]
}
