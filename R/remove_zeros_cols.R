#' Remove zero columns from a matrix
#'
#' @param matrix a numeric matrix
#' @return a matrix of the same type of \emph{matrix} with columns containing all zero removed
#' @export


remove_zeros_cols <- function(matrix) {
    matrix[, colSums(matrix) > 0]
}
