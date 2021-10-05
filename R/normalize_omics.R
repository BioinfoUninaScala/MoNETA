#' normalize_omics takes in input a Matrix that has to be centered and scaled along the columns
#'
#' @param matrix A numeric Matrix representing a particular omics, rows are genes and columns are samples
#' @return An object of the same type as mtx. The matrix mtx is normalize along the columns
#' @export

normalize_omics <- function(matrix) {

    scale(matrix)

}
