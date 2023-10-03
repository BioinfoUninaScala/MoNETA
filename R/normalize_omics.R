#' Center and scale an omic matrix along its columns
#'
#' @param matrix a numeric matrix representing a particular omics, features are on rows and samples are on columns
#' @return a matrix center and scaled along its columns
#' @export

normalize_omics <- function(matrix) {

    scale(matrix)

}
