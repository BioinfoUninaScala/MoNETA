#' PCA embedding of matrix
#'
#' @param matrix a squared n n x n numeric matrix with values in the range 0-1, with observations on columns
#' @param embedding_size size of the output embedding
#' @return a \emph{embedding_size} x n matrix having observations on columns
#' @export

get_pca_embedding <- function(matrix, embedding_size) {
    t(prcomp(t(matrix))$x)[1:embedding_size,]
}
