#' UMAP embedding of matrix
#'
#' @importFrom umap umap.knn umap.defaults
#' @param matrix a squared n n x n numeric matrix with values in the range 0-1, with observations on columns
#' @param embedding_size size of the output embedding
#' @return  a \emph{embedding_size} x n matrix having observations on columns
#' @export

get_umap_embedding <- function(matrix, embedding_size) {
    config = umap.defaults
    config$n_components = embedding_size

    embed = umap::umap(t(matrix), config = config)
    embed = t(embed$layout)
    embed
}
