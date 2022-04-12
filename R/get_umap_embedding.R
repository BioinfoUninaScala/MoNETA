#' Compute an embedding of the multi-omics matrix using umap
#'
#' @importFrom umap umap
#' @param matrix A squared numeric matrix with values in the range 0-1, with samples on columns
#' @param embedding_size Size of the output embedding
#' @return An embedding of the input Matrix: on the columns there are the samples, the number of rows are specified by embedding_size
#' @export

get_umap_embedding <- function(matrix, embedding_size) {
    config = umap.defaults
    config$n_components = embedding_size

    embed = umap::umap(t(matrix), config = config)
    embed = t(embed$layout)
    embed
}
