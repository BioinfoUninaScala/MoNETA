#' Compute an embedding of the multi-omics matrix using uwot umap
#'
#' @importFrom umap umap
#' @param matrix A squared numeric matrix with values in the range 0-1, with samples on columns
#' @param embedding_size Size of the output embedding
#' @param n_neighbors The size of local neighborhood. Larger values result in more global views of the manifold, while smaller values result in more local data being preserved. In general values should be in the range 2 to 100
#' @param n_threads Number of thread, NULL if single core
#' @param n_sgd_threads Number of thread for stochastic gradient descent, if > 1 it will not be reproducible anyway
#' @param grain_size Minimum amount of work to do on each thread
#' @return An embedding of the input Matrix: on the columns there are the samples, the number of rows are specified by embedding_size
#' @export

get_parallel_umap_embedding <- function(matrix, embedding_size, n_neighbors = 15, n_threads = NULL, n_sgd_threads = 0, grain_size = 1) {
    embed = uwot::umap(t(matrix), n_components = embedding_size, n_neighbors = n_neighbors,
                       n_threads = n_threads, n_sgd_threads = n_sgd_threads, grain_size = grain_size)
    embed = t(embed)
    colnames(embed) = colnames(matrix)
    embed
}
