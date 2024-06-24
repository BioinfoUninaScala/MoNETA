#' UMAP embedding of matrix using uwot umap
#'
#' @importFrom uwot tumap
#' @param matrix a squared n n x n numeric matrix with values in the range 0-1, with observations on columns
#' @param embedding_size size of the output embedding
#' @param n_neighbors size of local neighborhood. Larger values result in more global views of the manifold, while smaller values result in more local data being preserved. In general values should be in the range 2 to 100
#' @param n_threads number of threads, NULL if single core
#' @param n_sgd_threads number of thread for stochastic gradient descent, if > 1 it will not be reproducible
#' @param grain_size minimum amount of work to do on each thread
#' @param metric distance metric to use with umap
#' @return a \emph{embedding_size} x n matrix having observations on columns
#' @export

get_parallel_umap_embedding <- function(matrix, embedding_size, n_neighbors = 15, n_threads = NULL, n_sgd_threads = 0, grain_size = 1, metric = 'cosine') {
    embed = uwot::umap(t(matrix), n_components = embedding_size, n_neighbors = n_neighbors,
                       n_threads = n_threads, n_sgd_threads = n_sgd_threads, grain_size = grain_size,
                       target_metric = metric)
    embed = t(embed)
    colnames(embed) = colnames(matrix)
    embed
}
