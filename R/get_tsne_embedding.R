#' Compute an embedding of the multi-omics matrix using t-SNE
#'
#' @importFrom Rtsne Rtsne
#' @param matrix A squared numeric matrix with values in the range 0-1, with samples on columns
#' @param embedding_size Size of the output embedding
#' @param perplexity Perplexity parameter
#' @param max_iter Number of iterations
#' @param num_threads Number of threads for parallelization. It has to be positive integer or 0, if all available cores must be used. If it is equal to 1, no parallelization is not performed
#' @return An embedding of the input Matrix: on the columns there are the samples, the number of rows are specified by embedding_size
#' @export

get_tsne_embedding <- function(matrix, embedding_size, perplexity, max_iter, num_threads) {

    embed <- t(Rtsne::Rtsne(t(abs(1-matrix)), dims = embedding_size,
                     perplexity=perplexity, verbose=FALSE,
                     max_iter = max_iter, is_distance = T,
                     num_threads = num_threads)$Y)
    colnames(embed) = colnames(matrix)
    embed

}
