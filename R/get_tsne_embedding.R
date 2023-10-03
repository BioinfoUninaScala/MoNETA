#' TSNE embedding of matrix
#'
#' @importFrom Rtsne Rtsne
#' @param matrix a squared n n x n numeric matrix with values in the range 0-1, with observations on columns
#' @param embedding_size size of the output embedding
#' @param perplexity perplexity parameter
#' @param max_iter number of iterations
#' @param num_threads number of threads for parallelization. It has to be positive integer or 0, if all available cores must be used. If it is equal to 1, no parallelization is not performed
#' @return a \emph{embedding_size} x n matrix having observations on columns
#' @export

get_tsne_embedding <- function(matrix, embedding_size, perplexity, max_iter, num_threads) {

    embed <- t(Rtsne::Rtsne(t(abs(1-matrix)), dims = embedding_size,
                     perplexity=perplexity, verbose=FALSE,
                     max_iter = max_iter, is_distance = T,
                     num_threads = num_threads)$Y)
    colnames(embed) = colnames(matrix)
    embed

}
