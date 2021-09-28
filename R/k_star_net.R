#' Apply K Start algorithm in order to find knn nearest neighbors for each node, using Vantage-point tree
#'
#' @importFrom purrr array_branch map_dbl pmap_df map_chr
#' @import dplyr
#' @import tibble
#' @import magrittr
#' @importFrom stats median setNames
#' @importFrom doMC registerDoMC
#' @importFrom BiocParallel MulticoreParam
#' @importFrom BiocNeighbors buildIndex queryKNN VptreeParam KmknnParam
#' @importFrom foreach registerDoSEQ %dopar%
#' @param mat A numeric Matrix representing a particular omics, rows are genes and columns are samples
#' @param distFun A string that represents a distance function, two possible choices: "Manhattan" and "Euclidean"
#' @param sparsity A positive real
#' @param knn An integer, it is the number of neighbors to be considered
#' @param cores Number of threads for Parallelization. It has to be positive integer. If it is equal to 1, no parallelization is not performed
#' @return A Tibble with three columns : source node, destination node, weight of connection
#' @export

k_star_net <- function(mat, distFun = "Euclidean", sparsity = 1, knn = 25, cores = 1) {

    index <- BiocNeighbors::buildIndex(t(mat), BNPARAM = BiocNeighbors::VptreeParam(distance=distFun))

    if (cores > 1) {
        doMC::registerDoMC(cores)
        knns <- BiocNeighbors::queryKNN(BNINDEX = index, query = t(mat), k = knn + 1,
                                        BPPARAM = BiocParallel::MulticoreParam(cores))
    } else {
        foreach::registerDoSEQ()
        knns <- BiocNeighbors::queryKNN(BNINDEX = index, query = t(mat), k = knn + 1)
    }

    knn_elems <- foreach::foreach(i=1:dim(mat)[2], .export = c("search_k_star_nn")) %dopar%
        search_k_star_nn(colnames(mat)[knns$index[i,]], knns$distance[i,2:(knn+1)], sparsity = sparsity)

    if (cores > 1) {
        foreach::registerDoSEQ()
    }

    dplyr::bind_rows(knn_elems)

}


search_k_star_nn <- function(id, delta, sparsity = 1){

    beta <- delta * sparsity
    Sum_beta <- 0
    Sum_beta_square <- 0
    lambda <-  beta[1] + 1
    k <- 0

    while (lambda > beta[k + 1] && k < (length(beta))){
        k = k + 1;
        Sum_beta <- Sum_beta + beta[k]
        Sum_beta_square <- Sum_beta_square + beta[k] ^ 2

        lambda <- 1/k * ( Sum_beta + sqrt(k + Sum_beta ^ 2 - k * Sum_beta_square) )
    }

    edges <- tibble::tibble(source = id[1], dest = id[2:k], weight = delta[2:k])

    edges
}
