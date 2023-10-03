#' Build omic similarity network
#' Publication: https://papers.nips.cc/paper/2016/file/2c6ae45a3e88aee548c0714fad7f8269-Paper.pdf
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
#' @param matrix a numeric matrix representing a particular omics, features are on rows and samples are on columns
#' @param distFun a distance function chosen among: "Euclidean", "Manhattan", "Cosine"
#' @param sparsity a positive real value guiding the choice of the number of neighbors for each node based on the distribution of distances of its neighbors
#' @param knn an integer, it is the maximum number of neighbors to be considered as candidate neighbors for each node
#' @param k_star if TRUE  the number of neighbors id dynamically determined for each node using kstar algorithm ( https://papers.nips.cc/paper/2016/file/2c6ae45a3e88aee548c0714fad7f8269-Paper.pdf),
#' if FALSE k closest neighbors are chosen for each node
#' @param cores number of threads for parallelization. It has to be positive integer. If it is equal to 1, no parallelization is not performed
#' @param MAX_ASSOC number of maximum incoming edges that a node can have
#' @return a tibble containing the edge list of output network with the following columns: source, dest, weight.
#' Source and dest represent the edges and contain node ids, weight contains edge weights.
#' @export

k_star_net <- function(matrix, distFun = "Euclidean", sparsity = 1, knn = 25, k_star = TRUE, cores = 1, MAX_ASSOC = Inf) {

    index <- BiocNeighbors::buildIndex(t(matrix), BNPARAM = BiocNeighbors::VptreeParam(distance=distFun))


    if (cores > 1) {
        doMC::registerDoMC(cores)
        knns <- BiocNeighbors::queryKNN(BNINDEX = index, query = t(matrix), k = knn + 1,
                                        BPPARAM = BiocParallel::MulticoreParam(cores))
    } else {
        foreach::registerDoSEQ()
        knns <- BiocNeighbors::queryKNN(BNINDEX = index, query = t(matrix), k = knn + 1)
    }

    knn_elems <- foreach::foreach(i=1:dim(knns$index)[1], .export = c("search_k_star_nn")) %dopar% {
        ids = index@NAMES[knns$index[i,-which(knns$index[i,] == i)]]
        id = index@NAMES[knns$index[i,which(knns$index[i,] == i)]]
        delta = knns$distance[i,-which(knns$index[i,] == i)]

        if (k_star)
            search_k_star_nn(id = id, ids = ids, delta = delta, sparsity = sparsity)
        else
            tibble::tibble(source = id, dest = ids, weight = delta)

    }


    if (cores > 1) {
        foreach::registerDoSEQ()
    }


    knn_elems = dplyr::bind_rows(knn_elems)

    if (is.infinite(MAX_ASSOC)) {
        MAX_ASSOC = nrow(knn_elems)
    }

    knn_elems %>% group_by(dest) %>%
        arrange(weight, .by_group = TRUE) %>%
        dplyr::slice(1:MAX_ASSOC)

}


search_k_star_nn <- function(id, ids, delta, sparsity = 1){

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

    edges <- tibble::tibble(source = id, dest = ids[1:k], weight = delta[1:k])

    edges
}





########################################## OLD CODE






build_Vptreefrom_mat <- function (mat, distFun){
    build_Vptree(points = array_branch(mat,2),distFun = distFun)
}


build_Vptree <- function (points, distFun){

    if (length(points) == 0)
        stop('Points can not be empty.')

    if (length(points) == 1 )
        return (new_node(points[[1]], names(points)[1]))

    # Vantage point is point furthest from parent vp.
    vp_i <-  new_node(points[[1]], names(points)[1])
    points <-  points[-1]

    # Choose division boundary at median of distances.
    distances <-  points %>% map_dbl(distFun, vp_i$vp)
    median_d <- median(distances)
    vp_i$median_d <- median_d
    left_points <- list()
    right_points <- list()

    for (i in 1:length(points)){
        if (distances[i] >= median_d){
            vp_i$right_min = min(distances[i], vp_i$right_min)
            if (distances[i] > vp_i$right_max){
                vp_i$right_max = distances[i]
                right_points <- c(points[i],right_points) # put furthest first
            }else{
                right_points <- c(right_points,points[i])
            }
        }else{
            vp_i$left_min = min(distances[i], vp_i$left_min)
            if (distances[i] > vp_i$left_max){
                vp_i$left_max = distances[i]
                left_points <- c(points[i],left_points) # put furthest first
            }else{
                left_points <- c(left_points,points[i])
            }
        }
    }

    if (length(left_points) > 0)
        vp_i$left = build_Vptree(points=left_points, distFun=distFun)

    if (length(right_points) > 0)
        vp_i$right = build_Vptree(points=right_points, distFun=distFun)

    vp_i
}



new_node <- function(vp, vp_id){

    self <- list()
    self$vp <- vp
    self$left = NULL
    self$right = NULL
    self$left_min = Inf
    self$left_max = 0
    self$right_min = Inf
    self$right_max = 0
    self$curr_d = 0
    self$id <- vp_id

    self
}

get_neigh <- function(id, q, vp_t, k, distFun, sparsity = 1){
    # print(id)
    # print(head(q))
    top_vp <- get_n_nearest_neighbors(q = q, tree = vp_t, k = k + 1 , distFun = distFun)
    if (id %in% names(top_vp))
        top_vp <- top_vp[!names(top_vp) == id]

    search_k_star_nn_old(id = id, delta = top_vp, sparsity = sparsity)
}


get_n_nearest_neighbors <- function(q, tree, k, distFun){


    tau = Inf # furthest disteance collected when having at least k neighbours
    nodes_to_visit = list(tree)

    # fixed size array for nearest neightbors
    # sorted from closest to farthest neighbor
    neighbors=list(nodes=c(), keys=c())


    while (length(nodes_to_visit) > 0){
        node = nodes_to_visit[[1]]
        nodes_to_visit <- nodes_to_visit[-1]

        if(is.null(node))
            next

        d = distFun(q, node$vp)
        #print(paste(d,tau))
        if (length(neighbors) < k || d < tau){
            # store node.vp as a neighbor if it's closer than any other of the best k points
            # seen so far
            neighbors <- append_neigh(neighbors, list(node=node,key= d), k)

            # start shrinking tau when at least k neighbors have been found
            if (d < tau && length(neighbors) == k){
                farthest_nearest_neighbor = neighbors$nodes[[length(neighbors$nodes)]]
                tau = distFun(q, farthest_nearest_neighbor$vp)
            }
        }


        if(is_leaf(node))
            next
        # check for intersection between q-tau and vp-mu regions
        # and see which branches we absolutely must search

        if (d < node$median_d){
            if (d < node$median_d + tau)
                nodes_to_visit <- append(nodes_to_visit, list(node$left))
            if (d >= node$median_d - tau)
                nodes_to_visit <- append(nodes_to_visit, list(node$right))
        } else {
            if (d >= node$median_d - tau)
                nodes_to_visit <- append(nodes_to_visit, list(node$right))
            if (d < node$median_d + tau)
                nodes_to_visit <- append(nodes_to_visit, list(node$left))
        }
    }

    neighbors <- setNames(neighbors$keys, neighbors$nodes %>% map_chr(~.$id))
    neighbors
}


search_k_star_nn_old <- function(id , delta, sparsity = 1){

    #min-max normalise distances and multiply by 100
    #delta <- (delta-min(delta))/(max(delta)-min(delta))

    beta <- delta * sparsity
    Sum_beta <- 0
    Sum_beta_square <- 0
    lambda <-  beta[1] + 1
    k <- 0
    while (lambda > beta[k+1] && k < (length(beta)-1)){
        k = k+1;
        Sum_beta <- Sum_beta + beta[k]
        Sum_beta_square <- Sum_beta_square + beta[k]^2

        lambda <- 1/k * ( Sum_beta + sqrt(k + Sum_beta^2 - k*Sum_beta_square) )

    }


    edges <- tibble(source = id, dest = names(delta)[1:k], weight = delta[1:k])

    edges
}

append_neigh <- function(neighbors, new_n, max_k){

    neighbors$keys <- c(neighbors$keys, new_n$key)
    neighbors$nodes <- c(neighbors$nodes, list(new_n$node))

    order <- order(neighbors$keys)
    neighbors$keys <-  neighbors$keys[order]
    neighbors$nodes <-  neighbors$nodes[order]

    if (length(neighbors$keys) > max_k){
        neighbors$keys <- neighbors$keys[-length(neighbors$key)]
        neighbors$nodes <- neighbors$nodes[-length(neighbors$nodes)]
        #print(neighbors$keys)
    }

    neighbors
}

is_leaf <- function(vpt){
    is.null(vpt$left) && is.null(vpt$right)
}


