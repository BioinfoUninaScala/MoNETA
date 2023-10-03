#' Multi-omics matrix embedding
#'
#' Adapted from:
#' Research group:  Léo Pio-Lopez, Alberto Valdeolivas, Laurent Tichit, Élisabeth Remy, Anaïs Baudot
#' Publication:     https://arxiv.org/abs/2008.10085
#'
#' @import foreach
#' @import doParallel
#' @import bigmemory
#' @import bigstatsr
#' @import Matrix
#' @import wordspace
#' @importFrom doMC registerDoMC
#' @importFrom foreach registerDoSEQ %dopar%
#' @param matrix a n x n squared numeric similarity matrix with values in the range 0-1
#' @param embedding_size size of the output embedding
#' @param num_steps number of total epoches
#' @param cores number of threads for parallelization. It has to be positive integer. If it is equal to 1, no parallelization is performed
#' @return a \emph{embedding_size} x n numerical matrix
#' @export
#'
#'
get_embedding <- function(matrix, embedding_size, num_steps = 10 ^ 7, cores = 20) {

    cols <- base::colnames(matrix)
    base::colnames(matrix) <- 1:ncol(matrix)#as.character(1:ncol(RWR_mat_plot))
    base::rownames(matrix) <- 1:nrow(matrix)#as.character(1:nrow(RWR_mat_plot))
    embedding <- get_embed(matrix, embedding_size, num_steps, cores)
    embedding <- t(embedding)
    base::colnames(embedding) <- cols
    embedding
}


get_embed <- function(matrix, embedding_size, num_steps = 10 ^ 7, cores) {
    CLOSEST_NODES = 20
    NUM_SAMPLED = 3
    LEARNING_RATE = 0.01
    NB_CHUNK = 1
    CHUNK_SIZE = 10

    r_DistancematrixPPI <- matrix #read_csv("data/mat_1.csv")

    lst <- preprocess(r_DistancematrixPPI, CLOSEST_NODES)

    nodes <- lst$nodes
    neighborhood <- lst$neighborhood
    list_neighbours <- lst$list_neighbours
    reverse_data_DistancematrixPPI <- lst$reverse_data_DistancematrixPPI

    embeddings <- matrix(stats::rnorm(length(nodes) * embedding_size, mean = 0., sd = 1),
                         nrow = length(nodes), ncol = embedding_size)

    emb <- train(neighborhood, nodes, list_neighbours, num_steps, NUM_SAMPLED, LEARNING_RATE,
                 CLOSEST_NODES, CHUNK_SIZE, NB_CHUNK, embeddings, reverse_data_DistancematrixPPI,
                 cores)

    emb
}

preprocess <- function(r_DistancematrixPPI, CLOSEST_NODES) {
    rawdata_DistancematrixPPI <- t(as.matrix(r_DistancematrixPPI))
    node_size <- nrow(rawdata_DistancematrixPPI)


    #neighborhood : a list representing the size of neighborhoob of each node which are significantly similar
    neighborhood <- apply(rawdata_DistancematrixPPI, 1, function(x) sum(x > 1/node_size))


    #change the diagonal as the lowest positive score of each node
    rawdata_DistancematrixPPI <- as.matrix(r_DistancematrixPPI)
    diag(rawdata_DistancematrixPPI) <- apply(rawdata_DistancematrixPPI, 1, FUN = function(x) min(x[x > 0]))

    rawdata_DistancematrixPPI <- t(rawdata_DistancematrixPPI)

    # ?normalize.rows
    rawdata_DistancematrixPPI <- wordspace::normalize.rows(rawdata_DistancematrixPPI, method = "manhattan")


    nodes <- base::colnames(r_DistancematrixPPI)

    #list_neighbours : list of neighbours sorted from the most similar ones to each node
    #reverse_data : list of maximum similarity (length CLOSEST_NODES) scores for each node
    sort_genes <- apply(rawdata_DistancematrixPPI, 1, order, decreasing = TRUE)
    base::colnames(sort_genes) <-  NULL
    sort_genes <- apply(sort_genes, 1, function(x) nodes[x])
    list_neighbours <- sort_genes#[,1:CLOSEST_NODES]


    sort_values <- apply(rawdata_DistancematrixPPI, 1, sort, decreasing = TRUE)
    base::colnames(sort_values) <-  NULL
    reverse_data_DistancematrixPPI <- sort_values[, 1:CLOSEST_NODES]

    return(list("neighborhood" = neighborhood, "nodes" = nodes,
                "list_neighbours" = list_neighbours, "reverse_data_DistancematrixPPI" = reverse_data_DistancematrixPPI))
}

node_positive_weighted <- function(u, list_neighbours, CLOSEST_NODES, reverse_data_DistancematrixPPI) {
    #if the similarities of a node are all zeros, it takes one randomly
    #else it takes the one random element from the most similar with a probability higher for the most similar ones
    if (sum(reverse_data_DistancematrixPPI[u,]) == 0){
        as.integer(list_neighbours[u, sample(1:CLOSEST_NODES, 1)])
    } else {
        probas = reverse_data_DistancematrixPPI[u, 1:CLOSEST_NODES]
        as.integer(sample(list_neighbours[u, 1:CLOSEST_NODES], 1, prob = probas))
    }
}

sigmoid <- function(x) {
    1 / (1 + exp(-x))
}

update <- function(W_u, W_v, D, learning_rate, bias) {
    sim <- sigmoid(W_u %*% W_v - bias)
    gradient <- as.vector((D - sim) * learning_rate)
    W_u <- W_u + gradient *  W_v
    W_v <- W_v + gradient *  W_u
    list("W_u" = W_u, "W_v" = W_v)
}


train <- function(neighborhood, nodes, list_neighbours, NUM_STEPS, NUM_SAMPLED, LEARNING_RATE,
                  CLOSEST_NODES, CHUNK_SIZE, NB_CHUNK, embd, reverse_data_DistancematrixPPI,
                  cores) {

    embeddings <- bigstatsr::as_FBM(embd)

    if (cores > 1) {
        #cl <- parallel::makeCluster(cores)
        #doParallel::registerDoParallel(cl)
        doMC::registerDoMC(cores)
    } else {
        foreach::registerDoSEQ()
    }


    nb_nodes <- length(nodes)
    nce_bias <- log(nb_nodes)
    nce_bias_neg <- log(nb_nodes / NUM_SAMPLED)

    #lock <- tempfile()


    foreach (steps = 1:cores, .export = c("update", "sigmoid", "node_positive_weighted")) %dopar% {

        for (k in 1:floor(NUM_STEPS/cores)) {
            # select randomly CHUNK_SIZE nodes
            nodes_opt <- sample(1:nb_nodes, CHUNK_SIZE, replace=FALSE)
            #if (k %% floor(10**6/CHUNK_SIZE) == 0) {
            if(k %% 1000 == 0){ print(paste("Intermediate step: ", steps, "   Inner step: ", k)) }

            for (i in 1:CHUNK_SIZE) {
                u <- nodes_opt[i]
                v <- node_positive_weighted(u, list_neighbours, CLOSEST_NODES, reverse_data_DistancematrixPPI)
                tmp_lst <- update(embeddings[u,], embeddings[v,], 1, LEARNING_RATE, nce_bias)
                #locked <- flock::lock(lock)
                embeddings[u,] <- tmp_lst$W_u
                embeddings[v,] <- tmp_lst$W_v
                #flock::unlock(locked)

                for (j in 1:NUM_SAMPLED) {
                    v_neg <- as.integer(list_neighbours[u, sample(CLOSEST_NODES+1:nb_nodes-CLOSEST_NODES,
                                                                  1, replace = FALSE)])
                    tmp_lst <- update(embeddings[u,], embeddings[v_neg,], 0, LEARNING_RATE, nce_bias_neg)
                    #locked <- flock::lock(lock)
                    embeddings[u,] <- tmp_lst$W_u
                    embeddings[v_neg,] <- tmp_lst$W_v
                    #flock::unlock(locked)
                }
            }
        }
        #print(steps * k)
    }

    if (cores > 1) {
        #parallel::stopCluster(cl)
        foreach::registerDoSEQ()
    }

    embeddings[]


}
