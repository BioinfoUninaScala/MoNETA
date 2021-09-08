#' Create jump matrix
#'
#' @import doFuture
#' @import igraph
#' @import mclust
#' @import Matrix
#' @import kernlab
#' @import R.matlab
#' @import bc3net
#' @import optparse
#' @import parallel
#' @import tidyverse
#' @import doParallel
#' @param network a network
#' @param tau tau
#' @param restart restart
#' @param delta delta
#' @param cond_jump cond_jump
#' @param cores number of threads launched
#' @return RWRM_similarity
#' @export

#source("Functions_RWRM.R")

gen_sim_mat_M <- function(network, tau = NA, restart = 0.7, delta = 0.5, cond_jump = NULL, cores = 1){

    InputNetwork <- network

    if (ncol(InputNetwork) != 4){
        stop("The input network should has 4 columns. See help for format details",
             call. = FALSE)
    } else {
        base::colnames(InputNetwork) <- c("EdgeType", "source", "target", "weight")
        Number_Layers <- length(unique(InputNetwork$EdgeType))
        print(paste0("Input multiples network with ", Number_Layers, " Layers"))
        ## We check that the weights are numeric.
        if (!is.numeric(InputNetwork$weight)){
            stop("The weights in the input network should be numeric")
        }
    }

    if (restart > 1 || restart < 0){
        stop("Restart parameter should range between 0 and 1")
    }

    if (length(tau) == 1 && is.na(tau)){
        tau <- rep(1,Number_Layers)/Number_Layers
    }
    print(tau)

    # if (is.null(out)){
    #   stop("You need to specify a name to be used to generate the output files.")
    # }

    # MachineCores <- detectCores()
    #
    # if (cores < 1 || cores > MachineCores){
    #   stop("The number of cores should be between 1 and the total number of cores of
    #        your machine")
    # }


    #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### ####
    #### Data Transformation and associated calculations to apply RWR
    #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### ####

    ## We transform the input multiplex format to a L-length list of igraphs objects.
    ## We also scale the weigths of every layer between 1 and the minimun
    ## divided by the maximun weight.

    LayersNames <- unique(InputNetwork$EdgeType)
    Layers <- lapply(LayersNames, function(x) {
        Original_weight <- InputNetwork$weight[InputNetwork$EdgeType == x]
        if (all(!is.na(Original_weight)) && any(Original_weight > 0) && !length(unique(Original_weight)) == 1){
            b <- 1
            a <- min(Original_weight) / max(Original_weight)
            range01 <- (b - a) * (Original_weight - min(Original_weight)) / (max(Original_weight) - min(Original_weight)) + a
            InputNetwork$weight[InputNetwork$EdgeType == x] <- range01
        }else {
            InputNetwork$weight[InputNetwork$EdgeType == x] <- 1
        }

        igraph::simplify(igraph::graph_from_data_frame(InputNetwork[InputNetwork$EdgeType == x, 2:4],
                                               directed = FALSE), edge.attr.comb=mean)
    })
    names(Layers) <- LayersNames


    ## In case the network is monoplex, we have to be aware of isolated nodes.
    if (Number_Layers == 1){
        message("Dealing with isolated nodes...")
        IsolatedVertex <- igraph::V(Layers[[1]])$name[which(igraph::degree(Layers[[1]]) == 0)]
        mynetPrev <- Layers[[1]]
        Layers[[1]] <- igraph::delete.vertices(igraph::simplify(Layers[[1]]),
                                               igraph::degree(Layers[[1]]) == 0)
    }

    ## We prepate the data to compute RWR (Compute the adjacency matrix of the
    ## multiplex network and its normalalisation)

    MultiplexObject <- create.multiplex(Layers)
    AdjMatrix <- compute.adjacency.matrix2(MultiplexObject, delta, cond_jump)
    AdjMatrixNorm <- normalize.multiplex.adjacency(AdjMatrix)
    Allnodes <- MultiplexObject$Pool_of_Nodes
    numberNodes <- length(Allnodes)

    #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### ####
    #### Apply RWR on every node of the multiplex network (Parallelise version)
    #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### ####

    ## Apply Random Walk with Restart for each node of the network.
    ## (Parallelise version)

    ###### Version using the package doParallel

    message("Computing RWR for every network node...")

    cl <- parallel::makeForkCluster(mc <- getOption("cl.cores", cores))
    doParallel::registerDoParallel(cl)

    Results <- foreach::foreach(i = 1:length(Allnodes), .packages = c("Matrix") ) %dopar% {
        Random.Walk.Restart.Multiplex.default(AdjMatrixNorm, MultiplexObject,
                                              Allnodes[i], r = restart, DispResults = "Alphabetic", MeanType = "Sum", tau = tau)
    }

    stopCluster(cl)


    #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### ####
    #### Generation of the output Similarity Matrix
    #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### ####

    ## Generation of the RWR similarity matrix
    message("Building Similarity Matrix...")
    RWRM_similarity <- matrix(data = 0, nrow = numberNodes, ncol = numberNodes)

    for (j in seq(numberNodes)) {
        RWRM_similarity[,j] <- Results[[j]]$RWRM_Results$Score
    }


    base::rownames(RWRM_similarity) <- Allnodes
    base::colnames(RWRM_similarity) <- Allnodes

    ## Here comes the trick. If we had isolated nodes, which have been removed
    ## from the analyses because the RWR cannot deal with them, we add them into
    ## the sim matrix with 1 on the node, 0 elsewhere in the matrix. This only
    ## applies to Monoplex networks where the particle can get trapped.



    if (Number_Layers == 1){
        numberIso <- length(IsolatedVertex)
        if (numberIso > 0){
            print("Including isolated nodes in the Similarity Matrix...")
            TotalSizeMatrix <-  numberNodes + numberIso
            Index_IsoNodes <- IsolatedVertex

            NewMatrix <- matrix(data = 0, nrow = TotalSizeMatrix, ncol = TotalSizeMatrix)
            base::rownames(NewMatrix) <- c(base::rownames(RWRM_similarity), Index_IsoNodes)
            base::colnames(NewMatrix) <- c(base::colnames(RWRM_similarity), Index_IsoNodes)

            NewMatrix[1:numberNodes, 1:numberNodes] <- RWRM_similarity
            NewMatrix[(numberNodes + 1):TotalSizeMatrix, (numberNodes + 1):TotalSizeMatrix] <- diag(1, numberIso, numberIso)

            RWRM_similarity <- NewMatrix
        }
        ## Check how do we want to order
        # NewMatrix <- NewMatrix[order(base::rownames(NewMatrix)), order(base::colnames(NewMatrix))]

        # NewMatrix <- NewMatrix[order(as.integer(base::rownames(NewMatrix))), order(as.integer(base::colnames(NewMatrix)))]
    }

    ## check ! Be carefull with this order
    if (is.numeric(base::rownames(RWRM_similarity))){
        RWRM_similarity <- RWRM_similarity[order(as.numeric(base::rownames(RWRM_similarity))),
                                           order(as.numeric(base::colnames(RWRM_similarity)))]
    } else {
        RWRM_similarity <- RWRM_similarity[order(base::rownames(RWRM_similarity)),
                                           order(base::colnames(RWRM_similarity))]
    }

    # as.matrix(RWRM_similarity)
    # # colSums(RWRM_similarity)
    # Similarity_outputfile <- paste0(out,".rds")
    # message(paste0("Saving file ", Similarity_outputfile))
    # saveRDS(RWRM_similarity,file = Similarity_outputfile)
    RWRM_similarity
}






#ipak <- function(pkg){
#    new.pkg <- pkg[!(pkg %in% installed.packages()[, "Package"])]
#    if (length(new.pkg))
#        install.packages(new.pkg, dependencies = TRUE)
#    sapply(pkg, require, character.only = TRUE, quietly = TRUE)
#}



isMultiplex <- function(x)
{
    methods::is(x, "Multiplex")
}

isMultiplexHet <- function(x)
{
    methods::is(x, "MultiplexHet")
}





compute.adjacency.matrix2 <- function(x, delta = 0.5, cond_jump = NULL)
{
    if (!isMultiplex(x) & !isMultiplexHet(x)) {
        stop("Not a Multiplex or Multiplex Heterogeneous object")
    }
    if (delta > 1 || delta < 0) {
        stop("Delta should be between 0 and 1")
    }


    N <- x$Number_of_Nodes_Multiplex
    L <- x$Number_of_Layers

    ###############################################################################
    if (!is.null(cond_jump) && (any(base::colnames(cond_jump) != names(x)[seq(L)]) || any(base::rownames(cond_jump) != names(x)[seq(L)]))) {
        stop("cond_jump elements order differs from x")
    }
    ###############################################################################
    ## We impose delta=0 in the monoplex case.
    if (L  ==1){
        delta = 0
    }

    Layers_Names <- names(x)[seq(L)]

    ###############################################################################
    if (is.null(cond_jump))
        cond_jump <- matrix(1 / (L - 1), nrow = L, ncol = L, dimnames = list(Layers_Names, Layers_Names))
    ###############################################################################

    diag(cond_jump) <- 0
    ## IDEM_MATRIX.
    Idem_Matrix <- Matrix::Diagonal(N, x = 1)

    counter <- 0
    Layers_List <- lapply(x[Layers_Names], function(x){

        counter <<- counter + 1;
        if (igraph::is_weighted(x)) {
            Adjacency_Layer <- igraph::as_adjacency_matrix(x, sparse = TRUE,
                                                    attr = "weight")
        } else {
            Adjacency_Layer <- igraph::as_adjacency_matrix(x, sparse = TRUE)
        }

        if (is.numeric(base::rownames(Adjacency_Layer))){
            Adjacency_Layer <- Adjacency_Layer[order(as.numeric(base::rownames(Adjacency_Layer))),
                                               order(as.numeric(base::colnames(Adjacency_Layer)))]
        } else {
            Adjacency_Layer <- Adjacency_Layer[order(base::rownames(Adjacency_Layer)),
                                               order(base::colnames(Adjacency_Layer))]
        }

        base::colnames(Adjacency_Layer) <-
            paste0(base::colnames(Adjacency_Layer), "_", counter)
        base::rownames(Adjacency_Layer) <-
            paste0(base::rownames(Adjacency_Layer), "_", counter)
        Adjacency_Layer
    })

    MyColNames <- unlist(lapply(Layers_List, function(x) unlist(base::colnames(x))))
    MyRowNames <- unlist(lapply(Layers_List, function(x) unlist(base::rownames(x))))
    names(MyColNames) <- c()
    names(MyRowNames) <- c()
    SupraAdjacencyMatrix <- (1 - delta) * (Matrix::bdiag(unlist(Layers_List)))
    base::colnames(SupraAdjacencyMatrix) <- MyColNames
    base::rownames(SupraAdjacencyMatrix) <- MyRowNames

    #offdiag <- (delta/(L-1))*Idem_Matrix

    i <- seq_len(L)
    Position_ini_row <- 1 + (i - 1) * N
    Position_end_row <- N + (i - 1) * N
    j <- seq_len(L)
    Position_ini_col <- 1 + (j - 1) * N
    Position_end_col <- N + (j - 1) * N

    for (i in seq_len(L)){
        for (j in seq_len(L)){
            if (j != i){
                SupraAdjacencyMatrix[(Position_ini_row[i]:Position_end_row[i]),
                                     (Position_ini_col[j]:Position_end_col[j])] <- delta * cond_jump[i,j] * Idem_Matrix
            }
        }
    }

    SupraAdjacencyMatrix <- methods::as(SupraAdjacencyMatrix, "dgCMatrix")
    return(SupraAdjacencyMatrix)
}



###??????? USEMETHOD

## Create a multiplex object
create.multiplex <- function(...){
    UseMethod("create.multiplex")
}







normalize.multiplex.adjacency <- function(x)
{
    if (!is(x,"dgCMatrix")){
        stop("Not a dgCMatrix object of Matrix package")
    }

    Adj_Matrix_Norm <- t(t(x) / (Matrix::colSums(x, na.rm = FALSE, dims = 1,
                                               sparseResult = FALSE)))

    return(Adj_Matrix_Norm)
}





Random.Walk.Restart.Multiplex.default <-
    function(x, MultiplexObject, Seeds, r = 0.7, tau, MeanType = "Geometric",
             DispResults = "TopScores", ...){

        ### We control the different values.
        if (!is(x, "dgCMatrix")){
            stop("Not a dgCMatrix object of Matrix package")
        }

        if (!isMultiplex(MultiplexObject)) {
            stop("Not a Multiplex object")
        }

        L <- MultiplexObject$Number_of_Layers
        N <- MultiplexObject$Number_of_Nodes

        Seeds <- as.character(Seeds)
        if (length(Seeds) < 1 | length(Seeds) >= N){
            stop("The length of the vector containing the seed nodes is not
           correct")
        } else {
            if (!all(Seeds %in% MultiplexObject$Pool_of_Nodes)){
                stop("Some of the seeds are not nodes of the network")
            }
        }

        if (r >= 1 || r <= 0) {
            stop("Restart partameter should be between 0 and 1")
        }

        #base o rlang
        if(base::missing(tau)){
            tau <- rep(1, L) / L
        } else {
            tau <- as.numeric(tau)
            if (sum(tau) != 1) {
                stop("The sum of the components of tau divided by the number of
             layers should be 1")
            }
        }

        if(!(MeanType %in% c("Geometric","Arithmetic","Sum"))){
            stop("The type mean should be Geometric, Arithmetic or Sum")
        }

        if(!(DispResults %in% c("TopScores","Alphabetic"))){
            stop("The way to display RWRM results should be TopScores or
           Alphabetic")
        }

        ## We define the threshold and the number maximum of iterations for
        ## the random walker.
        Threeshold <- 1e-10
        NetworkSize <- ncol(x)

        ## We initialize the variables to control the flux in the RW algo.
        residue <- 1
        iter <- 1

        ## We compute the scores for the different seeds.
        Seeds_Score <- get.seed.scoresMultiplex(Seeds, L, tau)

        ## We define the prox_vector(The vector we will move after the first RWR
        ## iteration. We start from The seed. We have to take in account
        ## that the walker with restart in some of the Seed nodes, depending on
        ## the score we gave in that file).
        prox_vector <- matrix(0, nrow = NetworkSize, ncol = 1)

        prox_vector[which(base::colnames(x) %in% Seeds_Score[, 1])] <- (Seeds_Score[, 2])

        prox_vector  <- prox_vector / sum(prox_vector)
        restart_vector <-  prox_vector

        while(residue >= Threeshold){

            old_prox_vector <- prox_vector
            prox_vector <- (1 - r) * (x %*% prox_vector) + r * restart_vector
            residue <- sqrt(sum((prox_vector - old_prox_vector) ^ 2))
            iter <- iter + 1;
        }

        NodeNames <- character(length = N)
        Score = numeric(length = N)

        rank_global <- data.frame(NodeNames = NodeNames, Score = Score)
        rank_global$NodeNames <- gsub("_1", "", row.names(prox_vector)[seq_len(N)])

        if (MeanType == "Geometric"){
            rank_global$Score <- geometric.mean(as.vector(prox_vector[, 1]), L, N)
        } else {
            if (MeanType == "Arithmetic") {
                rank_global$Score <- regular.mean(as.vector(prox_vector[, 1]), L, N)
            } else {
                rank_global$Score <- sumValues(as.vector(prox_vector[, 1]), L, N)
            }
        }

        if (DispResults == "TopScores"){
            ## We sort the nodes according to their score.
            Global_results <-
                rank_global[with(rank_global, order(-Score, NodeNames)), ]

            ### We remove the seed nodes from the Ranking and we write the results.
            Global_results <-
                Global_results[which(!Global_results$NodeNames %in% Seeds), ]
        } else {
            Global_results <- rank_global
        }

        base::rownames(Global_results) <- c()

        RWRM_ranking <- list(RWRM_Results = Global_results, Seed_Nodes = Seeds)

        class(RWRM_ranking) <- "RWRM_Results"
        return(RWRM_ranking)
    }



get.seed.scoresMultiplex <- function(Seeds, Number_Layers, tau) {

    Nr_Seeds <- length(Seeds)

    Seeds_Seeds_Scores <- rep(tau / Nr_Seeds, Nr_Seeds)
    Seed_Seeds_Layer_Labeled <-
        paste0(rep(Seeds, Number_Layers), sep = "_", rep(seq(Number_Layers),
                                                    length.out = Nr_Seeds * Number_Layers, each = Nr_Seeds))

    Seeds_Score <- data.frame(Seeds_ID = Seed_Seeds_Layer_Labeled,
                              Score = Seeds_Seeds_Scores, stringsAsFactors = FALSE)

    return(Seeds_Score)
}


geometric.mean <- function(Scores, L, N) {

    FinalScore <- numeric(length = N)

    for (i in seq_len(N)){
        FinalScore[i] <- prod(Scores[seq(from = i, to = N * L, by = N)]) ^ (1 / L)
    }

    return(FinalScore)
}



regular.mean <- function(Scores, L, N) {

    FinalScore <- numeric(length = N)

    for (i in seq_len(N)){
        FinalScore[i] <- mean(Scores[seq(from = i, to = N * L, by = N)])
    }

    return(FinalScore)
}

sumValues <- function(Scores, L, N) {

    FinalScore <- numeric(length = N)

    for (i in seq_len(N)){
        FinalScore[i] <- sum(Scores[seq(from = i, to = N * L, by = N)])
    }

    return(FinalScore)
}
