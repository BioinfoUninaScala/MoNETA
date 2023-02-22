## ----install, eval=F, echo=T----------------------------------------------------------------------
#  
#  devtools::install_github("BioinfoUninaScala/MoNETA")
#  

## ----read-files-----------------------------------------------------------------------------------

library(MoNETA)

pdata <- readRDS("pdata.RDS")
GBM_mtx <- readRDS("GBM_mtx.RDS")


## ----dim-data-------------------------------------------------------------------------------------
message("Multi-omic \t # of patients \t # of features",
        "\nCNV\t\t\t\t", ncol(GBM_mtx$GliomaCNV_norm), "\t\t\t\t", nrow(GBM_mtx$GliomaCNV_norm),
        "\nMethylation\t\t", ncol(GBM_mtx$GliomaMethylation_norm), "\t\t\t\t", nrow(GBM_mtx$GliomaMethylation_norm),
        "\nExpression\t\t", ncol(GBM_mtx$GliomaExpression_norm), "\t\t\t\t", nrow(GBM_mtx$GliomaExpression_norm),
        "\nMutation\t\t\t", ncol(GBM_mtx$GliomaMUT_norm), "\t\t\t\t", nrow(GBM_mtx$GliomaMUT_norm))

## ----data-show------------------------------------------------------------------------------------

GBM_mtx$GliomaCNV_norm[1:5,1:5]


## ----preprocessing1-------------------------------------------------------------------------------

GBM_mtx$GliomaExpression_norm <- normalize_omics(GBM_mtx$GliomaExpression_norm)
GBM_mtx$GliomaCNV_norm <- normalize_omics(GBM_mtx$GliomaCNV_norm)
GBM_mtx$GliomaMethylation_norm <- normalize_omics(GBM_mtx$GliomaMethylation_norm)


## ----preprocessing2-------------------------------------------------------------------------------

GBM_mtx$GliomaMUT_norm = remove_zeros_cols(GBM_mtx$GliomaMUT_norm)


## ----preprocessing3-------------------------------------------------------------------------------

GBM_mtx <- get_intersection_matrices(GBM_mtx)


## ----dim-data2------------------------------------------------------------------------------------
message("Multi-omic \t # of patients \t # of features",
        "\nCNV\t\t\t\t", ncol(GBM_mtx$GliomaCNV_norm), "\t\t\t\t", nrow(GBM_mtx$GliomaCNV_norm),
        "\nMethylation\t\t", ncol(GBM_mtx$GliomaMethylation_norm), "\t\t\t\t", nrow(GBM_mtx$GliomaMethylation_norm),
        "\nExpression\t\t", ncol(GBM_mtx$GliomaExpression_norm), "\t\t\t\t", nrow(GBM_mtx$GliomaExpression_norm),
        "\nMutation\t\t\t", ncol(GBM_mtx$GliomaMUT_norm), "\t\t\t\t", nrow(GBM_mtx$GliomaMUT_norm))

## ----k-star-net-----------------------------------------------------------------------------------

net_list <- list(CNV_norm = k_star_net(matrix = GBM_mtx$GliomaCNV_norm,
                                       sparsity = .7,
                                       distFun = "Euclidean", cores = 50),
                 METH_norm = k_star_net(matrix = GBM_mtx$GliomaMethylation_norm,
                                        sparsity = .7,
                                        distFun = "Euclidean", cores = 50),
                 Expr_norm = k_star_net(matrix = GBM_mtx$GliomaExpression_norm,
                                        sparsity = .7,
                                        distFun = "Euclidean", cores = 50),
                 MUT_norm = k_star_net(matrix = GBM_mtx$GliomaMUT_norm,
                                       sparsity = .1,
                                       distFun = "Manhattan", knn = 10, cores = 50)
)


## ----show-k-star-net------------------------------------------------------------------------------

print(head(net_list$CNV_norm))


## ----multiplex------------------------------------------------------------------------------------

multiplex <-  create_multiplex(net_list, weighted = T)


## ----show-multiplex-------------------------------------------------------------------------------

print(head(multiplex))


## ----prune-multiplex------------------------------------------------------------------------------

multiplex <-  prune_multiplex_network(multiplex, 1000)


## ----layer-transition-----------------------------------------------------------------------------


layer_transition  <-  create_layer_transition_matrix(net_list)


## ----show-layer-transition------------------------------------------------------------------------


print(layer_transition)


## ----similarity-matrix----------------------------------------------------------------------------

RWR_mat   <-  gen_sim_mat_M(network = multiplex,
                         tau = NA, restart = 0.7,
                         jump_neighborhood = F, weighted_multiplex = F, cores = 50)


## ----show-similarity-matrix-----------------------------------------------------------------------

print(RWR_mat[1:5,1:5])


## ----embedding------------------------------------------------------------------------------------

multiverse_emb       = get_embedding(RWR_mat, 64, num_steps = 10 ^ 5, cores = 50)

pca_emb   = get_pca_embedding(RWR_mat, 2)

umap_emb  = get_umap_embedding(RWR_mat, 2)

pumap_emb  = get_parallel_umap_embedding(RWR_mat, 2, n_neighbors = 15, n_threads = 10, n_sgd_threads = 0, grain_size = 1)


## ----show-emb-------------------------------------------------------------------------------------

print(multiverse_emb[1:5,1:5])


## ----plot-net, out.width="100%"-------------------------------------------------------------------

plot_net(edgeList = net_list$MUT_norm, nodes_anno = pdata,
         interactive = FALSE, id_name = "Case",
         id_anno_color = "Supervised.DNA.Methylation.Cluster",
         id_anno_shape = "IDH.status", html = FALSE, wo_legend = FALSE, title = "MUT_norm")


## ----plot-2d-matrix, fig.width = 7, fig.height = 4------------------------------------------------

multiverse_pca_emb   = get_pca_embedding(multiverse_emb, 2)

plot_2D_matrix(coord = multiverse_pca_emb, nodes_anno = pdata, id_name = "Case", interactive = FALSE,
                   id_anno_color = "Supervised.DNA.Methylation.Cluster", id_anno_shape = "IDH.status",
                   wo_legend = FALSE, title = "Embedding MultiVERSE with first two PC")


## ----plot-2d-matrix2, fig.width = 7, fig.height = 4-----------------------------------------------

plot_2D_matrix(coord = pca_emb, nodes_anno = pdata, id_name = "Case", interactive = FALSE,
                   id_anno_color = "Supervised.DNA.Methylation.Cluster", id_anno_shape = "IDH.status",
                   wo_legend = FALSE, title = "First two PC")


