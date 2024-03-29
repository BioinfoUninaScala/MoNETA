## ----install, eval=F, echo=T---------------------------------------------------------------------------------------------------------------------------
#  
#  devtools::install_github("BioinfoUninaScala/MoNETA")
#  

## ----read-files----------------------------------------------------------------------------------------------------------------------------------------

library(MoNETA)

data(GBM_pdata)
data(GBM_mtx)



## ----dim-data------------------------------------------------------------------------------------------------------------------------------------------
message("Multi-omic \t # of patients \t # of features",
        "\nCNV\t\t\t\t", ncol(GBM_mtx$GliomaCNV_norm), "\t\t\t\t", nrow(GBM_mtx$GliomaCNV_norm),
        "\nMethylation\t\t", ncol(GBM_mtx$GliomaMethylation_norm), "\t\t\t\t", nrow(GBM_mtx$GliomaMethylation_norm),
        "\nExpression\t\t", ncol(GBM_mtx$GliomaExpression_norm), "\t\t\t\t", nrow(GBM_mtx$GliomaExpression_norm))

## ----data-show-----------------------------------------------------------------------------------------------------------------------------------------

GBM_mtx$GliomaCNV_norm[1:5,1:5]


## ----preprocessing1------------------------------------------------------------------------------------------------------------------------------------

GBM_mtx$GliomaExpression_norm <- normalize_omics(GBM_mtx$GliomaExpression_norm)
GBM_mtx$GliomaCNV_norm <- normalize_omics(GBM_mtx$GliomaCNV_norm)
GBM_mtx$GliomaMethylation_norm <- normalize_omics(GBM_mtx$GliomaMethylation_norm)


## ----k-star-net----------------------------------------------------------------------------------------------------------------------------------------

net_list <- list(CNV_norm = k_star_net(matrix = GBM_mtx$GliomaCNV_norm,
                                       sparsity = .7,
                                       distFun = "Euclidean", cores = 50),
                 METH_norm = k_star_net(matrix = GBM_mtx$GliomaMethylation_norm,
                                        sparsity = .7,
                                        distFun = "Euclidean", cores = 50),
                 Expr_norm = k_star_net(matrix = GBM_mtx$GliomaExpression_norm,
                                        sparsity = .7,
                                        distFun = "Euclidean", cores = 50)
)


## ----show-k-star-net-----------------------------------------------------------------------------------------------------------------------------------

print(head(net_list$CNV_norm))


## ----multiplex-----------------------------------------------------------------------------------------------------------------------------------------

multiplex <-  create_multiplex(net_list, weighted = T)


## ----show-multiplex------------------------------------------------------------------------------------------------------------------------------------

print(head(multiplex))


## ----prune-multiplex-----------------------------------------------------------------------------------------------------------------------------------

multiplex <-  prune_multiplex_network(multiplex, 1000)


## ----layer-transition----------------------------------------------------------------------------------------------------------------------------------


layer_transition  <-  create_layer_transition_matrix(net_list)


## ----show-layer-transition-----------------------------------------------------------------------------------------------------------------------------


print(layer_transition)


## ----similarity-matrix---------------------------------------------------------------------------------------------------------------------------------

RWR_mat   <-  gen_sim_mat_M(network = multiplex,
                         tau = NA, restart = 0.7,
                         jump_neighborhood = F, weighted_multiplex = F, cores = 50)


## ----show-similarity-matrix----------------------------------------------------------------------------------------------------------------------------

print(RWR_mat[1:5,1:5])


## ----embedding-----------------------------------------------------------------------------------------------------------------------------------------

tsne_emb = get_tsne_embedding(RWR_mat, 2, 70, 20000, 50)

multiverse_emb       = get_embedding(RWR_mat, 64, num_steps = 10 ^ 5, cores = 50)

pca_emb   = get_pca_embedding(RWR_mat, 2)

umap_emb  = get_umap_embedding(RWR_mat, 2)

pumap_emb  = get_parallel_umap_embedding(RWR_mat, 2, n_neighbors = 15, n_threads = 10, n_sgd_threads = 0, grain_size = 1)


## ----show-emb------------------------------------------------------------------------------------------------------------------------------------------

print(multiverse_emb[1:5,1:5])


## ----plot-net, out.width="100%"------------------------------------------------------------------------------------------------------------------------

plot_net(edgeList = net_list$METH_norm, nodes_anno = GBM_pdata,
         interactive = FALSE, id_name = "Case",
         id_anno_color = "Supervised.DNA.Methylation.Cluster",
         id_anno_shape = "IDH.status", html = FALSE, wo_legend = FALSE, title = "METH_norm")


## ----plot-2d-matrix, fig.width = 7, fig.height = 4-----------------------------------------------------------------------------------------------------


plot_2D_matrix(coord = tsne_emb, nodes_anno = GBM_pdata, id_name = "Case", interactive = FALSE,
                   id_anno_color = "Supervised.DNA.Methylation.Cluster", id_anno_shape = "IDH.status",
                   wo_legend = FALSE, title = "t-SNE embedding")


## ----plot-2d-matrix2, fig.width = 7, fig.height = 4----------------------------------------------------------------------------------------------------

plot_2D_matrix(coord = pca_emb, nodes_anno = GBM_pdata, id_name = "Case", interactive = FALSE,
                   id_anno_color = "Supervised.DNA.Methylation.Cluster", id_anno_shape = "IDH.status",
                   wo_legend = FALSE, title = "First two PC")


