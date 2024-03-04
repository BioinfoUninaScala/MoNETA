<p align="right">
 <img src="https://github.com/BioinfoUninaScala/MoNETA/blob/main/Moneta%20-%20icon.jpg" width="250" alt="MoNETA Logo">
</p>

## MoNETA
#### A new R package for the multi-omics data integration 

## Introduction
MoNETA is a tool designed to **compress multi-omic** data into a single matrix of reduced size. 
Given multiple data representing different omics, it identify relevant multi-omics relationships between samples, merging all the information in a low dimensional space where close samples are more similar.
In this way, you can perform with the output numerous tasks, such as: stratify samples, train ML algorithms, better human-visualization.

---------

### Installation 
In R console, run 

```r
library(devtools)
install_github("BioinfoUninaScala/MoNETA", 
               build_vignettes=FALSE, 
               repos=BiocManager::repositories(),
               dependencies=TRUE, type="source")
```
----------

### MoNETA shiny

The MoNETA pipeline is also implemented as a user-friendly shiny app, available both in MoNETA package and as docker image (https://hub.docker.com/r/bioinfouninascala/moneta). <br/>

#### MoNETA function
In R, you can run:
```r
library(MoNETA)
MoNETAshiny()
```

#### MoNETA docker image
To pull the docker image, run in the terminal:
```
$ docker pull bioinfouninascala/moneta
```

To run the docker image, map http port 3838 on the host port 8585:
```
$ docker run --rm -p 8585:3838 bioinfouninascala/moneta
```

Finally, use MoNETA shiny in a browser by typing http://localhost:8585 (if you are on a local machine) or https://server_name:8585 (if you are on a server machine).
<br/>

### Usage 

We provide an example of the workflow for our project. To access it, follow the link : [MoNETA use case](https://bioinfouninascala.github.io/MoNETA/inst/doc/MoNETA_use_case.html).

#### Main functions

* `normalize_omics`
* `get_intersection_matrices`
* `remove_zeros_cols`
* `prune_multiplex_network`

* `k_star_net`
* `create_multiplex`
* `create_layer_transition_matrix`
* `gen_sim_mat_M`

* `get_embedding`
* `get_pca_embedding`
* `get_tsne_embedding`
* `get_umap_embedding`
* `get_parallel_umap_embedding`

* `plot_net`
* `plot_2D_matrix`

#### Other functions

* `MoNETAshiny`

### Workflow 

The functions that are provided by the tool can be grouped into 4 modules: 



**1. Preprocessing**
Input data consist of the matrices with samples on columns

* `normalize_omics` : takes in input a Matrix that has to be centered and scaled along the columns

* `get_intersection_matrices` **(NOT MANDATORY)** : get intersection of columns (samples) of cancer matrices, returns a list of matrices with samples that appear in all of them

* `remove_zeros_cols` : cleans a Matrix removing the columns with all zeros

* `prune_multiplex_network` : returns network with weigths strictly less than a given threshold



**2. Generation of similarity matrix using multi-omics data**

First apply `k_star_net` function to each omics data. 
After with `create_multiplex` and `create_layer_transition_matrix` the networks will be combined. 
Finally use `gen_sim_mat_M` to generate the similarity matrix that describes the multi-omics integration.

* `k_star_net` : applies [K\*nn greedy algorithm](https://papers.nips.cc/paper/2016/file/2c6ae45a3e88aee548c0714fad7f8269-Paper.pdf) to dynamically find nearest neighbors for each node 

* `create_multiplex` : creates a multiplex, which corresponds to a DataFrame that comprises different omics with these information: EdgeType, source, target, weight

* `create_layer_transition_matrix` : creates jump matrix containing the probabilities to switch from one omics to another

* `gen_sim_mat_M` : generates similarity matrix using a random walker with restart



**3. Dimensionality reduction**

It is possible to choose different kind of dimensionality reduction algorithm.

* `get_embedding` : computes an embedding of the similarity matrix with [MultiVERSE algorithm](https://github.com/Lpiol/MultiVERSE) described by [Léo Pio-Lopez, et al.](https://arxiv.org/abs/2008.10085)

* `get_pca_embedding` : computes an embedding of the similarity matrix using Principal Component Analysis (PCA)

* `get_tsne_embedding` : computes an embedding of the similarity matrix using t-distributed stochastic neighbor embedding (t-SNE)

* `get_umap_embedding` : computes an embedding of the similarity matrix using Uniform Manifold Approximation and Projection for Dimension Reduction (UMAP)

* `get_parallel_umap_embedding` : computes an embedding of the similarity matrix using a parallelized version of UMAP



**4. Plots**

* `plot_net` : plot the networks created by `k_star_net`

* `plot_2D_matrix` : plot the a 2D matrix with two rows representing the coordinates

