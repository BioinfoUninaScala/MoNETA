<p align="right">
 <img src="https://github.com/BioinfoUninaScala/MoNETA/blob/main/Moneta%20-%20icon.jpg" width="250" alt="EpiStatProfiler Logo">
</p>

## MoNETA
#### A new R package for the multi-omics data integration 

## About
MoNETA is a tool designed to **compress multi-omic** data into a single matrix of reduced size. 

The main problem with modern biological data, representing a single patient, is both the **size and the number of data types**.
So the aim of this project is to identify relevant multi-omics relationships between samples and to represent it in a low dimensional space.
The shape of this undersized matrix is useful for numerous tasks, such as: stratify samples, train ML algorithms, better human-visualization.


## Functions
***normalize_omics*** : takes in input a Matrix that has to be centered and scaled along the columns

***get_intersection_matrices*** **(NOT MANDATORY)** : get intersection of columns (samples) of cancer matrices, returns a list of matrices with samples that appear in all of them


***k_star_net*** : Applies [K\*nn greedy algorithm](https://papers.nips.cc/paper/2016/file/2c6ae45a3e88aee548c0714fad7f8269-Paper.pdf) to dynamically find nearest neighbors for each node  

***create_multiplex*** : creates a multiplex, a DataFrame that comprises different omics with these information: EdgeType, source, target, weight

***create_jump_matrix*** : creates jump matrix containing the probabilities to switch from one omics to another

***gen_sim_mat_M*** : creates a “similarity” matrix of using Random Walk with restart from each node applied to the multiplex.


***get_embedding*** : compute an embedding of the multi-omics matrix using the [MultiVERSE algorithm](https://github.com/Lpiol/MultiVERSE)


## Plotting Functions

***plot_umap***
***plot_umap_db***
***plot_net***
