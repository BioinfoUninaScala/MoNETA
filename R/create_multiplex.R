#' Create a multiplex omic network
#'
#' @importFrom dplyr bind_rows select mutate
#' @importFrom purrr map
#' @param omics_list a named list of data frames containing adjacency lists for different omics networks sharing node ids.
#'  Each data frame in the list is required to have at least two columns containing nodes ids respectively named source and destination and one column named weight containing edge weights.
#' @param weighted logical. If TRUE then a weighted multiplex is returned, otherwise weights are all set to 1. FALSE by default.
#' @return A data frame representing the edge list of a multiplex network with the following four columns: EdgeType , source, target, weight.
#' EdgeType is a character representing the layers in the multiplex, source and target represent the edges and contain node ids, weight contains edge weights representing node similarity in the multiplex network.
#' @export


create_multiplex <- function(omics_list, weighted = FALSE) {

    # rbind omics matrices, add id column

    multiplex <- omics_list %>% purrr::map(dplyr::select, source, dest, weight) %>% dplyr::bind_rows(.id = "EdgeType") %>%
        dplyr::select(EdgeType, source = source, target = dest, weight = weight)
    if (! weighted) {
        multiplex <- multiplex %>% dplyr::mutate(weight = 1)
    }

    multiplex
}
