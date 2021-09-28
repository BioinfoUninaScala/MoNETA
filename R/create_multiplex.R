#' Create a multiplex, an object that comprises different omics
#'
#' @importFrom dplyr bind_rows select mutate
#' @importFrom purrr map
#' @param omics_list A list of DataFrames containings omics data, each with three columns : source, destination, weight
#' @return An object that comprises each DataFrame
#' @export


create_multiplex <- function(omics_list) {

    # rbind omics matrices, add id column

    multiplex <- omics_list %>% purrr::map(dplyr::select, source, dest) %>% dplyr::bind_rows(.id = "EdgeType") %>%
        dplyr::select(EdgeType, source = source, target = dest) %>%
        dplyr::mutate(weight = 1)

    multiplex
}
