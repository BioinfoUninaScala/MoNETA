#' Create a multiplex, a DataFrame that comprises different omics with these information: EdgeType, source, target, weight
#'
#' @importFrom dplyr bind_rows select mutate
#' @importFrom purrr map
#' @param omics_list A list of DataFrames containings omics data, each is a DataFrame with three columns: source, destination, weight
#' @return A DataFrame that comprises each DataFrame with these columns: EdgeType, source, target, weight
#' @export


create_multiplex <- function(omics_list) {

    # rbind omics matrices, add id column

    multiplex <- omics_list %>% purrr::map(dplyr::select, source, dest) %>% dplyr::bind_rows(.id = "EdgeType") %>%
        dplyr::select(EdgeType, source = source, target = dest) %>%
        dplyr::mutate(weight = 1)

    multiplex
}
