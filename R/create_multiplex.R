#' Create multiplex
#'
#' @importFrom dplyr bind_rows select mutate
#' @importFrom purrr map
#' @param omics_list a list of omics
#' @return multiplex
#' @export


create_multiplex <- function(omics_list) {

    # rbind omics matrices, add id column

    multiplex <- omics_list %>% purrr::map(dplyr::select, source, dest) %>% dplyr::bind_rows(.id = "EdgeType") %>%
        dplyr::select(EdgeType, source = source, target = dest) %>%
        dplyr::mutate(weight = 1)

    multiplex
}
