#' Get intersection between cancer matrices
#'
#' @param omics_list a list of omics
#' @export

get_intersection_matrices <- function(omics_list) {

    col_inter <- base::Reduce(intersect, base::lapply(omics_list, colnames))
    base::lapply(omics_list, function(x) x[,col_inter])
}
