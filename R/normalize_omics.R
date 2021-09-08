#' Normalize matrix using ???
#'
#' @param mtx a omics matrix
#' @return normalized omics matrix
#' @export

normalize_omics <- function(mtx) {

    # using standardNormalization for now, after by our

    #SNFtool::standardNormalization(mtx)
    #@importFrom SNFtool standardNormalization
    scale(mtx)

    # return a normalized matrix

}
