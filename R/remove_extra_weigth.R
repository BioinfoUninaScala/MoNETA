#' Remove extra weigth
#'
#' @param network a network with three columns
#' @param weight_remv the weight to be removed in the network
#' @return network with weigths strictly less than weight_remv
#' @export

remove_extra_weigth <- function(network, weight_remv = NULL) {

    if (is.null(weight_remv)) {
        weight_remv <- max(network[, "weight"])
    }

    network[network[, "weight"] < weight_remv,]

}
