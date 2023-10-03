#' Glioblastoma multiforme omics data
#'
#' GBM_mtx is a list of 3 matrices containing CNV, DNA methylation and RNA expression molecular data for 788 glioblastoma multiforme (GBM) samples.
#' Specifically, each matrix in the list contains features on rows and samples on columns and represents a different omics assy.
#'
#' @docType data
#'
#' @usage data(GBM_mtx)
#'
#' @format The matrices included in GBM_mtx dataset are:
#'   \describe{
#'   \item{GliomaCNV_norm}{A matrix containing normalized copy number variations (CNV) data for 788 GBM samples, with  24776 genes as features on rows}
#'   \item{GliomaMethylation_norm}{A matrix containing normalized DNA methylation data for 788 GBM samples, with 1300 individual CpG sites as features on rows}
#'   \item{GliomaExpression_norm}{A matrix containing normalized gene expression data for 788 GBM samples, with 12985 genes as features on rows}
#' }
#'
#' @keywords datasets
#'
#' @examples
#' data(GBM_mtx)
"GBM_mtx"
