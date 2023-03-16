#' Phenotype data of GBM use case
#'
#' GBM_mtx is a list of matrices containing molecular data on glioblastoma (GBM) samples.
#' Specifically, each matrix in the list contains features on rows and samples on columns,
#' and represents different molecular data types.
#'
#' @docType data
#'
#' @usage data(GBM_mtx)
#'
#' @format The matrices included in GBM_mtx are:
#'   \describe{
#'   \item{GliomaCNV_norm}{A matrix representing copy number variations (CNV) data for GBM samples, with genes as features on rows}
#'   \item{GliomaMethylation_norm}{A matrix representing DNA methylation data for GBM samples, with individual CpG sites as features on rows}
#'   \item{GliomaExpression_norm}{A matrix representing gene expression data for GBM samples, with genes as features on rows}
#' }
#'
#' @keywords datasets
#'
#' @examples
#' data(GBM_mtx)
"GBM_mtx"
