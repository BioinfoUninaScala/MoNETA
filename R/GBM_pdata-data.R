#' Glioblastoma multiforme sample annotation data
#'
#' This dataframe is designed to store phenotype data for glioblastoma (GBM) samples.
#' Specifically, it includes information on the histology, subtype, and IDH mutation status for each GBM sample.
#'
#' @docType data
#'
#' @usage data(GBM_pdata)
#'
#' @format A DataFrame with 788 rows and 4 columns, containing:
#'   \describe{
#'   \item{Case}{Sample ids}
#'   \item{Histology}{The classification of the GBM tumor based on its microscopic appearance and characteristics}
#'   \item{IDH.status}{A binary indicator of whether the GBM tumor has an IDH mutation or is wild type}
#'   \item{Supervised.DNA.Methylation.Cluster}{The categorization of the GBM tumor based on molecular features and methylation patterns}
#' }
#'
#' @keywords datasets
#'
#' @examples
#' data(GBM_pdata)
"GBM_pdata"
