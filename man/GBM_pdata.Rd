% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/GBM_pdata-data.R
\docType{data}
\name{GBM_pdata}
\alias{GBM_pdata}
\title{Glioblastoma multiforme sample annotation data}
\format{
A DataFrame with 788 rows and 4 columns, containing:
  \describe{
  \item{Case}{Sample ids}
  \item{Histology}{The classification of the GBM tumor based on its microscopic appearance and characteristics}
  \item{IDH.status}{A binary indicator of whether the GBM tumor has an IDH mutation or is wild type}
  \item{Supervised.DNA.Methylation.Cluster}{The categorization of the GBM tumor based on molecular features and methylation patterns}
}
}
\usage{
data(GBM_pdata)
}
\description{
This dataframe is designed to store phenotype data for glioblastoma (GBM) samples.
Specifically, it includes information on the histology, subtype, and IDH mutation status for each GBM sample.
}
\examples{
data(GBM_pdata)
}
\keyword{datasets}
