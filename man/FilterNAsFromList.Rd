% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/FilterSpectroOut.R
\name{FilterNAsFromList}
\alias{FilterNAsFromList}
\title{Filtering peptide and protein quantities based on NAs}
\usage{
FilterNAsFromList(SpectroList, annotS, infoCondition, maxNAperCondition)
}
\arguments{
\item{SpectroList}{A list of matrices, containing peptides/proteins
quantities. Rows represent features and columns refer to the samples. Can be
output from \code{ExtractDataFromSpectro}.}

\item{annotS}{A data.frame containing sample and condition annotation.
Rows are samples and must match to columns of SpectroList/QuantityMatrix.
Needs to contain column about different conditions/groups.}

\item{infoCondition}{A character string providing column name of
\code{annotS} in which condition is provided. Default is 'Condition'. If NA
filtering should not be applied on condition level provide only one value in
\code{infoCondition} in \code{annotS}.}

\item{maxNAperCondition}{A numeric value, defining maximal number of NAs
in an individual feature per condition. Default is '0'.}
}
\value{
list of NA filtered matrices with petide/protein quantities.
}
\description{
Function for removing peptides with too many NAs on one of the
peptide and protein quantity levels.
}
