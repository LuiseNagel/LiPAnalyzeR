% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/RunningModels.R
\name{createModelMatrix}
\alias{createModelMatrix}
\title{Creating model matrices to run RUV or contrast models}
\usage{
createModelMatrix(quantityList, formula, annotS, samples)
}
\arguments{
\item{quantityList}{A list of preprocessed matrices, containing quantities of
interest(e.g. peptide, modified peptide, precursor) and protein abundances.
Rows represent features and columns samples and should match between the
different matrices contained in the list.
Output from \code{preprocessQuantityMatrix} have to use the variable naming
'Y', 'XPep' and 'XProt'.}

\item{formula}{A formula providing structure of model matrices created in
this function}

\item{annotS}{A data.frame containing sample annotation. Must contain all
columns included in the RUV and contrast models. Rows are samples and must
match to columns of the matrices in \code{quantityList}. Must include
columns of any further variables used in \code{formulaRUV}.}

\item{samples}{A character vector providing sample names for the model.}
}
\value{
A list with model matrices for running the RUV or contrast models.
}
\description{
Creates one model matrices per peptide/protein which are
used in the RUV or contrast model function afterwards
}
