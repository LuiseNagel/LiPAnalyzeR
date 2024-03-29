% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/PlottingResults.R
\name{makeWoodsPlotSingleProtein}
\alias{makeWoodsPlotSingleProtein}
\title{Creating wood plots over a single proteins}
\usage{
makeWoodsPlotSingleProtein(sumDf, annotPP, coefCol="Coefficient",
pvalCol="Padj", nameProtQuant="Protein", startPosition="startPosition",
endPosition="endPosition", pvalCutoff=0.05, deltaColorIsTryptic=FALSE,
isTryptic="isTryptic", nameFT=c("Specific"),
nameHT=c("Specific-C", "Specific-N"), xlim=NULL, ylim=NULL, protName=NULL,
showPv=FALSE)
}
\arguments{
\item{sumDf}{A data.frame containing coefficients and p-values of the
condition of interest were rows corresponding to the quantity of interest
(e.g. peptides). Can be the output of \code{summarizeModelResults}.}

\item{annotPP}{A data.frame with peptides (/modified peptides/precursors) and
protein annotation. Rows are features and must contain the row names of
\code{sumDf}. Can be output of \code{getPepProtAnnot}. Must include all
columns needed for plotting the data:
\itemize{
  \item \code{nameProtQuant} protein IDs
  \item \code{startPosition} AA positions were peptides (/modified peptides
  /precursors) start in protein
  \item \code{endPosition} AA positions were peptides (/modified peptides
  /precursors) end in protein
  \item \code{isTryptic} optional column providing information if peptide is
  full-tryptic or half-tryptic, only needed if
  \code{deltaColorIsTryptic = 'TRUE'}
}}

\item{coefCol}{A character vector or numeric value defining column of
\code{sumDf} were coefficients of the quantity of interest (e.g. peptides)
are provided.
Default is 'Coefficient'.}

\item{pvalCol}{A character vector or numeric value defining column of
\code{sumDf} were (adjusted) p-values  of the quantity of interest
(e.g. peptides) are provided.
Default is 'Padj'.}

\item{nameProtQuant}{A character string giving column of \code{annotPP} which
protein names are provided. If a peptides matches to the same protein several
times, the protein name should be provided each time, separated by ','. If a
peptide maps to the multiple proteins, these different proteins can be
provided by separating them with ';'.
Default is 'Protein'.}

\item{startPosition}{A character string or numeric giving the column name or
column number in which start position from each peptide in its protein
sequence are provided in \code{annotPP}. If a peptides matches to the sam
protein several times, different start positions should be separated by ','.
If a peptide maps to the multiple proteins, the start positions should be
separated with ';'.
Default is 'startPosition'.}

\item{endPosition}{A character string or numeric giving the column name or
column number in which end position from each peptide in its protein
sequence are provided in \code{annotPP}. If a peptides matches to the sam
protein several times, different start positions should be separated by ','.
If a peptide maps to the multiple proteins, the start positions should be
separated with ';'.
Default is 'endPosition'.}

\item{pvalCutoff}{A numeric value, peptides with p-values below this cut-off
are considered significant.
Default is '0.05'.}

\item{deltaColorIsTryptic}{A boolean variable, if set to 'TRUE', full-tryptic
and half-tryptic peptides will be plotted in different colors, this requires
\code{isTryptic} column in \code{annotPP}
Default is FALSE'.}

\item{isTryptic}{A character string or numeric giving the column name or
column number of \code{annotPP} in which digest type can be found
Default is 'isTryptic'.}

\item{nameFT}{A character vector defining the name of full-tryptic peptides
provided in the \code{isTryptic} column in \code{annotPP}.
Default is c("Specific").}

\item{nameHT}{A character vector defining the name of half-tryptic peptides
provided in the \code{isTryptic} column in \code{annotPP}.
Default is c("Specific-C", "Specific-N").}

\item{xlim}{A numeric vector of the length two defining the limits of the
x-axis of the plot. If set to NULL, limits of x-axis are chosen based on
sequence positions of the peptides.
Default is NULL.}

\item{ylim}{A numeric vector of the length two defining the limits of the
y-axis of the plot. If set to NULL, limits of x-axis are chosen based on
the coefficients provided in \code{sumDf}.
 Default is NULL.}

\item{protName}{A character vector giving name of protein to be plotted. Has
to be included in the column of \code{nameProtQuant} of \code{annotPP}.}

\item{showPv}{A boolean value, if set to 'TRUE' p-values as defined in
\code{pvalCol} should be displayed in plots.
Default is 'FALSE'.}
}
\description{
Creating wood plots over a single proteins
}
