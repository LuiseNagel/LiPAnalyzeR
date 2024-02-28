#' @title Extracting model results and performing FDR correction on p-values.
#'
#' @description Function for extracting p-values and coefficients of interest
#' from the model output (\code{runModel}/\code{analyzeLiPPepData}/
#' \code{analyzeTrPPepData}/\code{analyzeTrPProtData}) output. The function
#' additionally performs FDR correction, either over all proteins or
#' protein-wise (meaning over all peptides matching to the same protein).
#'
#' @usage summarizeModelResults(resModel, evalCovariable="Condition_Contrast",
#' correctPval="all", pAdjust="fdr", annotPP=NULL, nameIDQuant="quantID",
#' nameProtQuant="Protein")
#'
#' @param resModel A list object, output from \code{runModel}\
#' \code{AnalyzeLiPPepData}\\code{AnalyzeTrPPepData}\\code{AnalyzeTrPProtData}.
#' Can only be run if the contrast model was run in the function.
#' @param evalCovariable A character string providing name of columns in which
#' coefficients and p-values of the variable of interest can be found.
#' Default is 'Condition_Contrast'.
#' @param correctPval A character string defining approach to adjust p-values.
#' \itemize{
#'   \item 'all': All p-values are FDR adjusted together.
#'   \item 'protein-wise': FDR adjustment is performed in smaller groups,
#'   meaning over all peptides matching to the same protein
#'   }
#' @param pAdjust A character string giving the correction method to use when
#' performing FDR Check '?p.adjust' for more information.
#' Default is 'fdr' (Benjamini & Hochberg (1995)).
#' @param annotPP A data.frame with peptides (/modified peptides/precursors) and
#' protein annotation that must be provided if \code{correctPval} = 'all'. Rows
#' are features. Data.frame must contain columns providing the
#' \code{nameIDQuant} and matching protein names in \code{nameProtQuant}.
#' The output from \code{getPepProtAnnot} can be given here.
#' Default is 'NULL'.
#' @param nameIDQuant A character string giving column name of \code{annotPP}
#' where quantity IDs (e.g. peptide identifies/sequences) used in row names of
#' \code{resModel}.
#' @param nameProtQuant A character string giving column of \code{annotPP} which
#' protein names/protein groups are provided. If \code{correctPval} is set to
#' 'protein-wise', peptides with the same \code{nameProtQuant} value will be
#' corrected together.
#' Default is 'Protein'.
#'
#' @return A data.frame were rows are features and columns are coefficient,
#' p-value and adjusted p-value of the condition of interest.
#'
#' @export
summarizeModelResults <- function(resModel, evalCovariable="Condition_Contrast",
                                  correctPval="all", pAdjust="fdr",
                                  annotPP=NULL, nameIDQuant="quantID",
                                  nameProtQuant="Protein"){
    modelCoef <- stats::setNames(resModel[[1]][, evalCovariable],
                                 row.names(resModel[[1]]))
    modelPval <- stats::setNames(resModel[[2]][, evalCovariable],
                                 row.names(resModel[[2]]))
    if(correctPval == "all"){
        modelPvalAdj <- stats::p.adjust(modelPval, method=pAdjust)
    }
    else if(correctPval == "protein-wise"){
        message("Preforming protein-wise FDR correction.")
        modelPvalAdj <- doProteinWiseFDR(modelPval, pAdjust, annotPP,
                                         nameIDQuant, nameProtQuant)
    }
    else{
        stop("No valid option choosen for correctPval. Please set correctPval to
             'all' or 'protein-wise'.")
    }
    modelSum <- data.frame(Coefficient=modelCoef[order(modelPvalAdj,
                                                       decreasing=FALSE)],
                           Pvalue=modelPval[order(modelPvalAdj,
                                                  decreasing=FALSE)],
                           Padj=modelPvalAdj[order(modelPvalAdj,
                                                   decreasing=FALSE)])
    return(modelSum)
}

#' @title Performing protein-wise FDR on peptide p-values
#' @param vecPval A named numeric vector with p-values, were names are the
#' corresponding features
#' @param pAdjust A character string giving the correction method to use when
#' performing FDR Check ?p.adjust for more information. Set to 'fdr' (Benjamini
#' & Hochberg (1995)) by default.
#' @param annotPP A data.frame with peptides (/modified peptides/precursors) and
#' protein annotation. Rows are features. Data.frame must contain
#' columns providing the \code{nameIDQuant} and matching protein names in
#' \code{nameProtQuant}. Per default, the output from \code{getPepProtAnnot}
#' can be given here.
#' @param nameIDQuant A character string giving column name of \code{annotPP}
#' where quantity IDs (e.g. peptide identifies/sequences) used in names of
#' \code{vecPval}.
#' @param nameProtQuant A character string giving column of \code{annotPP} which
#' protein names/protein groups are provided. If \code{correctPval} is set to
#' 'protein-wise', peptides with the same \code{nameProtQuant} value will be
#' corrected together.
#' Default is 'Protein'.
#'
#' @return vector of protein-wise adjusted P values

doProteinWiseFDR <- function(vecPval, pAdjust, annotPP, nameIDQuant,
                             nameProtQuant){
    if(is.null(annotPP)){
        stop("Please provide peptide and protein annotations in form of the
             annotPP matrix.\n")
    }
    if(length(intersect(names(vecPval), annotPP[,nameIDQuant])) == 0){
        stop("Please provide correct level on which the models were build. This
             level must be a column name in the annotPP file.\n")
    }
    matchPepProt <- annotPP[match(names(vecPval), annotPP[,nameIDQuant]),]
    listPval <- split(vecPval, matchPepProt[, nameProtQuant])

    vecAdj <- unlist(lapply(listPval, function(x){
        x <- stats::p.adjust(x, method=pAdjust)
        return(x)
    }), use.names=FALSE)
    names(vecAdj) <- unlist(lapply(listPval, function(x){
        names(x)
    }))
    vecAdj <- vecAdj[names(vecPval)]
    return(vecAdj)
}

