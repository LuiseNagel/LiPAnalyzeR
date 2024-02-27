globalVariables(names=c("Pval", "Yaxis", "aes",  "start", "end"))


#' @title Extracting model results and performing FDR correction on p-values.
#'
#' @description Function for extracting p-values and coefficients from
#' \code{RunModel} output. Also performs FDR correction, either over all
#' peptides/proteins or over all peptides of each protein individually
#'
#' @usage summarizeModelResults(resModel, evalCovariable="Condition_Contrast",
#' correctPval="all", pAdjust="fdr", annotPP=NULL, infoFeature="quantID",
#' infoProtein="Protein")
#'
#' @param resModel A list object, output from \code{runModel} or
#' \code{AnalyzeLiPPepData}\\code{AnalyzeTrpPepData}\\code{AnalyzeTrpProtData}.
#' @param evalCovariable A character string giving name of columns in which
#' coefficients and p-values of the variable of interest can be found. Default
#' is set to 'Condition_Contrast'.
#' @param correctPval A character string defining how p-values should be
#' adjusted. Per default all p-values are corrected together, and it is set to
#' 'all'. Alternatively, FDR can be performed protein-wise, please set the
#' variable to 'protein-wise' for this.
#' @param pAdjust A character string giving the correction method to use when
#' performing FDR Check ?p.adjust for more information. Set to 'fdr' (Benjamini
#' & Hochberg (1995)) by default.
#' @param annotPP A data.frame with peptide and protein annotation. Rows are
#' features and columns must include \code{infoFeature} and \code{infoProtein}.
#' Is only required if \code{correctPval} is set to 'protein-wise'.
#' @param infoFeature A character string providing column name of \code{annotPP}
#' where features which were fitted in the models are noted. Default is set to
#' 'Peptide'.
#' @param infoProtein A character string providing column name of \code{annotPP}
#' where protein names/groups are noted. If \code{correctPval} is defines as
#' 'protein-wise' peptides with the same \code{infoProtein} value will be
#' corrected together. Default is set to 'Protein'.
#'
#' @return returns a data.frame including all features of the model and
#' the corresponding coefficients, P values and adjusted P values.
#'
#' @export
#'
summarizeModelResults <- function(resModel, evalCovariable="Condition_Contrast",
                                  correctPval="all", pAdjust="fdr",
                                  annotPP=NULL, infoFeature="quantID",
                                  infoProtein="Protein"){
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
                                         infoFeature, infoProtein)
    }
    else{
        stop("No valid option choosen for correctPval. Please set correctPval to
             'all' or 'protein-wise'.\n")
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
#' @param vecPval A named numeric vector with p-values and corresponding
#' features
#' @param pAdjust A character string giving the correction method to use when
#' performing FDR Check ?p.adjust for more information. Set to 'fdr' (Benjamini
#' & Hochberg (1995)) by default.
#' @param annotPP A data.frame with peptide and protein annotation. Rows are
#' features and columns must include \code{infoFeature} and \code{infoProtein}.
#' Is only required if \code{correctPval} is set to 'protein-wise'.
#' @param infoFeature A character string providing column name of \code{annotPP}
#' where features which were fitted in the models are noted. Default is set to
#' 'Peptide'.
#' @param infoProtein A character string providing column name of \code{annotPP}
#' where protein names/groups are noted. If \code{correctPval} is defines as
#' 'protein-wise' peptides with the same \code{infoProtein} value will be
#' corrected together. Default is set to 'Protein'.
#'
#' @return vector of protein-wise adjusted P values

doProteinWiseFDR <- function(vecPval, pAdjust, annotPP, infoFeature,
                             infoProtein){
    if(is.null(annotPP)){
        stop("Please provide peptide and protein annotations in form of the
             annotPP matrix.\n")
    }
    if(length(intersect(names(vecPval), annotPP[,infoFeature])) == 0){
        stop("Please provide correct level on which the models were build. This
             level must be a column name in the annotPP file.\n")
    }
    matchPepProt <- annotPP[match(names(vecPval), annotPP[,infoFeature]),]
    listPval <- split(vecPval, matchPepProt[, infoProtein])

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

