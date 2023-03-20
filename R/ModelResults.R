globalVariables(names=c("Pval", "Yaxis", "aes",  "start", "end"))


#' @title Extracting model results and performing FDR correction on p-values.
#'
#' @description Function for extracting p-values and coefficients from
#' \code{RunModel} output. Also performs FDR correction, either over all
#' peptides/proteins or over all peptides of each protein individually
#'
#' @usage summarizeModelResults(resModel, evalCovariable="Condition_OLS",
#' correctPval="all", pAdjust="fdr", annotPP=NULL, infoFeature="Peptide",
#' infoProtein="Protein")
#'
#' @param resModel A list object, output from \code{runModel} or
#' \code{AnalyzeLiPPepData}\\code{AnalyzeTrpPepData}\\code{AnalyzeTrpProtData}.
#' @param evalCovariable A character string giving name of columns in which
#' coefficients and p-values of the variable of interest can be found. Default
#' is set to 'Condition_OLS'.
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
summarizeModelResults <- function(resModel, evalCovariable="Condition_OLS",
                                  correctPval="all", pAdjust="fdr",
                                  annotPP=NULL, infoFeature="Peptide",
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
    if(length(intersect(names(vecPval), annotPP[,infoFeature]))==0){
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


#' #' @title Plotting peptides over the protein sequence, visualizing effect sizes and P value
#' #'
#' #' @usage PlottingPepsPerProt(modelSum, infoPval="Padj", infoCoef="Coefficient",
#' #' annotPP, infoFeature="Peptide", infoProtein="Protein",
#' #' sigOnly=TRUE, pvalCut=0.1, savePlots=FALSE, filePath=NULL)
#' #'
#' #' @param modelSum a data.frame object, including coefficients and P values to be plotted
#' #' @param infoPval a character string providing column name of P values to use for plotting.
#' #' Output from \code{summarizeModelResults} provides a data.frame in required format.
#' #' Default is set to 'Padj'.
#' #' @param infoCoef a character string providing column name of coefficients to use for plotting.
#' #' Output from \code{summarizeModelResults} provides a data.frame in required format.
#' #' Default is set to 'Coefficient'.
#' #' @param annotPP a data.frame with peptide and protein annotation.
#' #' Rows are features and columns must include \code{infoFeature} and \code{infoProtein}.
#' #' @param infoFeature a character string providing column name of \code{annotPP} where
#' #' features which were fitted in the models are noted.
#' #' Default is set to 'Peptide'.
#' #' @param infoProtein a character string providing column name of \code{annotPP} where
#' #' protein names/groups are noted. \code{InfoFeature} with the same \code{infoProtein}
#' #' will be plotted together.
#' #' Default is set to 'Protein'.
#' #' @param sigOnly a boolean value, if set to TRUE only proteins with at least
#' #' one significant feature will be plotted.
#' #' Default is set to 'TRUE".
#' #' @param pvalCut a numeric variable, defining significance cut-off.
#' #' Default is set to '0.1'.
#' #' @param savePlots a boolean value, if set to TRUE plots are written as
#' #' .pdf file to location defined in \code{filePath}.
#' #' Default is set to 'FALSE'.
#' #' @param filePath a character string giving location to save plots to,
#' #' if \code{savePlots} is set to TRUE.
#' #'
#' #' @export
#' PlottingPepsPerProt <- function(modelSum, infoPval="Padj", infoCoef="Coefficient",
#'                                 annotPP, infoFeature="Peptide", infoProtein="Protein",
#'                                 sigOnly=TRUE, pvalCut=0.1, savePlots=FALSE, filePath=NULL){
#'     if(savePlots){
#'         if(is.null(filePath)){
#'             stop("Please provide information were to save plots.\n")
#'         }
#'     }
#'     plotData <- CreateData4PlotList(modelSum, infoPval, infoCoef, annotPP, infoFeature, infoProtein)
#'     if(sigOnly){
#'         plotData <- plotData[unlist(lapply, function(x){
#'             any(x[, infoPval] <- pvalCut)
#'         })]
#'     }
#'     plotList <- CreatePepProtPlot(plotData)
#'     if(savePlots){
#'         SavePlotsAsPDF(plotList, filePath)
#'         return(NULL)
#'     }
#'     else{
#'         return(plotList)
#'     }
#' }
#'
#' #' @title Create a peptide list per protein for plotting
#' #'
#' #' @keywords internal
#'
#' CreateData4PlotList <- function(modelSum, infoPval, infoCoef, annotPP, infoFeature, infoProtein){
#'     if(all(row.names(modelSum %in% annotPP[, infoFeature]))){
#'         stop("Not all features of 'modelSum' present in 'annotPP'. \n")
#'     }
#'     plotData <- annotPP[match(row.names(modelSum), annotPP[infoFeature]),]
#'     plotData$Yaxis <- modelSum[, infoCoef]
#'     plotData$Pval <- vapply(modelSum[, infoPval], function(x){
#'         ifelse(x<=0.1, "significant", "not significant")
#'     }, character(1))
#'     plotData <- split(plotData, plotData[, infoProtein])
#'     return(plotData)
#' }
#'
#' #' @title Creating a list of plots
#' #'
#' #' @keywords internal
#' #'
#' CreatePepProtPlot <- function(plotData, Yaxis="Yaxis", start="start", end="end", color="Pval", nameY="Coefficient", showPv=FALSE){
#'     allProt <- names(plotData)
#'     orderSig <- unlist(lapply(plotData), function(x){
#'         min(x[, "Pval"])
#'     })
#'     plotList <- plotList[order(orderSig, decreasing=FALSE)] # order proteins so that proteins with most significant peptides are plotted first
#'     plotList <- lapply(names(plotData), function(i){
#'         ggplotPepProt(df=plotData[[i]],
#'                       title=i,
#'                       nameY=nameY,
#'                       showPv=showPv)
#'
#'     })
#'     return(plotList)
#' }
#'
#' #' @title Plotting peptides protein-wise
#' #'
#' #' @keywords internal
#' #'
#' ggplotPepProt <- function(df, title, nameY, showPv){
#'     if(is.null(df)){
#'         myPlot <- NULL
#'     }
#'     if(showPv){
#'         message("Bad luck, not implemented yet.")
#'     }
#'     else{
#'         myPlot <- ggplot2::ggplot(df, aes(x=start, xend=end, y=Yaxis, yend=Yaxis, col=Pval)) +
#'             ggplot2::geom_segment(size=4) +
#'             ggplot2::geom_hline(yintercept=0, size=1, col="black") +
#'             ggplot2::scale_color_manual(values=c("#000000", "#c11509")) +
#'             ggplot2::ggtitle(title) +
#'             ggplot2::ylab(nameY) +
#'             ggplot2::theme_bw() +
#'             ggplot2::theme(legend.title=ggplot2::element_text(size=10), legend.text=ggplot2::element_text(size=8),
#'                            axis.title=ggplot2::element_text(size=12), axis.text=ggplot2::element_text(size=9),
#'                            plot.title=ggplot2::element_text(size=12, face="bold"))
#'     }
#'     return(myPlot)
#' }
#'
#' #' @title Writing plots as .pdf files
#' #'
#' #' @keywords internal
#'
#' SavePlotsAsPDF <- function(plotList, filePath){
#'     remain <- length(plotList) %% 3
#'     if(remain == 2){
#'         plotList[[length(plotList)+1]] <- NULL
#'     }
#'     if(remain == 1){
#'         plotList[[length(plotList)+1]] <- NULL
#'     }
#'     grDevices::pdf(paste0(filePath, ".pdf"))
#'     sapply(seq(1, length(plotList)/3, 3), function(i){
#'         ggpubr::ggarrange(plotList[[1]], plotList[[2]], plotList[[3]], nrow=3, common.legend=TRUE, legend="right")
#'     })
#'     grDevices::dev.off()
#' }
#'
#'
