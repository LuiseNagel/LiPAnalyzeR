#' @title Creating volcano plot
#'
#' @description Creating volcano plot with model coefficients & corresponding
#' p-values
#'
#' @usage makeVolcanoPlot(sumDf, coefCol="Coefficient", pvalCol="Padj",
#' pvalCutoff=0.05,xlim=NULL, ylim=NULL)
#'
#' @param sumDf A data.frame containing coefficients and p-values of the
#' condition of interest were rows corresponding to the quantity of interest
#' (e.g. peptides). Can be the output of \code{summarizeModelResults}.
#' @param coefCol A character vector or numeric value defining column of
#' \code{sumDf} were coefficients of the quantity of interest (e.g. peptides)
#' are provided.
#' Default is 'Coefficient'.
#' @param pvalCol A character vector or numeric value defining column of
#' \code{sumDf} were (adjusted) p-values  of the quantity of interest
#' (e.g. peptides) are provided.
#' Default is 'Padj'.
#' @param pvalCutoff A numeric value, peptides with p-values below this cut-off
#' are considered significant.
#' Default is '0.05'.
#' @param xlim A numeric vector of the length two defining the limits of the
#' x-axis of the plot. If set to NULL, limits of x-axis are chosen based on
#' coefficients provided in \code{sumDf}.
#' @param ylim A numeric vector of the length two defining the limits of the
#' y-axis of the plot. If set to NULL, the lower limit starts at 1 and the upper
#' limit is defined to the smallest provided p-value in \code{pvalCol}.
#'
#' @return Returns a volcano plot (ggplot) of coefficients and matching p-values
#'
#' @export
#'
makeVolcanoPlot <- function(sumDf, coefCol="Coefficient", pvalCol="Padj",
                            pvalCutoff=0.05, xlim=NULL, ylim=NULL){

    ## creating data.frame for plotting
    plotData <- data.frame(Coef=sumDf[, coefCol],
                           log10Pv=-log10(sumDf[,pvalCol]),
                           sig=sumDf[,pvalCol]< pvalCutoff)
    plotData <- plotData[order(plotData$sig), ]

    if(is.null(xlim)){
        xlim <- c(min(stats::na.omit(plotData[, "Coef"])),
                  max(stats::na.omit(plotData[, "Coef"])))
    }
    if(is.null(ylim)){
        ylim <- c(0, max(stats::na.omit(plotData[, "log10Pv"])))
    }

    ## creating volcano plot
    volcanoPlot <- ggplot2::ggplot(mapping=ggplot2::aes(
        x=plotData[, "Coef"],
        y=plotData[, "log10Pv"],
        color=plotData[, "sig"],
        alpha=plotData[, "sig"])) +
        ggplot2::geom_point(size=2) +
        ggplot2::xlim(xlim) +
        ggplot2::ylim(ylim) +
        ggplot2::xlab("Coefficient") +
        ggplot2::ylab("log10(p-value)") +
        ggplot2::scale_color_manual(name=paste0("Pval < ", pvalCutoff),
                                    values=c("#000000", "#c11509")) +
        ggplot2::scale_alpha_manual(name=paste0("Pval < ", pvalCutoff),
                                    values=c(0.3, 1)) +
        ggplot2::theme_bw(base_size=15)

    return(volcanoPlot)
}


#' @title Creating PCA plot
#'
#' @description Estimating and plotting PCA from LiP data, TrP data or RUV
#' residuals.
#'
#' @usage makePCA(dataMat, annotS, infoCondition="Condition",
#' typeCondition="categorical", runProt=FALSE, annotPP=NULL,
#' nameProtQuant="Protein", Xaxis="PC1", Yaxis="PC2", nPCs=5)
#'
#' @param dataMat A matrix or data.frame providing data for estimating the PCA.
#' Rows are features and columns are samples. Can be a LiP or TrP peptide or
#' protein matrix (e.g. one element of the \code{preprocessQuantityMatrix}
#' output) or a matrix of residuals from the RUV model estimated with
#' \code{runModel}\\code{analyzeLiPPepData}\\code{analyzeTrPPepData}\
#' \code{analyzeProtData}.
#' @param annotS A data.frame containing sample annotation. Rows are samples and
#' must match to columns of the matrices in \code{quantityList}. Must include
#' column matching to \code{formulaContrast}.
#' @param infoCondition A character string providing column name of
#' \code{annotS} were condition of interest is provided. PCA will be colored
#' according to condition.
#' Default is 'Condition'.
#' @param typeCondition A character string providing information if
#' \code{inofCondition} is a 'categorical' or 'continuous' variable.
#' Default is 'categorical'.
#' @param runProt A boolean variable, if set to 'TRUE', function will remove
#' duplicated rows in \code{dataMat} and match row names of \code{dataMat} to
#' \code{nameProtQuant} in \code{annotPP}.
#' Default is 'FALSE'.
#' @param annotPP A data.frame with peptides (/modified peptides/precursors) and
#' protein annotation. Rows are features and row.names of \code{dataMat} must be
#' found here. \code{annotPP} must include a column named \code{nameProtQuant}
#' providing protein (group) names. The output from \code{getPepProtAnnot} can
#' be given here.
#' @param nameProtQuant A character string giving column of \code{annotPP} were
#' protein names are provided.
#' Default is 'Protein'.
#' @param Xaxis Numeric variable giving which PC should be visualized on the
#' X-axis. Cannot be be higher than the number of PCs (\code{nPCs}) being
#' estimated.
#' Default is '1'.
#' @param Yaxis  Numeric variable giving which PC should be visualized on the
#' Y-axis. Cannot be be higher than the number of PCs (\code{nPCs}) being
#' estimated.
#' Default is '2'.
#' @param nPCs Numierc variable defining number of PCs to be estimated in PCA.
#' Default is '5'.
#'
#' @return PCA plot colored by condition
#'
#' @export

makePCA <- function(dataMat, annotS, infoCondition="Condition",
                    typeCondition="categorical", runProt=FALSE, annotPP=NULL,
                    nameProtQuant="Protein", Xaxis="PC1", Yaxis="PC2", nPCs=5){

    ## Testing input
    if(!inherits(dataMat, c("matrix", "data.frame"))){
        stop("'dataMat' has to be a data.frame or a matrix.")
    }
    if(!infoCondition %in% colnames(annotS)){
        stop("'infoCondition' has to be present in colnames of 'annotS'.")
    }

    ## matching peptides to proteins if necessary
    if(runProt){
        if(is.null(annotPP)){
            stop("'runProt' is set to 'TRUE', please add 'annotPP' to match the
row names of the 'dataMat' to proteins and remove duplicates. If you the rows
of 'dataMat' are already proteins, please set 'runProt' to 'FALSE'.")
        }
        else{
            Peps2Prots <- split(row.names(dataMat), annotPP[row.names(dataMat),
                                                            nameProtQuant])
            Peps2Prots <- unlist(lapply(Peps2Prots, \(x) (x[1])))
            dataMat <- dataMat[Peps2Prots,]
            row.names(dataMat) <- names(Peps2Prots)
        }
    }

    ## Adding color
    if(tolower(typeCondition) == "categorical"){
        myColTheme <- ggplot2::scale_color_viridis_d(name=infoCondition,
                                                     begin=0.9, end=0)
    }
    else if(tolower(typeCondition) == "continuous"){
        myColTheme <- ggplot2::scale_color_viridis_c(name=infoCondition,
                                                     begin=0.9, end=0)
    }
    else{
        stop("'typeCondition' has to be set to 'categorical' or 'continuous'.")
    }

    ## Estimating PCA & preparing plotting data
    PCAdata <- FactoMineR::PCA(t(dataMat), ncp=nPCs, graph=FALSE)

    plotData <- PCAdata[["ind"]][["coord"]]
    colnames(plotData) <- paste0("PC", seq(1:ncol(plotData)))
    myCol <- annotS[row.names(plotData), infoCondition]

    # plot PCA
    plot <- ggplot2::ggplot(mapping=ggplot2::aes(x=plotData[, Xaxis],
                                                 y=plotData[, Yaxis],
                                                 color=myCol)) +
        ggplot2::geom_point(size=3) +
        ggplot2::xlab(Xaxis) +
        ggplot2::ylab(Yaxis) +
        myColTheme +
        ggplot2::theme_bw(base_size=15)

    return(plot)
}


#' @title Creating wood plots for multiple proteins
#'
#' @description
#'
#' @usage makeWoodsPlotProteinList(sumDf, annotPP, coefCol="Coefficient",
#' pvalCol="Padj", nameProtQuant="Protein", startPosition="startPosition",
#' endPosition="endPosition", pvalCutoff=0.05, deltaColorIsTryptic=FALSE,
#' isTryptic="isTryptic", nameFT=c("Specific"),
#' nameHT=c("Specific-C", "Specific-N"), sigProt=TRUE, protVector=NULL,
#' showPv=FALSE, export=FALSE, file=NULL)
#'
#' @param sumDf A data.frame containing coefficients and p-values of the
#' condition of interest were rows corresponding to the quantity of interest
#' (e.g. peptides). Can be the output of \code{summarizeModelResults}.
#' @param annotPP A data.frame with peptides (/modified peptides/precursors) and
#' protein annotation. Rows are features and must contain the row names of
#' \code{sumDf}. Can be output of \code{getPepProtAnnot}. Must include all
#' columns needed for plotting the data:
#' \itemize{
#'   \item \code{nameProtQuant} protein IDs
#'   \item \code{startPosition} AA positions were peptides (/modified peptides
#'   /precursors) start in protein
#'   \item \code{endPosition} AA positions were peptides (/modified peptides
#'   /precursors) end in protein
#'   \item \code{isTryptic} optional column providing information if peptide is
#'   full-tryptic or half-tryptic, only needed if
#'   \code{deltaColorIsTryptic = 'TRUE'}
#' }
#' @param coefCol A character vector or numeric value defining column of
#' \code{sumDf} were coefficients of the quantity of interest (e.g. peptides)
#' are provided.
#' Default is 'Coefficient'.
#' @param pvalCol A character vector or numeric value defining column of
#' \code{sumDf} were (adjusted) p-values  of the quantity of interest
#' (e.g. peptides) are provided.
#' Default is 'Padj'.
#' @param nameProtQuant A character string giving column of \code{annotPP} which
#' protein names are provided. If a peptides matches to the same protein several
#' times, the protein name should be provided each time, separated by ','. If a
#' peptide maps to the multiple proteins, these different proteins can be
#' provided by separating them with ';'.
#' Default is 'Protein'.
#' @param startPosition A character string or numeric giving the column name or
#' column number in which start position from each peptide in its protein
#' sequence are provided in \code{annotPP}. If a peptides matches to the sam
#' protein several times, different start positions should be separated by ','.
#' If a peptide maps to the multiple proteins, the start positions should be
#' separated with ';'.
#' Default is 'startPosition'.
#' @param endPosition A character string or numeric giving the column name or
#' column number in which end position from each peptide in its protein
#' sequence are provided in \code{annotPP}. If a peptides matches to the sam
#' protein several times, different start positions should be separated by ','.
#' If a peptide maps to the multiple proteins, the start positions should be
#' separated with ';'.
#' Default is 'endPosition'.
#' @param pvalCutoff A numeric value, peptides with p-values below this cut-off
#' are considered significant.
#' Default is '0.05'.
#' @param deltaColorIsTryptic A boolean variable, if set to 'TRUE', full-tryptic
#' and half-tryptic peptides will be plotted in different colors, this requires
#' \code{isTryptic} column in \code{annotPP}
#' Default is FALSE'.
#' @param isTryptic A character string or numeric giving the column name or
#' column number of \code{annotPP} in which digest type can be found
#' Default is 'isTryptic'.
#' @param nameFT A character vector defining the name of full-tryptic peptides
#' provided in the \code{isTryptic} column in \code{annotPP}.
#' Default is c("Specific").
#' @param nameHT A character vector defining the name of half-tryptic peptides
#' provided in the \code{isTryptic} column in \code{annotPP}.
#' Default is c("Specific-C", "Specific-N").
#' @param sigProt A boolean value, if set to TRUE, only proteins with at least
#' one significant peptides are plotted.
#' Default is 'TRUE'.
#' @param protVector A character vector providing names of proteins to
#' be plotted. If set to 'NULL' \code{sigProt} has to be set to 'TRUE'.
#' Default is 'NULL'.
#' @param showPv A boolean value, if set to 'TRUE' p-values as defined in
#' \code{pvalCol} should be displayed in plots.
#' Default is 'FALSE'.
#' @param export A boolean value defining if plots should directly be exported.
#' If set to 'TRUE', plots will be written into a .pdf file, displaying three
#' plots per sheet. Additionally, a list with the plots will still be returned
#' in R. If set to 'TRUE', \code{file} has to be provided.
#' Default is 'FALSE'.
#' @param file A character string providing export path and file name of pdf
#' with plots, if \code{export} is set to 'TRUE'.
#'
#' @export

makeWoodsPlotProteinList <- function(sumDf, annotPP, coefCol="Coefficient",
                                     pvalCol="Padj", nameProtQuant="Protein",
                                     startPosition="startPosition",
                                     endPosition="endPosition", pvalCutoff=0.05,
                                     deltaColorIsTryptic=FALSE,
                                     isTryptic="isTryptic",
                                     nameFT=c("Specific"),
                                     nameHT=c("Specific-C", "Specific-N"),
                                     sigProt=TRUE, protVector=NULL,
                                     showPv=FALSE, export=FALSE, file=NULL){

    ## checking input
    if(length(intersect(row.names(sumDf), row.names(annotPP))) == 0){
        stop("Row names of sumDf and annotPP do not match.")
    }
    if(sigProt == FALSE &is.null(protVector)){
        stop("No proteins chosen. Either set 'sigProt' to TRUE or provide a
character vector giving proteins you want plotted in form of 'protVector'.")
    }
    if(export&is.null(file)){
        stop("'export' is set to TRUE, but path and file name are not defined in
'file'.")
    }

    ## creating data.frame for plotting
    plotData <- getPlottingFormat(sumDf, annotPP, coefCol, pvalCol,
                                  nameProtQuant, startPosition, endPosition,
                                  pvalCutoff, deltaColorIsTryptic, isTryptic,
                                  nameFT, nameHT)

    if(sigProt){
        message("Plotting proteins with at least one significant peptide.")
        plotData <- plotData[unlist(lapply(plotData, \(x){
            any(x$Pval <= pvalCutoff)
        }))]
    }
    else{
        plotData <- plotData[protVector]
    }

    if(all(unlist(lapply(plotData, is.null)))){
        stop("No proteins remaining to be plotted.")
    }

    # creating list of protein plots
    plotList <- lapply(plotData, \(x){
        plottingProtein(plotData=x, deltaColorIsTryptic=deltaColorIsTryptic,
                        showPv=showPv, export=export)
    })

    ## creating and writing pdf if export=TRUE
    if(export){
        exportWoodsPlots(plotList, length(plotList), file)
    }

    return(plotList)
}

#' @title Creating wood plots over a single proteins
#'
#' @description
#'
#' @usage makeWoodsPlotSingleProtein(sumDf, annotPP, coefCol="Coefficient",
#' pvalCol="Padj", nameProtQuant="Protein", startPosition="startPosition",
#' endPosition="endPosition", pvalCutoff=0.05, deltaColorIsTryptic=FALSE,
#' isTryptic="isTryptic", nameFT=c("Specific"),
#' nameHT=c("Specific-C", "Specific-N"), xlim=NULL, ylim=NULL, protName=NULL,
#' showPv=FALSE)
#'
#' @param sumDf A data.frame containing coefficients and p-values of the
#' condition of interest were rows corresponding to the quantity of interest
#' (e.g. peptides). Can be the output of \code{summarizeModelResults}.
#' @param annotPP A data.frame with peptides (/modified peptides/precursors) and
#' protein annotation. Rows are features and must contain the row names of
#' \code{sumDf}. Can be output of \code{getPepProtAnnot}. Must include all
#' columns needed for plotting the data:
#' \itemize{
#'   \item \code{nameProtQuant} protein IDs
#'   \item \code{startPosition} AA positions were peptides (/modified peptides
#'   /precursors) start in protein
#'   \item \code{endPosition} AA positions were peptides (/modified peptides
#'   /precursors) end in protein
#'   \item \code{isTryptic} optional column providing information if peptide is
#'   full-tryptic or half-tryptic, only needed if
#'   \code{deltaColorIsTryptic = 'TRUE'}
#' }
#' @param coefCol A character vector or numeric value defining column of
#' \code{sumDf} were coefficients of the quantity of interest (e.g. peptides)
#' are provided.
#' Default is 'Coefficient'.
#' @param pvalCol A character vector or numeric value defining column of
#' \code{sumDf} were (adjusted) p-values  of the quantity of interest
#' (e.g. peptides) are provided.
#' Default is 'Padj'.
#' @param nameProtQuant A character string giving column of \code{annotPP} which
#' protein names are provided. If a peptides matches to the same protein several
#' times, the protein name should be provided each time, separated by ','. If a
#' peptide maps to the multiple proteins, these different proteins can be
#' provided by separating them with ';'.
#' Default is 'Protein'.
#' @param startPosition A character string or numeric giving the column name or
#' column number in which start position from each peptide in its protein
#' sequence are provided in \code{annotPP}. If a peptides matches to the sam
#' protein several times, different start positions should be separated by ','.
#' If a peptide maps to the multiple proteins, the start positions should be
#' separated with ';'.
#' Default is 'startPosition'.
#' @param endPosition A character string or numeric giving the column name or
#' column number in which end position from each peptide in its protein
#' sequence are provided in \code{annotPP}. If a peptides matches to the sam
#' protein several times, different start positions should be separated by ','.
#' If a peptide maps to the multiple proteins, the start positions should be
#' separated with ';'.
#' Default is 'endPosition'.
#' @param pvalCutoff A numeric value, peptides with p-values below this cut-off
#' are considered significant.
#' Default is '0.05'.
#' @param deltaColorIsTryptic A boolean variable, if set to 'TRUE', full-tryptic
#' and half-tryptic peptides will be plotted in different colors, this requires
#' \code{isTryptic} column in \code{annotPP}
#' Default is FALSE'.
#' @param isTryptic A character string or numeric giving the column name or
#' column number of \code{annotPP} in which digest type can be found
#' Default is 'isTryptic'.
#' @param nameFT A character vector defining the name of full-tryptic peptides
#' provided in the \code{isTryptic} column in \code{annotPP}.
#' Default is c("Specific").
#' @param nameHT A character vector defining the name of half-tryptic peptides
#' provided in the \code{isTryptic} column in \code{annotPP}.
#' Default is c("Specific-C", "Specific-N").
#' @param xlim A numeric vector of the length two defining the limits of the
#' x-axis of the plot. If set to NULL, limits of x-axis are chosen based on
#' sequence positions of the peptides.
#' Default is NULL.
#' @param ylim A numeric vector of the length two defining the limits of the
#' y-axis of the plot. If set to NULL, limits of x-axis are chosen based on
#' the coefficients provided in \code{sumDf}.
#'  Default is NULL.
#' @param protName A character vector giving name of protein to be plotted. Has
#' to be included in the column of \code{nameProtQuant} of \code{annotPP}.
#' @param showPv A boolean value, if set to 'TRUE' p-values as defined in
#' \code{pvalCol} should be displayed in plots.
#' Default is 'FALSE'.
#'
#' @export

makeWoodsPlotSingleProtein <- function(sumDf, annotPP, coefCol="Coefficient",
                                       pvalCol="Padj", nameProtQuant="Protein",
                                       startPosition="startPosition",
                                       endPosition="endPosition",
                                       pvalCutoff=0.05,
                                       deltaColorIsTryptic=FALSE,
                                       isTryptic="isTryptic",
                                       nameFT=c("Specific"),
                                       nameHT=c("Specific-C", "Specific-N"),
                                       xlim=NULL, ylim=NULL, protName=NULL,
                                       showPv=FALSE){

    ## checking input
    if(length(intersect(row.names(sumDf), row.names(annotPP))) == 0){
        stop("Row names of sumDf and annotPP do not match.")
    }
    if(is.null(protName)){
        stop("Please provide the name of the protein you want plotted in
         'protName'.")
    }
    filterProt <- vapply(annotPP[, nameProtQuant], \(x){
        x <- unlist(strsplit(x, ";"))
        protName %in% x
    }, logical(1))

    if(sum(filterProt) == 0){
        stop("'protName' is not present in 'annotPP'.")
    }

    ## get only peptides from protein of interest from annotPP
    annotPP <- annotPP[vapply(annotPP[, nameProtQuant], \(x){
        x <- unlist(strsplit(x, ";"))
        protName %in% x
    }, logical(1)) ,]

    ## creating data.frame for plotting
    plotData <- getPlottingFormat(sumDf, annotPP, coefCol, pvalCol,
                                  nameProtQuant, startPosition, endPosition,
                                  pvalCutoff, deltaColorIsTryptic, isTryptic,
                                  nameFT, nameHT)
    plotData <- plotData[[protName]]

    ## creating plot of protein
    plotProt <- plottingProtein(plotData, deltaColorIsTryptic, xlim, ylim,
                                showPv, export=FALSE)
    return(plotProt)
}

#' @title Function for creating data.frame to create wood plot(s)
#'
#' @usage getPlottingFormat(sumDf, annotPP, coefCol, pvalCol, nameProtQuant,
#' startPosition, endPosition, pvalCutoff, deltaColorIsTryptic, isTryptic,
#' nameFT, nameHT)
#'
#' @param sumDf A data.frame containing coefficients and p-values of the
#' condition of interest were rows corresponding to the quantity of interest
#' (e.g. peptides).
#' @param annotPP A data.frame with peptides (/modified peptides/precursors) and
#' protein annotation. Rows are features and must contain the row names of
#' \code{sumDf}. Must include all columns needed for plotting the data:
#' \itemize{
#'   \item \code{nameProtQuant} protein IDs
#'   \item \code{startPosition} AA positions were peptides (/modified peptides
#'   /precursors) start in protein
#'   \item \code{endPosition} AA positions were peptides (/modified peptides
#'   /precursors) end in protein
#'   \item \code{isTryptic} optional column providing information if peptide is
#'   full-tryptic or half-tryptic, only needed if
#'   \code{deltaColorIsTryptic = 'TRUE'}
#' }
#' @param coefCol A character vector or numeric value defining column of
#' \code{sumDf} were coefficients of the quantity of interest (e.g. peptides)
#' are provided.
#' @param pvalCol A character vector or numeric value defining column of
#' \code{sumDf} were (adjusted) p-values  of the quantity of interest
#' (e.g. peptides) are provided.
#' @param nameProtQuant A character string giving column of \code{annotPP} which
#' protein names are provided. If a peptides matches to the same protein several
#' times, the protein name should be provided each time, separated by ','. If a
#' peptide maps to the multiple proteins, these different proteins can be
#' provided by separating them with ';'.
#' @param startPosition A character string or numeric giving the column name or
#' column number in which start position from each peptide in its protein
#' sequence are provided in \code{annotPP}. If a peptides matches to the sam
#' protein several times, different start positions should be separated by ','.
#' If a peptide maps to the multiple proteins, the start positions should be
#' separated with ';'.
#' @param endPosition A character string or numeric giving the column name or
#' column number in which end position from each peptide in its protein
#' sequence are provided in \code{annotPP}. If a peptides matches to the sam
#' protein several times, different start positions should be separated by ','.
#' If a peptide maps to the multiple proteins, the start positions should be
#' separated with ';'.
#' @param pvalCutoff A numeric value, peptides with p-values below this cut-off
#' are considered significant.
#' @param deltaColorIsTryptic A boolean variable, if set to 'TRUE', full-tryptic
#' and half-tryptic peptides will be plotted in different colors, this requires
#' \code{isTryptic} column in \code{annotPP}.
#' @param isTryptic A character string or numeric giving the column name or
#' column number of \code{annotPP} in which digest type can be found
#' @param nameFT A character vector defining the name of full-tryptic peptides
#' provided in the \code{isTryptic} column in \code{annotPP}.
#' @param nameHT A character vector defining the name of half-tryptic peptides
#' provided in the \code{isTryptic} column in \code{annotPP}.
#'
#' @return Data.frame based on which woods plots are created
getPlottingFormat <- function(sumDf, annotPP, coefCol, pvalCol, nameProtQuant,
                              startPosition, endPosition, pvalCutoff,
                              deltaColorIsTryptic,
                              isTryptic, nameFT, nameHT){


    if(length(intersect(row.names(sumDf), row.names(annotPP))) == 0){
        stop("Chosen protein(s) have no overlapping peptides present in sumDf
and annotPP.")
    }

    sumDf <- sumDf[intersect(row.names(sumDf), row.names(annotPP)), ]
    annotPP <- annotPP[intersect(row.names(sumDf), row.names(annotPP)), ]

    ## Create data frame for plotting the data
    plotData <- sumDf[, c(coefCol, pvalCol)]
    colnames(plotData) <- c("Coef", "Pval")
    plotData$Protein <- annotPP[row.names(plotData), nameProtQuant]
    plotData$quantID <- row.names(plotData)

    ## add information on IsTryptic
    if(deltaColorIsTryptic){
        plotData$isTryptic <- annotPP[row.names(plotData), isTryptic]
    }

    ## add peptide positions
    plotData$start <- annotPP[row.names(plotData), startPosition]
    plotData$end <- annotPP[row.names(plotData), endPosition]

    ## Split protein names for peptides in multiple proteins
    if(sum(grepl(";", plotData$Protein))>0){
        plotData <- separateProteins(plotData, deltaColorIsTryptic)
    }

    ## Split protein names for peptides present multiple times in one protein
    if(sum(grepl(",", plotData$start))>0){
        plotData <- separatePeptides(plotData, deltaColorIsTryptic)
    }

    ## Rename FT and HT column
    if(deltaColorIsTryptic){
        plotData$isTryptic <- vapply(plotData$isTryptic, \(x){
            if(x %in% nameFT){
                x <- "FT"
            }
            else if(x %in% nameHT){
                x <- "HT"
            }
            else{
                x <- as.character(NA)
            }
        }, character(1))
        plotData$Color <- paste0(plotData$isTryptic, "_",
                                 ifelse(plotData[, "Pval"] <= pvalCutoff,
                                        "Sig", "nonSig"))
    }

    else{
        plotData$Color <- ifelse(plotData[, "Pval"] <= pvalCutoff,
                                 "Sig", "nonSig")
    }

    plotData$start <- as.numeric(plotData$start)
    plotData$end <- as.numeric(plotData$end)
    plotData <- plotData[!apply(plotData, 1, \(x) any(is.na(x))),]
    row.names(plotData) <- seq(1, nrow(plotData))

    ## Create a list where every element is one Protein
    protList <- split(plotData, plotData$Protein)

    return(protList)
}

#' @title Function for duplicating peptides if they occur in multiple proteins
#'
#' @description Duplicatig rows of \code{plotData} based on how often the
#' respective quantID occurs in different proteins. Adjusting rest of
#' \code{plotData} format accordingly.
#'
#' @usage separateProteins(plotData, deltaColorIsTryptic)
#'
#' @param plotData Data.frame for plotting as created within
#' \code{getPlottingFormat}.
#' @param deltaColorIsTryptic A boolean variable, if set to 'TRUE', full-tryptic
#' and half-tryptic peptides will be plotted in different colors, this requires
#' \code{isTryptic} column in \code{annotPP}.
#'
#' @return \code{plotData} data.frame
separateProteins <- function(plotData, deltaColorIsTryptic){

    message("At least one quantity of interest (e.g. peptide) matches to more
then one single protein based on the protein information provided in 'annotPP'.
Peptide(s) is/are being matched to all positions and proteins matched in
'annotPP' and will therefore be plotted multiple times. If not all proteins are
plotted, this peptide may not be present in the proteins chosen for plotting.")

    plotData <- do.call(rbind, lapply(seq(1, nrow(plotData)), \(i){
        x <- plotData[i,]
        if(sum(grepl(";", x$Protein)) == 0){
            pepDf <- as.data.frame(x)
        }
        else{
            dPeps <- length(unlist(strsplit(x$Protein, ";")))
            pepDf <- data.frame(Coef=rep(x$Coef, dPeps),
                                Pval=rep(x$Pval, dPeps),
                                Protein=unlist(strsplit(x$Protein, ";")),
                                quantID=rep(x$quantID, dPeps))

            if(dPeps == length(unlist(strsplit(x$start, ";")))){
                pepDf$start <- unlist(strsplit(x$start, ";"))
                pepDf$end <- unlist(strsplit(x$end, ";"))
            }
            else{
                allStart <- unique(unlist(strsplit(x$start, ";")))
                if(length(allStart)==1){
                    allEnd <- unique(unlist(strsplit(x$end, ";")))
                    pepDf$start <- allStart
                    pepDf$end <- allEnd
                }
                else{
                    pepDf$start <- NA
                    pepDf$end <- NA
                }
            }

            if(deltaColorIsTryptic){
                pepDf <- addTrypticInfo(pepDf, x, ";")
            }
            return(pepDf)
        }
    }))

    return(plotData)
}

#' @title Function for duplicating peptides if they occur multiple times in the
#' same protein
#'
#' @description Duplicatig rows of \code{plotData} based on how many start
#' positions the respective quantID as in the same protein. Adjusting rest of
#' \code{plotData} format accordingly.
#'
#' @usage separatePeptides(plotData, deltaColorIsTryptic)
#'
#' @param plotData Data.frame for plotting as created within
#' \code{getPlottingFormat}.
#' @param deltaColorIsTryptic A boolean variable, if set to 'TRUE', full-tryptic
#' and half-tryptic peptides will be plotted in different colors, this requires
#' \code{isTryptic} column in \code{annotPP}.
#'
#' @return \code{plotData} data.frame
separatePeptides <- function(plotData, deltaColorIsTryptic){
    message("At least one peptide matches to more then one position in the
respective protein provided in 'annotPP'. Peptide(s) is/are being matched to
all positions and proteins matched in 'annotPP' and will therefore be plotted
multiple times.If not all proteins are plotted, this peptide may not be present
in the proteins chosen for plotting.")

    plotData <- do.call(rbind, lapply(seq(1, nrow(plotData)), \(i){
        x <- plotData[i,]
        if(sum(grepl(",", x$start)) == 0){
            pepDf <- as.data.frame(x)
        }
        else{
            dPeps <- length(unlist(strsplit(x$start, ",")))
            pepDf <- data.frame(Coef=rep(x$Coef, dPeps),
                                Pval=rep(x$Pval, dPeps),
                                Protein=rep(x$Protein, dPeps),
                                quantID=rep(x$quantID, dPeps),
                                start=unlist(strsplit(x$start, ",")),
                                end=unlist(strsplit(x$end, ",")))
        }

        if(deltaColorIsTryptic){
            pepDf <- addTrypticInfo(pepDf, x, ",")
        }
        return(pepDf)
    }))
    return(plotData)
}


#' @title Function for duplicating peptides if they occur multiple times in the
#' same protein
#'
#' @description Function splits up \code{isTryptic} information in case
#' quantities mapping to multiple positions
#'
#' @usage addTrypticInfo(pepDf, x, sign)
#'
#' @param pepDf Data.frame for plotting protein in woods plot as created within
#' \code{separateProteins} or \code{separatePeptides}
#' @param x Data.frame for plotting protein in woods plot not yet altered by
#' \code{separateProteins} or \code{separatePeptides}
#' @param sign A character variable to separate the 'isTryptic' column in
#' \code{x}
#'
#' @return \code{plotData} data.frame
addTrypticInfo <- function(pepDf, x, sign){
    pepDf$isTryptic <- if(grepl(sign, x$isTryptic)){
        unlist(strsplit(x$isTryptic, sign))
    }
    else{
        rep(x$isTryptic, nrow(pepDf))
    }
    return(pepDf)
}


#' Function for plotting woods plot for single protein
#'
#' @usage plottingProtein(plotData, deltaColorIsTryptic, xlim=NULL, ylim=NULL,
#' showPv, export)
#'
#' @param plotData A data.frame to create the woods plot for one protein from,
#' created using the function \code{getPlottingFormat}.
#' @param deltaColorIsTryptic A boolean variable, if set to 'TRUE', full-tryptic
#' and half-tryptic peptides will be plotted in different colors, this requires
#' \code{isTryptic} column in \code{annotPP}.
#' @param xlim A numeric vector of the length two defining the limits of the
#' x-axis of the plot. If set to NULL, limits of x-axis are chosen based on
#' sequence positions of the peptides.
#' @param ylim A numeric vector of the length two defining the limits of the
#' y-axis of the plot. If set to NULL, limits of x-axis are chosen based on
#' the coefficients provided in \code{sumDf}.
#' @param showPv A boolean value, if set to 'TRUE' p-values as defined in
#' \code{pvalCol} should be displayed in plots.
#' @param export A boolean value defining if plots should directly be exported.
#' If set to 'TRUE', plots will be written into a .pdf file, displaying three
#' plots per sheet. Additionally, a list with the plots will still be returned
#' in R.
#'
#' @return Woods plot for peptides over a protein

plottingProtein <- function(plotData, deltaColorIsTryptic, xlim=NULL, ylim=NULL,
                            showPv, export){

    ## set sizes of plot based on if it is exported as pdf or not
    if(export){
        myTheme <- ggplot2::theme_bw(base_size=12)
        sz <- 4
    }
    else{
        myTheme <- ggplot2::theme_bw(base_size=15)
        sz <- 5
    }

    ## setting colors
    if(deltaColorIsTryptic){
        myColor <- ggplot2::scale_color_manual(name="",
                                               breaks=c("FT_nonSig", "FT_Sig",
                                                        "HT_nonSig",
                                                        "HT_Sig"),
                                               values=c("#000000", "#c11509",
                                                        "#a6a6a6", "#f8a19d"))
    }
    else{
        myColor <- ggplot2::scale_color_manual(name="",
                                               breaks=c("nonSig", "Sig"),
                                               values=c("#000000", "#c11509"))
    }

    plotName <- plotData$Protein[1]

    ## adding xlim and ylim if set to NULL
    if(is.null(xlim)){
        xlim <- c(min(stats::na.omit(plotData[ ,"start"])),
                  max(stats::na.omit(plotData[, "end"])))
    }
    if(is.null(ylim)){
        ylim <- c(min(stats::na.omit(plotData[, "Coef"])),
                  max(stats::na.omit(plotData[, "Coef"])))
    }

    plotProt <- ggplot2::ggplot(mapping=ggplot2::aes(
        x=plotData[, "start"],
        xend=plotData[, "end"],
        y=plotData[, "Coef"],
        yend=plotData[, "Coef"],
        col=plotData[, "Color"])) +
        ggplot2::geom_segment(size=sz) +
        myColor +
        ggplot2::xlim(xlim) +
        ggplot2::ylim(ylim) +
        ggplot2::ggtitle(plotName) +
        ggplot2::ylab("Coefficient") +
        ggplot2::xlab("Protein positions") +
        myTheme

    if(ylim[1]<=0 & ylim[2] >=0){
        plotProt <- plotProt +
            ggplot2::geom_hline(yintercept=0, linewidth=1, col="black")
    }

    if(showPv){
        plotProt <- plotProt +
            ggplot2::geom_label(mapping=ggplot2::aes(x=plotData[, "start"],
                                                     y=plotData[, "Coef"],
                                                     label=round(plotData
                                                                 [, "Pval"],
                                                                 3)),
                                show.legend=F)
    }

    return(plotProt)
}

#' @title Function for exporting woods plots and writing them into a .pdf file
#'
#' @usage exportWoodsPlots(plotList, lList, file)
#'
#' @param plotList A list of woods plots created with the function
#' \code{plottingProtein} were the variable \code{export} was set to 'TRUE.
#' @param lList A numeric variable providing the length of \code{plotList}.
#' @param file A character string providing export path and file name of pdf
#' with plots.
#'
#' @return 'NULL'
exportWoodsPlots <- function(plotList, lList, file){
    message("Writing plots into pdf file.")

    if(!grepl("[.]pdf$", file)){
        file <- paste0(file, ".pdf")
    }

    grDevices::pdf(file)
    lapply(as.list(seq(1, lList, 3)), \(i){
        p1 <- plotList[[i]]
        if(lList < i+1){
            p2 <- patchwork::plot_spacer()
        }
        else{
            p2 <- plotList[[i+1]]
        }
        if(lList < i+2){
            p3 <- patchwork::plot_spacer()
        }
        else{
            p3 <- plotList[[i+2]]
        }

        print(patchwork::wrap_plots(p1, p2, p3, nrow = 3))
    })
    grDevices::dev.off()

    return(NULL)
}
