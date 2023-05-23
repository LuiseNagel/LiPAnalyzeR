globalVariables(names=c("Coef", "Color"))

#' @title Creating wood plots over multiple proteins
#' @usage makeWoodsPlotProteinList(sumDf, annotPP, coefCol="Coefficient",
#' pvalCol="Padj", infoProtein="Protein", infoStart="startPosition",
#' infoEnd="endPosition", pvalCutoff=0.05, deltaColorIsTryptic=FALSE,
#' infoTryptic="isTryptic", nameFT=c("Specific"),
#' nameHT=c("Specific-C", "Specific-N"), sigProt=TRUE, protList=NULL,
#' showPv=FALSE, export=FALSE, file=NULL)
#'
#' @param sumDf A data.frame containing coefficients and p-values of
#' peptides, with rows corresponding to the peptides. Can be the output of
#' \code{summarizeModelResults}
#' @param annotPP A data.frame with peptide and protein annotation. Rows are
#' peptides and must match to the row.names of \code{sumDf}. Can be output of
#' \code{getPepProtAnnot}.Must include all columns needed for plotting the data
#' (Protein name, start and stop position if each peptide, if
#' \code{deltaColorIsTryptic} is set to TRUE, information if each peptide is
#' full or half tryptic).
#' @param coefCol A character vector or numeric value defining column of
#' \code{sumDf} in which coefficients of the peptides are provided. Set to
#' 'Coefficient' per default.
#' @param pvalCol A character vector or numeric value defining column of
#' \code{sumDf} in which p-values of the peptides are provided. Set to
#' 'Padj' per default.
#' @param infoProtein A character vector or numeric value defining column of
#' \code{annotPP} in which the protein name of the peptides are provided. Set to
#' 'Protein' per default. If a peptide is mapped and should be plotted to
#' multiple proteins, these should be seperated within the column by ';'.
#' @param infoStart A character vector or numeric value defining column of
#' \code{annotPP} in which the start position the peptides are provided. Set to
#' 'startPosition' per default. If a peptide is mapped and should be plotted to
#' multiple proteins, these should be seperated within the column by ';'. If a
#' peptide matches to multiple positions within the same protein, seperate by
#' ','.
#' @param infoEnd A character vector or numeric value defining column of
#' \code{annotPP} in which the end position the peptides are provided. Set to
#' 'startPosition' per default. If a peptide is mapped and should be plotted to
#' multiple proteins, these should be seperated within the column by ';'. If a
#' peptide matches to multiple positions within the same protein, seperate by
#' ','.
#' @param pvalCutoff A numeric value, peptides with p-values below this cut-off
#' are considered to be significant. Default is set to '0.05'.
#' @param deltaColorIsTryptic A boolean variable, if set to 'TRUE', full and
#' half tryptic peptides will be plotted in different colors. Default is
#' 'FALSE'.
#' @param infoTryptic A character string or numeric giving the column name or
#' column number of \code{annotPP} in which digest type can be found
#' Default is 'isTryptic'.
#' @param nameFT A character vector defining the name of fully tryptic peptides
#' provided in the \code{isTryptic} column in \code{annotPP}. Default is set to
#' c("Specific").
#' @param nameHT A character vector defining the name of half tryptic peptides
#' provided in the \code{isTryptic} column in \code{annotPP}. Default is set to
#' c("Specific-C", "Specific-N").
#' @param sigProt A boolean value, if set to TRUE, only proteins with
#' significant peptides are plotted. Default is set to 'TRUE'.
#' @param protList A character vector providing names of proteins to be plotted.
#' @param showPv A boolean value defining if p-values should be displayed in
#' plots. Default is set to 'FALSE'.
#' @param export A boolean value defining if plots should directly be exported.
#' If set to 'TRUE', plots will be written into a pdf file, displaying three
#' plots per sheet. Additionally, a list with the plots will still be returned
#' in R. If set to 'TRUE', \code{file} has to be provided, Default is set to
#' 'FALSE'.
#' @param file A character string providing export path and file name of pdf
#' with plots, if \code{export} is set to 'TRUE'.
#'
#' @export

makeWoodsPlotProteinList <- function(sumDf, annotPP, coefCol="Coefficient",
                                     pvalCol="Padj", infoProtein="Protein",
                                     infoStart="startPosition",
                                     infoEnd="endPosition",
                                     pvalCutoff=0.05,
                                     deltaColorIsTryptic=FALSE,
                                     infoTryptic="isTryptic",
                                     nameFT=c("Specific"),
                                     nameHT=c("Specific-C", "Specific-N"),
                                     sigProt=TRUE, protList=NULL, showPv=FALSE,
                                     export=FALSE, file=NULL){

    ## checking input
    if(length(intersect(row.names(sumDf), row.names(annotPP))) == 0){
        stop("Row names of sumDf and annotPP do not match.")
    }
    if(sigProt == FALSE &is.null(protList)){
        stop("No proteins chosen. Either set 'sigProt' to TRUE or provide a
    character vector giving proteins you want plotted in form of 'protList'.")
    }
    if(export&is.null(file)){
        stop("'export' is set to TRUE, but path and file name are not defined in
         'file'.")
    }

    ## creating data.frame for plotting
    plotData <- getPlottingFormat(sumDf, annotPP, coefCol, pvalCol, infoProtein,
                                  infoStart, infoEnd, pvalCutoff,
                                  deltaColorIsTryptic, infoTryptic, nameFT,
                                  nameHT)

    if(sigProt){
        message("Plotting proteins with at least one significant peptide.")
        plotData <- plotData[unlist(lapply(plotData, \(x){
            any(x$Pval <= pvalCutoff)
        }))]
    }
    else{
        plotData <- plotData[protList]
    }

    if(all(unlist(lapply(plotData, is.null)))){
        stop("No proteins remaining to be plotted.")
    }

    # creating list of protein plots
    plotList <- lapply(plotData, \(x){
        plottingProtein(plotData = x, deltaColorIsTryptic = deltaColorIsTryptic,
                        showPv = showPv, export = export)
    })

    ## creating and writing pdf if export = TRUE
    if(export){
        exportWoodsPlots(plotList, length(plotList), file)
    }

    return(plotList)
}

#' @title Creating wood plots over a single protein
#' @usage makeWoodsPlotSingleProtein(sumDf, annotPP, coefCol="Coefficient",
#' pvalCol="Padj", infoProtein="Protein", infoStart="startPosition",
#' infoEnd="endPosition", pvalCutoff=0.05, deltaColorIsTryptic=FALSE,
#' infoTryptic="isTryptic", nameFT=c("Specific"),
#' nameHT=c("Specific-C", "Specific-N"), xlim=NULL, ylim=NULL, protName=NULL,
#' showPv=FALSE)
#'
#' @param sumDf A data.frame containing coefficients and p-values of
#' peptides, with rows corresponding to the peptides. Can be the output of
#' \code{summarizeModelResults}
#' @param annotPP A data.frame with peptide and protein annotation. Rows are
#' peptides and must match to the row.names of \code{sumDf}. Can be output of
#' \code{getPepProtAnnot}.Must include all columns needed for plotting the data
#' (Protein name, start and stop position if each peptide, if
#' \code{deltaColorIsTryptic} is set to TRUE, information if each peptide is
#' full or half tryptic).
#' @param coefCol A character vector or numeric value defining column of
#' \code{sumDf} in which coefficients of the peptides are provided. Set to
#' 'Coefficient' per default.
#' @param pvalCol A character vector or numeric value defining column of
#' \code{sumDf} in which p-values of the peptides are provided. Set to
#' 'Padj' per default.
#' @param infoProtein A character vector or numeric value defining column of
#' \code{annotPP} in which the protein name of the peptides are provided. Set to
#' 'Protein' per default. If a peptide is mapped and should be plotted to
#' multiple proteins, these should be seperated within the column by ';'.
#' @param infoStart A character vector or numeric value defining column of
#' \code{annotPP} in which the start position the peptides are provided. Set to
#' 'startPosition' per default. If a peptide is mapped and should be plotted to
#' multiple proteins, these should be seperated within the column by ';'. If a
#' peptide matches to multiple positions within the same protein, seperate by
#' ','.
#' @param infoEnd A character vector or numeric value defining column of
#' \code{annotPP} in which the end position the peptides are provided. Set to
#' 'startPosition' per default. If a peptide is mapped and should be plotted to
#' multiple proteins, these should be seperated within the column by ';'. If a
#' peptide matches to multiple positions within the same protein, seperate by
#' ','.
#' @param pvalCutoff A numeric value, peptides with p-values below this cut-off
#' are considered to be significant. Default is set to '0.05'.
#' @param deltaColorIsTryptic A boolean variable, if set to 'TRUE', full and
#' half tryptic peptides will be plotted in different colors. Default is
#' 'FALSE'.
#' @param infoTryptic A character string or numeric giving the column name or
#' column number of \code{annotPP} in which digest type can be found
#' Default is 'isTryptic'.
#' @param nameFT A character vector defining the name of fully tryptic peptides
#' provided in the \code{isTryptic} column in \code{annotPP}. Default is set to
#' c("Specific").
#' @param nameHT A character vector defining the name of half tryptic peptides
#' provided in the \code{isTryptic} column in \code{annotPP}. Default is set to
#' c("Specific-C", "Specific-N").
#' @param xlim A numeric vector of the length two defining the limits of the
#' x-axis of the plot. If set to NULL, limits of x-axis are chosen based on
#' sequence positions of the peptides. Default is set to NULL.
#'
#' @param ylim A numeric vector of the length two defining the limits of the
#' y-axis of the plot. If set to NULL, limits of x-axis are chosen based on
#' the coefficients of the peptides. Default is set to NULL.
#' @param protName A character vector giving name of protein to be plotted.
#' @param showPv A boolean value defining if p-values should be displayed in
#' plots. Default is set to 'FALSE'.
#'
#' @export

makeWoodsPlotSingleProtein <- function(sumDf, annotPP, coefCol="Coefficient",
                                       pvalCol="Padj", infoProtein="Protein",
                                       infoStart="startPosition",
                                       infoEnd="endPosition",
                                       pvalCutoff=0.05,
                                       deltaColorIsTryptic=FALSE,
                                       infoTryptic="isTryptic",
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
    filterProt <- vapply(annotPP[, infoProtein], \(x){
        x <- unlist(strsplit(x, ";"))
        protName %in% x
    }, logical(1))

    if(sum(filterProt) == 0){
        stop("'protName' is not present in 'annotPP'.")
    }

    ## get only peptides from protein of interest from annotPP
    annotPP <- annotPP[vapply(annotPP[, infoProtein], \(x){
        x <- unlist(strsplit(x, ";"))
        protName %in% x
    }, logical(1)) ,]

    ## creating data.frame for plotting
    plotData <- getPlottingFormat(sumDf, annotPP, coefCol, pvalCol, infoProtein,
                                  infoStart, infoEnd, pvalCutoff,
                                  deltaColorIsTryptic, infoTryptic, nameFT,
                                  nameHT)
    plotData <- plotData[[protName]]

    ## creating plot of protein
    plotProt <- plottingProtein(plotData, deltaColorIsTryptic, xlim, ylim,
                                showPv, export = FALSE)
    return(plotProt)
}


getPlottingFormat <- function(sumDf, annotPP, coefCol, pvalCol, infoProtein,
                              infoStart, infoEnd, pvalCutoff,
                              deltaColorIsTryptic,
                              infoTryptic, nameFT, nameHT){


    if(length(intersect(row.names(sumDf), row.names(annotPP))) == 0){
        stop("Chosen protein(s) have no overlapping peptides present in sumDf and
         annotPP.")
    }

    sumDf <- sumDf[intersect(row.names(sumDf), row.names(annotPP)), ]
    annotPP <- annotPP[intersect(row.names(sumDf), row.names(annotPP)), ]


    ## Create data frame for plotting the data
    plotData <- sumDf[, c(coefCol, pvalCol)]
    colnames(plotData) <- c("Coef", "Pval")
    plotData$Protein <- annotPP[row.names(plotData), infoProtein]
    plotData$Peptide <- row.names(plotData)

    ## add information on IsTryptic
    if(deltaColorIsTryptic){
        plotData$isTryptic <- annotPP[row.names(plotData), infoTryptic]
    }

    ## add peptide positions
    plotData$start <- annotPP[row.names(plotData), infoStart]
    plotData$end <- annotPP[row.names(plotData), infoEnd]

    ## Split protein names for peptides in multiple proteins
    if(sum(grepl(";", plotData$Protein))>0){
        plotData <- seperateProteins(plotData, deltaColorIsTryptic)
    }

    ## Split protein names for peptides present multiple times in one protein
    if(sum(grepl(",", plotData$start))>0){
        plotData <- seperatePeptides(plotData, deltaColorIsTryptic)
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

## Function for adding peptides as often as they are presented in different
## proteins or in different positions in the data frame
seperateProteins <- function(plotData, deltaColorIsTryptic){
    message("At least one peptide matches to more then one single protein based
    on the protein information provided in 'annotPP'. Peptide(s) is/are being
    matched to all positions and proteins matched in 'annotPP' and will
    therefore be plotted multiple times.")

    plotData <- do.call(rbind, lapply(seq(1, nrow(plotData)), \(i){
        x <- plotData[i,]
        if(sum(grepl(";", x$Protein)) == 0){
            pepDf <- as.data.frame(x)
        }
        else{
            dPeps <- length(unlist(strsplit(x$Protein, ";")))
            pepDf <- data.frame(Coef = rep(x$Coef, dPeps),
                                Pval = rep(x$Pval, dPeps),
                                Protein = unlist(strsplit(x$Protein, ";")),
                                Peptide = rep(x$Peptide, dPeps))

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

## Function for adding peptides as often as they are presented in different
## proteins or in different positions in the data frame
seperatePeptides <- function(plotData, deltaColorIsTryptic){
    message("At least one peptide matches to more then one position in the
  respective protein provided in 'annotPP'. Peptide(s) is/are being matched to
  all positions and proteins matched in 'annotPP' and will therefore be plotted
  multiple times.")

    plotData <- do.call(rbind, lapply(seq(1, nrow(plotData)), \(i){
        x <- plotData[i,]
        if(sum(grepl(",", x$start)) == 0){
            pepDf <- as.data.frame(x)
        }
        else{
            dPeps <- length(unlist(strsplit(x$start, ",")))
            pepDf <- data.frame(Coef = rep(x$Coef, dPeps),
                                Pval = rep(x$Pval, dPeps),
                                Protein = rep(x$Protein, dPeps),
                                Peptide = rep(x$Peptide, dPeps),
                                start = unlist(strsplit(x$start, ",")),
                                end = unlist(strsplit(x$end, ",")))
        }

        if(deltaColorIsTryptic){
            pepDf <- addTrypticInfo(pepDf, x, ",")
        }
        return(pepDf)
    }))
    return(plotData)
}

## Function for splitting up IsTryptic information in case peptides map to
## multiple positions
addTrypticInfo <- function(pepDf, x, sign){
    pepDf$isTryptic <- if(grepl(sign, x$isTryptic)){
        unlist(strsplit(x$isTryptic, sign))
    }
    else{
        rep(x$isTryptic, nrow(pepDf))
    }
    return(pepDf)
}

## Function for plotting wood plot pf single protein
plottingProtein <- function(plotData, deltaColorIsTryptic, xlim=NULL, ylim=NULL,
                            showPv, export){

    ## set sizes of plot based on if it is exported as pdf or not
    if(export){
        myTheme <- ggplot2::theme_bw(base_size = 12)
        sz <- 4
    }
    else{
        myTheme <- ggplot2::theme_bw(base_size = 15)
        sz <- 5
    }

    ## setting colors
    if(deltaColorIsTryptic){
        myColor <- ggplot2::scale_color_manual(name = "",
                                               breaks = c("FT_nonSig", "FT_Sig",
                                                          "HT_nonSig",
                                                          "HT_Sig"),
                                               values = c("#000000", "#c11509",
                                                          "#a6a6a6", "#f8a19d"))
    }
    else{
        myColor <- ggplot2::scale_color_manual(name = "",
                                               breaks = c("nonSig", "Sig"),
                                               values = c("#000000", "#c11509"))
    }

    plotName <- plotData$Protein[1]

    ## adding xlim and ylim if set to NULL
    if(is.null(xlim)){
        xlim <- c(min(stats::na.omit(plotData$start)),
                  max(stats::na.omit(plotData$end)))
    }
    if(is.null(ylim)){
        ylim <- c(min(stats::na.omit(plotData$Coef)),
                  max(stats::na.omit(plotData$Coef)))
    }

    plotProt <- ggplot2::ggplot(plotData, ggplot2::aes(x = start, xend = end,
                                                       y = Coef, yend = Coef,
                                                       col = Color)) +
        ggplot2::geom_segment(size = sz) +
        ggplot2::geom_hline(yintercept = 0, size = 1, col = "black") +
        myColor +
        ggplot2::xlim(xlim) +
        ggplot2::ylim(ylim) +
        ggplot2::ggtitle(plotName) +
        ggplot2::ylab("Coefficient") +
        ggplot2::xlab("Protein positions") +
        myTheme

    if(showPv){
        plotProt <- plotProt +
            ggplot2::geom_label(data = plotData, aes(x = start, y = Coef,
                                                     label = round(plotData[, "Pval"], 3)),
                                show.legend = F)
    }

    return(plotProt)
}

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
        print(p1 + p2 + p3 + patchwork::plot_layout(nrow = 3))
    })
    grDevices::dev.off()

    return(NULL)
}


