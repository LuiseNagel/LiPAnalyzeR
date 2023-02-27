
#' @title Preprocessing peptide and protein quantities
#'
#' @description Function for preprocessing data before fitting models. Different
#' filtering can be applied based on various features such as digest type,
#' proteotypic status or number of missed cleavages.
#'
#' @usage PreprocessQuantityMatrix(SpectroList=NULL, QuantityMatrix=NULL,
#' annotPP=NULL, annotS=NULL, logT=TRUE, filterMinLogQuant=TRUE,
#' thresholdMinLogQuant=10, filterNA=TRUE, maxNAperCondition=0,
#' infoCondition="Condition", filterTryptic=TRUE, infoTryptic="isTryptic",
#' nameFT="Specific", filterProteotryptic=FALSE,
#' infoProteotypic="isProteotypic", nameProteotypic="True",
#' filterMissedCleaved=FALSE, maxMissedCleave=2,
#' infoMissedCleave="NMissedCleavages")
#'
#' @param SpectroList A list of matrices, containing peptides/proteins
#' quantities. Rows represent features and columns refer to the samples. Can
#' simply be output from \code{ExtractDataFromSpectro}.
#' @param QuantityMatrix single matrix with peptide/protein quantities. Rows
#' represent features and columns refer to the samples.
#' @param annotPP A data.frame with peptide and protein annotation. Rows are
#' features and must match to the row.names of SpectroList/QuantityMatrix.
#' Must include all columns needed for filtering of data.
#' @param annotS A data.frame object containing sample annotation. Rows are
#' samples and must match to columns of \code{SpectroList}. If \code{filterNA}
#' is set to 'TRUE', it needs to contain column about different
#' conditions/groups.
#' @param logT A boolean value defining if data should be log-transformed.
#' Default is 'TRUE'.
#' @param filterMinLogQuant A boolean value defining small quantities should be
#' set to NA. Default is 'TRUE'.
#' @param thresholdMinLogQuant A numeric value, quantities below this threshold
#' will be set to NA if \code{filterMinLogQuant} is 'TRUE'. Default is set to
#' 10. Setting it to 12 would be an optional stricter choice.
#' @param filterNA A boolean value defining peptides should be filtered based on
#' number of NAs. Default is 'TRUE'.
#' @param maxNAperCondition A numeric value, defining maximal number of NAs
#' in an individual feature per condition. Default is '0'.
#' @param infoCondition A character string providing column name of
#' \code{annotS} in which condition is provided. Default is 'Condition'. If NA
#' filtering should not be applied on condition level provide only one value in
#' \code{infoCondition} in \code{annotS}.
#' @param filterTryptic A boolean value defining if peptides should be filtered
#' based digest type. Default is 'TRUE', removing all not fully tryptic peptides
#' from the data matrices.
#' @param infoTryptic A character string giving column of \code{annotPP} in
#' which digest annotation is provided. Default is 'isTryptic'.
#' @param nameFT A character string defining which peptide in \code{infoTryptic}
#' should remain in the data. Default is 'Specific'.
#' @param filterProteotryptic A boolean value, defining if peptides should be
#' filtered based on if they are proteotypic. Default is 'FALSE'.
#' @param infoProteotypic A character string giving column of \code{annotPP} in
#' which annotation if peptides is proteotypic is provided. Default is
#' 'isProteotypic'.
#' @param nameProteotypic A character string defining which peptide in
#' \code{infoProteotypic} is proteotypic and will be filtered out if
#' \code{filterProteotryptic} is set to 'TRUE'. Default is 'True'.
#' @param filterMissedCleaved A boolean value defining if peptides should be
#' filtered based on number of missed cleavages. Default is 'FALSE'.
#' @param maxMissedCleave A numeric value, defining maximal number of missed
#' cleavages allowed per peptide if \code{filterMissedCleaved} is set to 'TRUE'.
#' Default is '2'.
#' @param infoMissedCleave a character string giving column of \code{annotPP} in
#' which number of missed cleavage are provided. Default is 'NMissedCleavages'.
#'
#' @return Returns a list of matrices containing preprocessed and filtered
#' peptide and protein quantities  OR a single matrix containing preprocessed
#' and filtered peptide or protein quantities. Rows represent features and
#' columns refer to the samples.
#'
#' @export
PreprocessQuantityMatrix <- function(SpectroList=NULL, QuantityMatrix=NULL,
                                     annotPP=NULL, annotS=NULL, logT=TRUE,
                                     filterMinLogQuant=TRUE,
                                     thresholdMinLogQuant=10, filterNA=TRUE,
                                     maxNAperCondition=0, infoCondition="Condition",
                                     filterTryptic=TRUE,
                                     infoTryptic="isTryptic", nameFT="Specific",
                                     filterProteotryptic=FALSE,
                                     infoProteotypic="isProteotypic",
                                     nameProteotypic="True",
                                     filterMissedCleaved=FALSE,
                                     maxMissedCleave=2,
                                     infoMissedCleave="NMissedCleavages"){


    # check input data
    if(is.null(SpectroList)){
        if(is.null(QuantityMatrix)){
            stop("Please provide an input with peptide or protein quantities
                 using 'SpectroList' or 'QuantityMatrix'")
        }
        else{
            SpectroList <- as.list(QuantityMatrix)
        }
    }
    else{
        message("Preprocessing SpectroList.")
    }
    if((filterTryptic|filterProteotryptic|filterMissedCleaved)&is.null(annotPP)){
        stop("Please provide annotPP.")
    }
    if(filterNA&is.null(annotS)){
        stop("Please provide annotS.")
    }

    # adjusting peptides and samples to be identical in all data
    peps <- Reduce(intersect, lapply(SpectroList, row.names))
    samples <- Reduce(intersect, lapply(SpectroList, colnames))
    if(!is.null(annotPP)){
        peps <- intersect(peps, row.names(annotPP))
        annotPP <- annotPP[peps,]
    }
    if(!is.null(annotS)){
        samples <- intersect(samples, row.names(annotS))
        annotS <- annotS[samples, ]
    }
    SpectroList <- lapply(SpectroList, function(x)
    {x[peps, samples]})

    # Filter FT, proteotypic and missed cleavages
    SpectroList <- lapply(SpectroList, function(x){
        if(logT){
            x <- log2(x)
        }
        if(filterMinLogQuant){
            x[x<thresholdMinLogQuant] <- NA
        }
        if(filterTryptic){
            x <- x[annotPP[, infoTryptic] == nameFT, ]
        }
        if(filterProteotryptic){
            x[annotPP[, infoProteotypic] == nameProteotypic, ]
        }
        if(filterMissedCleaved){
            x[annotPP[, infoMissedCleave] == maxMissedCleave, ]
        }
        return(x)
    })

    # Filter NAs
    if(filterNA){
        SpectroList <- FilterNAsFromList(SpectroList, annotS, infoCondition,
                                         maxNAperCondition)
    }
    if(is.null(SpectroList)){
        SpectroList <- unlist(SpectroList)
    }
    return(SpectroList)
}

#' @title Preprocessing peptide and protein quantities
#'
#' @description Function for preprocessing data before fitting models.
#' Filtering can be applied based on different features such as digest type,
#' proteotypic or number of missed cleavages.
#'
#' @param SpectroList  A list ist of matrices, containing peptides/proteins
#' quantities. Rows represent features and columns refer to the samples. Can be
#' output from \code{ExtractDataFromSpectro}.
#' @param annotS  A data.frame containing sample and condition annotation.
#' Rows are samples and must match to columns of SpectroList/QuantityMatrix.
#' Needs to contain column about different conditions/groups.
#' @param maxNAperCondition A numeric value, defining maximal number of NAs
#' in an individual feature per condition. Default is '0'.
#' @param infoCondition A character string providing column name of
#' \code{annotS} in which condition is provided. Default is 'Condition'. If NA
#' filtering should not be applied on condition level provide only one value in
#' \code{infoCondition} in \code{annotS}.
#'
#' @return list of NA filtered matrices with petide/protein quantities.
#'
FilterNAsFromList <- function(SpectroList, annotS, infoCondition,
                              maxNAperCondition){
    # add matrices together to get NA if any of them is NA
    matNA <- Reduce('+', SpectroList)
    if(maxNAperCondition == 0){
        matNA <- !is.na(rowSums(matNA))
    }

    else{
        # splitting data based on conditions
        matNA <- split(as.data.frame(t(matNA)), annotS[, infoCondition])
        matNA <- lapply(matNA, function(x){
            apply(t(x), 1, function(y){
                sum(is.na(y))<=maxNAperCondition
            })
        })
        matNA <- Reduce("&", matNA)
    }

    SpectroList <- lapply(SpectroList, function(x){
        x <- x[matNA,]
    })
    return(SpectroList)
}
