
#' @title Preprocessing quantity matrices
#'
#' @description Function for preprocessing quantity matrices before fitting RUV
#' and/or contrast models. Samples and quantity IDs (e.g. peptide
#' names/sequences) between different matrices in the quantityList are matched.
#' Samples or quantity IDs no present in all matrices will be removed. Data is
#' log-transformed, filtered for NAs and based on settings various properties
#' such as digest type, number of missed cleavages or proteotypic state may also
#' be used for further filtering the data. If run in HTonly mode, half-tryptic
#' LiP peptides are matched to full-tryptic TrP peptides based on correlation.
#'
#' @usage preprocessQuantityMatrix(quantityList=NULL, quantityMatrix=NULL,
#' annotPP=NULL, annotS=NULL, logT=TRUE, filterMinLogQuant=TRUE,
#' thresholdMinLogQuant=10, filterNA=TRUE, maxNAperCondition=0,
#' infoCondition="Condition", runHT = FALSE,  runHTonly = FALSE,
#' namePepQuant="Peptide", nameProtQuant="Protein", filterTryptic=TRUE,
#' infoTryptic="isTryptic", nameFT="Specific", filterProteotryptic=FALSE,
#' infoProteotypic="isProteotypic", nameProteotypic="True",
#' filterMissedCleaved=FALSE, maxMissedCleave=2,
#' infoMissedCleave="NMissedCleavages")
#'
#' @param quantityList A list of matrices, containing peptides/proteins
#' quantities. Rows represent features and columns refer to the samples. Can
#' simply be output from \code{ExtractDataFromSpectro}.
#' @param quantityMatrix A single matrix with peptide/protein quantities. Rows
#' represent features and columns refer to the samples.
#' @param annotPP A data.frame with peptide and protein annotation. Rows are
#' features and must match to the row.names of quantityList/quantityMatrix.
#' Must include all columns needed for filtering of data.
#' @param annotS A data.frame object containing sample annotation. Rows are
#' samples and must match to columns of \code{quantityList}. If \code{filterNA}
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
#' @param infoCondition A character string or numeric giving the column name
#' or column number of \code{annotS} in which condition is provided. Default is
#' 'Condition'. If NA filtering should not be applied on condition level provide
#' only one value in \code{infoCondition} in \code{annotS}.
#' @param runHT A boolean value defining if the function and further analysis
#' steps should be run including the half-tryptic LiP peptides. If set to true,
#' the TrPPep matrix will be removed, resulting in only the LiPPep and TrPProt
#' matrices remaining in the data and \code{filterTryptic} will automatically be
#' set to 'FALSE'. Need to provide \code{annotPP} to match protein quantities to
#' half-tryptic peptides in LiP data. \code{runHT} is set to 'FALSE' per
#' default.
#' @param runHTonly have to add description
#' @param namePepQuant A character string giving column of \code{annotPP} which
#' matches to the row.names of the peptide quantities.
#' Set to 'Peptide' by default.
#' @param nameProtQuant A character string giving column of \code{annotPP} which
#' matches to the row.names of the protein quantities.
#' Set to 'Protein' by default.
#' @param filterTryptic A boolean value defining if peptides should be filtered
#' based digest type. Default is 'TRUE', removing all not fully tryptic peptides
#' from the data matrices.
#' @param infoTryptic A character string or numeric giving the column name or
#' column number of \code{annotPP} in which digest type can be found
#' Default is 'isTryptic'.
#' @param nameFT A character string defining which annotation in
#' \code{infoTryptic} should remain in the data. Default is 'Specific'.
#' @param filterProteotryptic A boolean value, defining if peptides should be
#' filtered based on if they are proteotypic. Default is 'FALSE'.
#' @param infoProteotypic A character string or numeric giving the column name
#' or column number of \code{annotPP} in which annotation if peptides is
#' proteotypic is provided. Default is 'isProteotypic'.
#' @param nameProteotypic A character string defining which peptide in
#' \code{infoProteotypic} is proteotypic and will be filtered out if
#' \code{filterProteotryptic} is set to 'TRUE'. Default is 'True'.
#' @param filterMissedCleaved A boolean value defining if peptides should be
#' filtered based on number of missed cleavages. Default is 'FALSE'.
#' @param maxMissedCleave A numeric value, defining maximal number of missed
#' cleavages allowed per peptide if \code{filterMissedCleaved} is set to 'TRUE'.
#' Default is '2'.
#' @param infoMissedCleave A character string or numeric giving the column name
#' or column number of \code{annotPP} in which number of missed cleavage are
#' provided. Default is 'NMissedCleavages'.
#'
#' @return Returns a list of matrices containing preprocessed and filtered
#' peptide and protein quantities  OR a single matrix containing preprocessed
#' and filtered peptide or protein quantities. Rows represent features and
#' columns refer to the samples.
#'
#' @export
preprocessQuantityMatrix <- function(quantityList=NULL, quantityMatrix=NULL,
                                     annotPP=NULL, annotS=NULL, logT=TRUE,
                                     filterMinLogQuant=TRUE,
                                     thresholdMinLogQuant=10, filterNA=TRUE,
                                     maxNAperCondition=0,
                                     infoCondition="Condition",
                                     runHT=FALSE, runHTonly=FALSE,
                                     namePepQuant="Peptide",
                                     nameProtQuant="Protein",
                                     filterTryptic=TRUE,
                                     infoTryptic="isTryptic", nameFT="Specific",
                                     filterProteotryptic=FALSE,
                                     infoProteotypic="isProteotypic",
                                     nameProteotypic="True",
                                     filterMissedCleaved=FALSE,
                                     maxMissedCleave=2,
                                     infoMissedCleave="NMissedCleavages"){


    ## check input and data
    if(is.null(quantityList)){
        if(is.null(quantityMatrix)){
            stop("Please provide an input with peptide or protein quantities
                 using 'quantityList' or 'quantityMatrix'")
        }
        else{
            quantityList <- as.list(quantityMatrix)
        }
    }
    else{
        message("Preprocessing quantityList")
    }
    if((filterTryptic|filterProteotryptic|filterMissedCleaved|runHT)&
       is.null(annotPP)){
        stop("Please provide annotPP.")
    }
    if(filterNA&is.null(annotS)){
        stop("Please provide annotS.")
    }

    ## removing TrPPep and mapping TrPProt to HT peptides in case the function
    ## is run in 'runHT' mode
    if(runHT){
        filterTryptic <- FALSE
        if(paste(names(quantityList), collapse = "") ==
           c("LiPPepTrPPepTrPProt")){
            quantityList <- filterForRunHT(quantityList, annotPP,
                                           namePepQuant, nameProtQuant)
        }
        else if(paste(names(quantityList), collapse = "") != c("LiPPepLiPProt")){
            stop("Function is run in 'runHT' mode and names of the quantityList
                 do not meet the expectations. Please provide either a list with
                 three matrices defined as 'LiPPep', 'TrPPep' and 'TrPProt' or a
                 list with three matrices defined as 'LiPPep' and 'LiPProt'.")
        }
    }

    ## Optionally perfrom log2 transformation and set very low peptide
    ## quantities to NA
    quantityList <- lapply(quantityList, function(x){
        if(logT){
            x <- log2(x)
        }
        if(filterMinLogQuant){
            x[x<thresholdMinLogQuant] <- NA
        }
        return(x)
    })

    ## getting most correlated TrPPep and mapping TrPProt to HT peptides,
    ## in case the functio is run in 'runHTonly' mode
    if(runHTonly){
        filterTryptic <- FALSE
        if(paste(names(quantityList), collapse = "") !=
           c("LiPPepTrPPepTrPProt")){
            stop("Function is run in 'runHTonly' mode and names of the
                 quantityList do not meet the expectations. Please provide a list
                 with three matrices defined as 'LiPPep', 'TrPPep' and
                 'TrPProt'.")
        }
        quantityList <- filterForRunHTonly(quantityList, annotPP, annotS,
                                           filterNA, maxNAperCondition,
                                           infoCondition, namePepQuant,
                                           nameProtQuant, infoTryptic, nameFT)

    }

    ## adjusting peptides and samples to be identical in all data
    peps <- Reduce(intersect, lapply(quantityList, row.names))
    samples <- Reduce(intersect, lapply(quantityList, colnames))
    if(!is.null(annotPP)){
        peps <- intersect(peps, row.names(annotPP))
        annotPP <- annotPP[peps,]
    }
    if(!is.null(annotS)){
        samples <- intersect(samples, row.names(annotS))
        annotS <- annotS[samples, ]
    }
    quantityList <- lapply(quantityList, function(x)
    {x[peps, samples]})

    ## filter FT, proteotypic and missed cleavages
    quantityList <- lapply(quantityList, function(x){
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

    ## filter NAs
    if(filterNA){
        if(length(quantityList)==1){
            quantityList <- quantityList[[1]]
        }
        quantityList <- filterNAsFromList(quantityList, annotS, infoCondition,
                                          maxNAperCondition)
    }

    if(is.null(quantityList)){
        quantityList <- unlist(quantityList)
    }

    return(quantityList)
}

#' @title Filtering data for running runHT mode of package
#'
#' @description Function for preprocessing data before fitting models.
#' Filtering can be applied based on different features such as digest type,
#' proteotypic or number of missed cleavages.
#'
#' @param quantityList A list of matrices, containing peptides/proteins
#' quantities. Rows represent features and columns refer to the samples. Can
#' simply be output from \code{ExtractDataFromSpectro}.
#' @param annotPP A data.frame with peptide and protein annotation. Rows are
#' features and must match to the row.names of quantityList.
#' Must include all columns needed for filtering of data.
#' @param namePepQuant A character string giving column of \code{annotPP} which
#' matches to the row.names of the peptide quantities.
#' Set to 'Peptide' by default.
#' @param nameProtQuant A character string giving column of \code{annotPP} which
#' matches to the row.names of the protein quantities.
#' Set to 'Protein' by default.
#'
#' @return A list containing LiPPep and TrPProt quantities

filterForRunHT <- function(quantityList, annotPP, namePepQuant, nameProtQuant){
    message("Running function in 'runHT' mode, removing TrPPep data and matching
            TrPProt data to all LiP peptides.")
    LiPPP <- annotPP[match(row.names(quantityList$LiPPep),
                           annotPP[, namePepQuant]),]
    TrPPP <- annotPP[match(row.names(quantityList$TrPProt),
                           annotPP[, namePepQuant]), nameProtQuant]

    ## only keep peptides from proteins with TrPProt quantity available
    TrPProt <-  split(quantityList$TrPProt, TrPPP)
    TrPPP <- names(TrPProt)
    TrPProt <- do.call(rbind, lapply(TrPProt, \(x){x[1,]}))
    row.names(TrPProt) <- TrPPP
    LiPPep <- quantityList$LiPPep[LiPPP[, nameProtQuant] %in% TrPPP, ]
    LiPPP <- LiPPP[match(row.names(LiPPep),
                         LiPPP[, namePepQuant]), nameProtQuant]

    ## matching TrPProt to LiPPep quantities
    TrPProt <- TrPProt[LiPPP,]
    row.names(TrPProt) <- row.names(LiPPep)
    quantityList <- list(LiPPep = LiPPep,
                        TrPProt = TrPProt)

    return(quantityList)
}

#' @title Filtering data for running runHTonly mode of package
#'
#' @description Function for preprocessing data before fitting models.
#' Filtering can be applied based on different features such as digest type,
#' proteotypic or number of missed cleavages.
#'
#' @param quantityList A list of matrices, containing peptides/proteins
#' quantities. Rows represent features and columns refer to the samples. Can
#' simply be output from \code{ExtractDataFromSpectro}.
#' @param annotPP A data.frame with peptide and protein annotation. Rows are
#' features and must match to the row.names of quantityList.Must include all
#' columns needed for filtering of data.
#' @param filterNA A boolean value defining peptides should be filtered based on
#' number of NAs. Default is 'TRUE'.
#' @param maxNAperCondition A numeric value, defining maximal number of NAs
#' in an individual feature per condition. Default is '0'.
#' @param infoCondition A character string or numeric giving the column name
#' or column number of \code{annotS} in which condition is provided. Default is
#' 'Condition'. If NA filtering should not be applied on condition level provide
#' only one value in \code{infoCondition} in \code{annotS}.
#' @param annotS A data.frame object containing sample annotation. Rows are
#' samples and must match to columns of \code{quantityList}. If \code{filterNA}
#' is set to 'TRUE', it needs to contain column about different
#' conditions/groups.
#' @param namePepQuant A character string giving column of \code{annotPP} which
#' matches to the row.names of the peptide quantities.
#' Set to 'Peptide' by default.
#' @param nameProtQuant A character string giving column of \code{annotPP} which
#' matches to the row.names of the protein quantities.
#' Set to 'Protein' by default.
#' @param infoTryptic A character string or numeric giving the column name or
#' column number of \code{annotPP} in which digest type can be found
#' Default is 'isTryptic'.
#' @param nameFT A character string defining which annotation in
#' \code{infoTryptic} should remain in the data. Default is 'Specific'.
#'
#' @return A list containing LiPPep and TrPProt quantities

filterForRunHTonly <- function(quantityList, annotPP, annotS, filterNA,
                               maxNAperCondition, infoCondition, namePepQuant,
                               nameProtQuant, infoTryptic, nameFT){

    message("Running function in 'runHTonly' mode, finding full-tryptic TrPPep
matching half-tryptic LiP peptides - this step might take a few moments.")

    ## match sample names
    samples <- Reduce(intersect, lapply(quantityList, colnames))
    quantityList <- lapply(quantityList, \(x){
        x[, samples]
        return(x)
    })

    ## Filtering LiP and TrP peptide matrix
    LiPPep <- filterNAsFromList(quantityList$LiPPep, annotS, infoCondition,
                                maxNAperCondition)
    TrPPep <- filterNAsFromList(quantityList$TrPPep, annotS, infoCondition,
                                maxNAperCondition)

    ## annotation files for LiP and TrP data seperated
    LiPPP <- annotPP[match(row.names(LiPPep),
                           annotPP[, namePepQuant]),]
    LiPPP <- LiPPP[LiPPP[, infoTryptic] != nameFT, ] ## Keeping only HT
    TrPPP <- annotPP[match(row.names(TrPPep),
                           annotPP[, namePepQuant]), ]

    ## Getting List of TrPPep matrices, every element holds all peptides
    ## present in a specific protein (not protein group)
    TrPPP <- do.call(rbind, apply(TrPPP, 1, \(x){
        if(grepl(";", x[nameProtQuant])){
            allProts <- unlist(strsplit(x[nameProtQuant], ";" ))
            extendPP <- as.data.frame(t(replicate(length(allProts), x)))
            extendPP[, nameProtQuant] <- allProts
        }
        else{
            extendPP <- as.data.frame(t(x))
        }
        return(extendPP)
    }))
    TrPPP <- split(TrPPP, TrPPP[, nameProtQuant])

    ## Finding TrPPep with the best correlation to the HT LiPPep from the
    ## matching protein(s)
    myTrPPep <- sapply(row.names(LiPPP), \(x){
        myLiPPP <- LiPPP[LiPPP[, namePepQuant]==x,]
        myProts <- unlist(strsplit(myLiPPP[, nameProtQuant], ";" ))
        myTrPOpt <- do.call(rbind, TrPPP[myProts])
        if(is.null(myTrPOpt)){
            return(NA) # return NA if no matching TrP peptide
        }
        else{myTrPOpt <- myTrPOpt[!duplicated(myTrPOpt[, namePepQuant]),]
        myTrPCor <- sapply(myTrPOpt[, namePepQuant], \(y){
            stats::cor(as.numeric(LiPPep[x,]), as.numeric(TrPPep[y,]),
                       use = "pairwise.complete.obs")
        })
        myHighCor <- names(myTrPCor)[unname(myTrPCor==max(myTrPCor))][1]
        return(myHighCor)}
    })

    ## removing peptides with no matching TrP peptide
    myTrPPep <- myTrPPep[!is.na(myTrPPep)]

    ## Filtering data
    LiPPep <- LiPPep[names(myTrPPep),]
    TrPPep <- TrPPep[unname(myTrPPep),]
    TrPProt <- quantityList$TrPProt[unname(myTrPPep),]

    row.names(TrPPep) <- row.names(LiPPep)
    row.names(TrPProt) <- row.names(LiPPep)

    quantityList <- list(LiPPep = LiPPep,
                        TrPPep = TrPPep,
                        TrPProt = TrPProt)
    return(quantityList)
}


#' @title Filtering peptide and protein quantities based on NAs
#'
#' @description Function for removing peptides with too many NAs on one of the
#' peptide and protein quantity levels.
#'
#' @param quantityList  A list of matrices, containing peptides/proteins
#' quantities. Rows represent features and columns refer to the samples. Can be
#' output from \code{ExtractDataFromSpectro}.
#' @param annotS  A data.frame containing sample and condition annotation.
#' Rows are samples and must match to columns of quantityList.
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
filterNAsFromList <- function(quantityList, annotS, infoCondition,
                              maxNAperCondition){

    ## if list input matrices together to get NA if any of them is NA
    ## else set matNA to be the input matrix
    if(inherits(quantityList, "list")){
        matNA <- Reduce('+', quantityList)
    }
    else if(inherits(quantityList, "data.frame")|
            inherits(quantityList, "matrix")){
        matNA <- quantityList
    }
    else{
        stop("Wrong input format. Expecting list, data.frame or matrix as
'quantityList' input.")
    }

    if(maxNAperCondition == 0){
        matNA <- !is.na(rowSums(matNA))
    }

    else{
        ## splitting data based on conditions
        matNA <- split(as.data.frame(t(matNA)), annotS[, infoCondition])
        matNA <- lapply(matNA, function(x){
            apply(t(x), 1, function(y){
                sum(is.na(y))<=maxNAperCondition
            })
        })
        matNA <- Reduce("&", matNA)
    }

    ## Filter rows with to many NAs
    if(inherits(quantityList, "list")){
        quantityList <- lapply(quantityList, function(x){
            x <- x[matNA,]
        })
    }
    else{
        quantityList <- quantityList[matNA,]
    }

    return(quantityList)
}
