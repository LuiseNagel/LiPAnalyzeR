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
#' annotPP=NULL, annotS=NULL, mode="default", logT=TRUE, filterMinLogQuant=TRUE,
#' thresholdMinLogQuant=10, filterNA="all", maxNAs=0, infoCondition="Condition",
#' nameIDQuant="quantID", nameProtQuant="Protein", filterTryptic=TRUE,
#' infoTryptic="isTryptic", nameFT="Specific", filterProteotryptic=FALSE,
#' infoProteotypic="isProteotypic", nameProteotypic="True",
#' filterMissedCleaved=FALSE, maxMissedCleave=2,
#' infoMissedCleave="NMissedCleavages")
#'
#' @param quantityList A list of matrices, containing quantities of interest
#' (e.g. peptide, modified peptide, precursor) and protein abundances. Rows
#' represent features and columns samples. Sample names between the different
#' matrices contained in the list should match. Not matching samples will be
#' removed.
#' Output from \code{extractMSData} or \code{extractSpectroData} is in
#' the correct format (except of potentially not matching sample names, then
#' these need to be adjusted manually.)
#' \code{quantityList} may include the following matrices
#' \itemize{
#'   \item 'LiPPep': LiP peptide quantities (or alternatively the modified
#'   peptide/ precursor/other quantities, dependent on \code{quantName} and
#'   \code{quantValue})
#'   \item 'TrPPep': TrP peptide quantities (or alternatively the modified
#'   peptide/precursor/other quantities, dependent on \code{quantName} and
#'   \code{quantValue})
#'   \item 'TrPProt': TrP protein quantities
#'   \item 'LiPProt': LiP protein quantities
#'   }
#' @param quantityMatrix A single matrix with peptide/protein quantities. Rows
#' represent features and columns refer to the samples. Single matrix filtering
#' can only be run in \code{mode='default'}. Do not provide \code{quantityList}
#' if you want to filter \code{quantityMatrix}
#' @param annotPP A data.frame with peptides (/modified peptides/precursors) and
#' protein annotation. Rows are features and must match to the row.names of
#' \code{quantityList} or \code{quantityMatrix}. Must include all columns
#' required for filtering data according to settings provided to this function.
#' The output from \code{getPepProtAnnot} can be given here.
#' @param annotS A data.frame object containing sample annotation. Rows are
#' samples and must match to columns of \code{quantityList}.
#' If \code{filterNA='byCondition'} it needs to contain column providing the
#' conditions/groups.
#' @param mode A character variable defining mode in which function is run. Can
#' be set to c('default', 'HTonly', 'FTHTjoin' or 'LiPonly').
#' \itemize{
#'   \item 'default': Filtering all matrices to contain the same quantity IDs
#'   in rows, e.g. peptide names. If no HT peptides are quantified in the TrP
#'   data, all HT peptides will be removed.
#'   \item 'HTonly': Full-tryptic (FT) LiP peptides will be removed from the
#'   data. All half-tryptic (HT) LiP peptides will be matched to the best
#'   correlating FT TrP peptide from the same protein. All other FT TrP peptides
#'   will be removed from the data. The row names of the TrP data will
#'   subsequently be adjusted to the IDs of the half-tryptic LiP quantites, such
#'   as half-tryptic AA sequences.
#'   \item 'FTHTjoin': The TrPPep matrix will be removed from the data,
#'   matching all FT and HT LiP peptides to the matching TrPProt quantities.
#'   \code{filterTryptic} will automatically be set to 'FALSE'. \code{annotPP}
#'   need to be provided to match TrP protein quantities to half-tryptic
#'   peptides in LiP data.
#'   \item 'LiPonly': Filtering all matrices to contain the same quantity IDs
#'   in rows, e.g. peptide names. If you which to retain half-tryptic peptides
#'   in the data, please set  \code{filterTryptic} to 'FALSE'.
#'   }
#' @param logT A boolean value defining if data should be log-transformed.
#' Default is 'TRUE', meaning the data will be log2 transforming.
#' @param filterMinLogQuant A boolean value defining if small quantities should
#' be set to NA.
#' Default is 'TRUE'.
#' @param thresholdMinLogQuant A numeric value, quantities below this threshold
#' will be set to NA if \code{filterMinLogQuant} is 'TRUE'.
#' Default is set to 10. Setting it to 12 would be an optional stricter choice.
#'
#' @param filterNA A character variable defining type of NA filtering applied.
#' Can be set to c('all', 'byCondition' or 'none').
#' \itemize{
#'   \item 'all': Filtering of NAs will be applied on all samples, independent
#'   of their condition. \code{maxNA} defines the number of NAs allowed per
#'   quantity over all samples.
#'   \item 'byCondition': Filtering of NAs will be applied on samples based on
#'   their condition.  \code{maxNA} defines the number of NAs allowed per
#'   quantity over all samples of each individual condition. It is advised to
#'   only use 'byCondition', if the condition you are providing is a
#'   categorical variable and not continuous.
#'   \item 'none': No NA filtering will be applied. This my cause issues if
#'   running RUV and/or contrast models subsequently.
#'   }
#' Default is 'all'.
#' @param maxNAs A numeric value, defining maximal number of NAs allowed.
#' Default is '0'.
#' @param infoCondition A character string or numeric giving the column name
#' or column number of \code{annotS} in which condition is provided. Required,
#' if \code{filterNA = 'byCondition'}.
#' Default is 'Condition'.
#' @param nameIDQuant A character string giving column of \code{annotPP} which
#' matches to the row.names of \code{quantityList} or \code{quantityMatrix}.
#' Defailt is 'quantID'.
#' @param nameProtQuant A character string giving column of \code{annotPP} which
#' protein names are provided. If a peptides matches to the same protein several
#' times, the protein name should be provided each time, separated by ','. If a
#' peptide maps to the multiple proteins, these different proteins can be
#' provided by separating them with ';'.
#' Default is 'Protein'.
#' @param filterTryptic A boolean value defining if peptides should be filtered
#' based digest type.
#' Default is 'TRUE', keeping only full-tryptic peptides in data matrices.
#' @param infoTryptic A character string or numeric giving the column name or
#' column number of \code{annotPP} in which digest type is provided.
#' Default is 'isTryptic'.
#' @param nameFT A character string defining how full-tryptic peptides are
#' annotated in the \code{infoTryptic} colunm in \code{annotPP}. Only peptides
#' annotated accordingly will remain in the data if \code{filterTryptic} is set
#' to 'TRUE'.
#' Default is 'Specific'.
#' @param filterProteotryptic A boolean value, defining if peptides should be
#' filtered based on if they are proteotypic.
#' Default is 'FALSE'.
#' @param infoProteotypic A character string or numeric giving the column name
#' or column number of \code{annotPP} in which annotation which peptides are
#' proteotypic is provided.
#' Default is 'isProteotypic'.
#' @param nameProteotypic A character string defining proteotypic peptides are
#' annotated in the \code{infoProteotypic} colunm in \code{annotPP}. Onl
#' peptides annotated accordingly will remain in the data if
#' \code{filterProteotryptic} is set to 'TRUE'.
#' Default is 'True'.
#' @param filterMissedCleaved A boolean value defining if peptides should be
#' filtered based on the number of missed cleavages.
#' Default is 'FALSE'.
#' @param maxMissedCleave A numeric value, defining maximal number of missed
#' cleavages allowed per peptide if \code{filterMissedCleaved} is set to 'TRUE'.
#' Default is '2'.
#' @param infoMissedCleave A character string or numeric giving the column name
#' or column number of \code{annotPP} in which number of missed cleavage are
#' provided.
#' Default is 'NMissedCleavages'.
#'
#' @return Returns a list of matrices containing preprocessed and filtered
#' quantities of interest(e.g. peptide, modified peptide, precursor) and protein
#' quantities (OR a single matrix containing preprocessed quantities of interest
#' if a 'quantityMatrix' was provided. Feature and samples names (i.e., row and
#' column names) of the different matrices in the list match.
#'
#' If a "quantityMatrix' is provided the function returns a single matrix
#' containing preprocessed quantities of interest.
#'
#' Rows represent features and columns samples.
#'
#' @export

preprocessQuantityMatrix <- function(quantityList=NULL, quantityMatrix=NULL,
                                     annotPP=NULL, annotS=NULL, mode="default",
                                     logT=TRUE, filterMinLogQuant=TRUE,
                                     thresholdMinLogQuant=10, filterNA="all",
                                     maxNAs=0,
                                     infoCondition="Condition",
                                     nameIDQuant="quantID",
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
    ## check mode setting
    if(!tolower(mode) %in% c("default", "htonly", "fthtjoin", "liponly")){
        stop("Please set 'mode' to one of the allowed options: c('default',
 'HTonly', 'FTHTjoin', 'LiPonly'.")
    }

    ## check quantityList and quantityMatrix input
    if(is.null(quantityList)){
        if(is.null(quantityMatrix)){
            stop("Please provide an input with peptide or protein quantities
using 'quantityList' or 'quantityMatrix'")
        }
        if(tolower(mode) != 'default'){
            stop("'quantityMatrix' can only be filtered if mode is set to
default. Either set mode='default' or provide quantityList as input.")
        }
        else{
            if(!inherits(quantityMatrix, c("matrix", "data.frame"))){
                stop("'quantityMatrix' has to be a data.frame or a matrix.")

            }
            else{
                quantityList <- as.list(quantityMatrix)
            }
        }
    }

    else{
        if(!all(unlist(lapply(quantityList, \(x)
                              inherits(x, c("matrix","data.frame")))))){
            stop("Elements of 'quantityList' have to be data.frames or
matrices.")
        }
        if(!is.null(quantityMatrix)){
            message("'quantityList' and 'quantitMatrix' are provided. Only using
'quantityList'.")
        }
    }


    ## check if annotPP is required and provided
    if((tolower(mode) %in% c("htonly", "fthtjoin")|filterTryptic!="none"|
       filterProteotryptic|filterMissedCleaved)&is.null(annotPP)){
        stop("Please provide annotPP.")
    }

    ## check input to filterTryptic
    if(!(is.logical(filterTryptic)&length(filterTryptic)==1)){
        stop("'filterTryptic' is not set correctly. Please set it to either
'TRUE' or 'FALSE'.")
    }

    ## check input filterNA
    if(!tolower(filterNA) %in% c("all", "bycondition", "none")){
        stop("Please set 'filterNA' to one of the allowed options: c('all',
'byCondition', 'none'.")
    }

    if(tolower(filterNA) == 'bycondition'&is.null(annotS)){
        stop("Please provide annotS for filtering NAs by conidition.")
    }

    message("Preprocessing quantityList.")

    ## removing TrPPep and mapping TrPProt to HT peptides in case the function
    ## is run in 'runHT' mode
    if(tolower(mode) == "fthtjoin"){
        if(filterTryptic){
            message("'mode' is set to 'FTHTjoin', setting 'filterTryptic' to
'FALSE'.")
            filterTryptic <- FALSE
        }

        if(paste(names(quantityList), collapse="") ==
           c("LiPPepTrPPepTrPProt")){
            quantityList <- filterForFTHTjoin(quantityList, annotPP,
                                              nameIDQuant, nameProtQuant)
        }

        else if(paste(names(quantityList), collapse="") != c("LiPPepLiPProt")){
            stop("Function is run in 'FTHTjoin' mode and names of the
quantityList do not meet the expectations. Please provide either a list with
matrices defined as 'LiPPep', 'TrPPep' and 'TrPProt' or a list with two
matrices defined as 'LiPPep' and 'LiPProt'.")
        }
    }

    ## Optional: perform log2 transformation and set low pep quantities to NA
    quantityList <- lapply(quantityList, function(x){
        if(logT){
            x <- log2(x)
        }
        if(filterMinLogQuant){
            x[x<thresholdMinLogQuant] <- NA
        }
        return(x)
    })

    ## if mode is 'HTonly': correlate TrPPep to HT LiPPeps from same protein
    if(tolower(mode) == "htonly"){
        filterTryptic <- FALSE

        if(paste(names(quantityList), collapse="") !=
           c("LiPPepTrPPepTrPProt")){
            stop("Function is run in 'HTonly' mode and the names list elements
in 'quantityList' do not meet the expectations. Please provide a list with three
matrices named 'LiPPep', 'TrPPep' and 'TrPProt'.")
        }
        quantityList <- filterForHTonly(quantityList, annotPP, annotS, filterNA,
                                        maxNAs, infoCondition,
                                        nameIDQuant, nameProtQuant,
                                        infoTryptic, nameFT)
    }

    ## adjusting peptides and samples to be identical in all data
    peps <- Reduce(intersect, lapply(quantityList, row.names))
    nPeps <- length(peps)
    samples <- Reduce(intersect, lapply(quantityList, colnames))

    if(length(peps) == 0){
        message("No overlapping identifiers (row.names) between the provided
data matrices.")
    }

    if(!is.null(annotPP)){
        peps <- intersect(peps, row.names(annotPP))
        annotPP <- annotPP[peps,]
        if(length(peps) == 0){
            message("No overlapping identifiers (row.names) between the provided
data matrices and the 'annotPP' file.")
        }

    }
    if(!is.null(annotS)){
        samples <- intersect(samples, row.names(annotS))
        annotS <- annotS[samples, ]
    }

    quantityList <- lapply(quantityList, function(x){
        x[peps, samples]})

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
    if(tolower(filterNA)!="none"){
        if(length(quantityList)==1){
            quantityList <- quantityList[[1]]
        }
        quantityList <- filterNAsFromList(quantityList, filterNA, annotS,
                                          infoCondition, maxNAs)
    }


    message(paste0(nrow(quantityList[[1]]),
" different quantities (e.g. peptides) from originally ", nPeps,
" remaining in filtered data matrix.
", ncol(quantityList[[1]]), " samples are included in the data."))

    if(is.null(quantityList)){
        quantityList <- unlist(quantityList)
    }


    return(quantityList)
}

#' @title Filtering data for running FHTHjoin mode
#'
#' @description Function for preprocessing data before fitting models.
#' Filtering can be applied based on different features such as digest type,
#' proteotypic or number of missed cleavages.
#'
#' @param quantityList A list of matrices, containing quantities of interest
#' (e.g. peptide, modified peptide, precursor) and protein abundances. Rows
#' represent features and columns samples. Sample names between the different
#' matrices contained in the list should match. Not matching samples will be
#' removed.
#' Output from \code{extractMSData} or \code{extractSpectroData} is in
#' the correct format (except of potentially not matching sample names, then
#' these need to be adjusted manually.)
#' \code{quantityList} myst include the following matrices
#' \itemize{
#'   \item 'LiPPep': LiP peptide quantities (or alternatively the modified
#'   peptide/ precursor/other quantities, dependent on \code{quantName} and
#'   \code{quantValue})
#'   \item 'TrPPep': TrP peptide quantities (or alternatively the modified
#'   peptide/precursor/other quantities, dependent on \code{quantName} and
#'   \code{quantValue})
#'   \item 'TrPProt': TrP protein quantities
#'   }
#' @param annotPP A data.frame with peptides (/modified peptides/precursors) and
#' protein annotation. Rows are features and should be the same identifiers as
#' in the row names of the \code{quantityList}. Data.frame must contain columns
#' providing the \code{nameIDQuant} and matching protein names in
#' \code{nameProtQuant}.
#' @param nameIDQuant A character string giving column of \code{annotPP} which
#' matches to the row.names of \code{quantityList} or \code{quantityMatrix}.
#' @param nameProtQuant A character string giving column of \code{annotPP} which
#' protein names are provided. If a peptides matches to the same protein several
#' times, the protein name should be provided each time, separated by ','. If a
#' peptide maps to the multiple proteins, these different proteins can be
#' provided by separating them with ';'.
#'
#' @return Returns a list of matrices containing full-tryptic and half-tryptic
#' LiP quantities of interest (e.g. peptide, modified peptide, precursor) and
#' corresponding TrP protein quantities. The row names in all matrices
#' correspond to the LiP quantity ID.

filterForFTHTjoin <- function(quantityList, annotPP, nameIDQuant,
                              nameProtQuant){
    message("Running function in 'FTHTjoin' mode, removing TrPPep data and
matching TrPProt data to all FT and HT LiP peptides.")
    LiPPP <- annotPP[match(row.names(quantityList$LiPPep),
                           annotPP[, nameIDQuant]),]
    TrPPP <- annotPP[match(row.names(quantityList$TrPProt),
                           annotPP[, nameIDQuant]), nameProtQuant]

    ## only keep peptides from proteins with TrPProt quantity available
    TrPProt <-  split(quantityList$TrPProt, TrPPP)
    TrPPP <- names(TrPProt)
    TrPProt <- rbindlist(lapply(TrPProt, \(x){x[1,]}))
    row.names(TrPProt) <- TrPPP
    LiPPep <- quantityList$LiPPep[LiPPP[, nameProtQuant] %in% TrPPP, ]
    LiPPP <- LiPPP[match(row.names(LiPPep),
                         LiPPP[, nameIDQuant]), nameProtQuant]

    ## matching TrPProt to LiPPep quantities
    TrPProt <- TrPProt[LiPPP,]
    row.names(TrPProt) <- row.names(LiPPep)
    quantityList <- list(LiPPep=LiPPep,
                         TrPProt=TrPProt)

    return(quantityList)
}

#' @title Filtering data for running 'HTonly' mode
#'
#' @description Function for preprocessing data for 'HTonly' mode prior to
#' running the RUV models. Half-tryptic (HT) peptides in the LiP data will be
#' matched to the best correlating full-tryptic (FT) TrP peptide originating
#' from the same protein(s). If the peptide matches to multiple proteins, TrP
#' peptides from all of these proteins can be matched. All full-tryptic LiP
#' peptides are subsequently removed and TrP peptides and protein qunatity
#' names will be altered to match the names of the corresponding half-tryptic
#' peptides.
#'
#' @param quantityList A list of matrices, containing quantities of interest
#' (e.g. peptide, modified peptide, precursor) and protein abundances. Rows
#' represent features and columns samples. Sample names between the different
#' matrices contained in the list should match. Not matching samples will be
#' removed.
#' Output from \code{extractMSData} or \code{extractSpectroData} is in
#' the correct format (except of potentially not matching sample names, then
#' these need to be adjusted manually.)
#' \code{quantityList} must include the following matrices
#' \itemize{
#'   \item 'LiPPep': LiP peptide quantities (or alternatively the modified
#'   peptide/ precursor/other quantities, dependent on \code{quantName} and
#'   \code{quantValue})
#'   \item 'TrPPep': TrP peptide quantities (or alternatively the modified
#'   peptide/precursor/other quantities, dependent on \code{quantName} and
#'   \code{quantValue})
#'   \item 'TrPProt': TrP protein quantities
#'   }
#' @param annotPP A data.frame with peptides (/modified peptides/precursors) and
#' protein annotation. Rows are features and must match to the row.names of
#' \code{quantityList}. Must include all columns required for filtering
#' data according to settings provided to this function. The output from
#' \code{getPepProtAnnot} can be given here.
#' @param annotS A data.frame object containing sample annotation. Rows are
#' samples and must match to columns of \code{quantityList}.
#' If \code{filterNA='byCondition'} it needs to contain column providing the
#' conditions/groups.
#' @param filterNA A character variable defining type of NA filtering applied.
#' Can be set to c('all', 'byCondition' or 'none').
#' \itemize{
#'   \item 'all': Filtering of NAs will be applied on all samples, independent
#'   of their condition. \code{maxNA} defines the number of NAs allowed per
#'   quantity over all samples.
#'   \item 'byCondition': Filtering of NAs will be applied on samples based on
#'   their condition.  \code{maxNA} defines the number of NAs allowed per
#'   quantity over all samples of each individual condition. It is advised to
#'   only use 'byCondition', if the condition you are providing is a
#'   categorical variable and not continuous.
#'   \item 'none': No NA filtering will be applied. This my cause issues if
#'   running RUV and/or contrast models subsequently.
#'   }
#' @param maxNAs A numeric value, defining maximal number of NAs allowed.
#' @param infoCondition A character string or numeric giving the column name
#' or column number of \code{annotS} in which condition is provided. Required,
#' if \code{filterNA = 'byCondition'}.
#' @param nameIDQuant A character string giving column of \code{annotPP} which
#' matches to the row.names of \code{quantityList} or \code{quantityMatrix}.
#' @param nameProtQuant A character string giving column of \code{annotPP} which
#' protein names are provided. If a peptides matches to the same protein several
#' times, the protein name should be provided each time, separated by ','. If a
#' peptide maps to the multiple proteins, these different proteins can be
#' provided by separating them with ';'.
#' @param infoTryptic A character string or numeric giving the column name or
#' column number of \code{annotPP} in which digest type is provided.
#' @param nameFT A character string defining how full-tryptic peptides are
#' annotated in the \code{infoTryptic} colunm in \code{annotPP}. Only peptides
#' annotated accordingly will remain in the data if \code{filterTryptic} is set
#' to 'TRUE'.
#'
#' @return Returns a list of matrices containing half-tryptic LiP quantities of
#' interest (e.g. peptide, modified peptide, precursor), the best matching
#' full-tryptic TrP quantities of interest (e.g. peptide, modified peptide,
#' precursor) and corresponding TrP protein quantities. The row names in all
#' matrices correspond to the half-tryptic LiP quantity ID.

filterForHTonly <- function(quantityList, annotPP, annotS, filterNA, maxNAs,
                            infoCondition, nameIDQuant, nameProtQuant,
                            infoTryptic, nameFT){

    message("Running function in 'HTonly' mode, finding full-tryptic TrPPep
matching half-tryptic LiP peptides by correlation - this step might take a
few moments.")

    ## match sample names
    samples <- Reduce(intersect, lapply(quantityList, colnames))
    quantityList <- lapply(quantityList, \(x){
        x[, samples]
        return(x)
    })

    ## Filtering LiP and TrP peptide matrix
    if(tolower(filterNA)!="none"){
        LiPPep <- filterNAsFromList(quantityList$LiPPep, filterNA, annotS,
                                    infoCondition, maxNAs)
        TrPPep <- filterNAsFromList(quantityList$TrPPep, filterNA, annotS,
                                    infoCondition, maxNAs)
    }

    else{
        LiPPep <- quantityList$LiPPep
        TrPPep <- quantityList$TrPPep
    }

    ## annotation files for LiP and TrP data seperated
    LiPPP <- annotPP[match(row.names(LiPPep),
                           annotPP[, nameIDQuant]),]
    LiPPP <- LiPPP[LiPPP[, infoTryptic] != nameFT, ] ## Keeping only HT
    TrPPP <- annotPP[match(row.names(TrPPep),
                           annotPP[, nameIDQuant]), ]

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
        myLiPPP <- LiPPP[LiPPP[, nameIDQuant]==x,]
        myProts <- unlist(strsplit(myLiPPP[, nameProtQuant], ";" ))
        myTrPOpt <- do.call(rbind, TrPPP[myProts])
        if(is.null(myTrPOpt)){
            return(NA) # return NA if no matching TrP peptide
        }
        else{
            myTrPOpt <- myTrPOpt[!duplicated(myTrPOpt[, nameIDQuant]),]
            myTrPCor <- stats::cor(t(LiPPep[x,]), t(TrPPep[myTrPOpt
                                                           [, nameIDQuant],]),
                                   use="pairwise.complete.obs")
            myHighCor <- colnames(myTrPCor)[which.max(c(myTrPCor))][1]
            return(myHighCor)
        }
        myHighCor <- names(myTrPCor)[unname(myTrPCor==max(myTrPCor))][1]
        return(myHighCor)
    })

    ## removing peptides with no matching TrP peptide
    myTrPPep <- myTrPPep[!is.na(myTrPPep)]

    ## Filtering data
    LiPPep <- LiPPep[names(myTrPPep),]
    TrPPep <- TrPPep[unname(myTrPPep),]
    TrPProt <- quantityList$TrPProt[unname(myTrPPep),]

    row.names(TrPPep) <- row.names(LiPPep)
    row.names(TrPProt) <- row.names(LiPPep)

    quantityList <- list(LiPPep=LiPPep,
                        TrPPep=TrPPep,
                        TrPProt=TrPProt)

    message("Matching half-tryptic LiP peptides to full-tryptic TrP peptides
completed.")

    return(quantityList)
}


#' @title Filtering quantities based on number of NAs
#'
#' @description Function for removing quantities containing more NAs than
#' defined by the user on a joined quantity level. This reflects the number
#' of samples that cannot be modeled in the RUV model later.
#'
#' @param quantityList A list of matrices, containing quantities of interest
#' (e.g. peptide, modified peptide, precursor) and protein abundances. Rows
#' represent features and columns samples. Sample names between the different
#' matrices contained in the list should match. Not matching samples will be
#' removed.
#' Output from \code{extractMSData} or \code{extractSpectroData} is in
#' the correct format (except of potentially not matching sample names, then
#' these need to be adjusted manually.)
#' \code{quantityList} may include the following matrices
#' \itemize{
#'   \item 'LiPPep': LiP peptide quantities (or alternatively the modified
#'   peptide/ precursor/other quantities, dependent on \code{quantName} and
#'   \code{quantValue})
#'   \item 'TrPPep': TrP peptide quantities (or alternatively the modified
#'   peptide/precursor/other quantities, dependent on \code{quantName} and
#'   \code{quantValue})
#'   \item 'TrPProt': TrP protein quantities
#'   \item 'LiPProt': LiP protein quantities
#'   }
#' @param filterNA A character variable defining type of NA filtering applied.
#' Can be set to c('all', 'byCondition' or 'none').
#' \itemize{
#'   \item 'all': Filtering of NAs will be applied on all samples, independent
#'   of their condition. \code{maxNA} defines the number of NAs allowed per
#'   quantity over all samples.
#'   \item 'byCondition': Filtering of NAs will be applied on samples based on
#'   their condition.  \code{maxNA} defines the number of NAs allowed per
#'   quantity over all samples of each individual condition. It is advised to
#'   only use 'byCondition', if the condition you are providing is a
#'   categorical variable and not continuous.
#'   \item 'none': No NA filtering will be applied. This my cause issues if
#'   running RUV and/or contrast models subsequently.
#'   }
#' @param annotS A data.frame object containing sample annotation. Rows are
#' samples and must match to columns of \code{quantityList}.
#' If \code{filterNA='byCondition'} it needs to contain column providing the
#' conditions/groups.
#' @param maxNAs A numeric value, defining maximal number of NAs allowed.
#' @param infoCondition A character string or numeric giving the column name
#' or column number of \code{annotS} in which condition is provided. Required,
#' if \code{filterNA = 'byCondition'}.
#'
#' @return \code{quantityList} filtered based on number of NAs per quantity ID

filterNAsFromList <- function(quantityList, filterNA, annotS, infoCondition,
                              maxNAs){

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
'quantityList' input to 'filterNAsFromList'.")
    }

    if(tolower(filterNA) == "all"){
        matNA <- apply(matNA, 1, \(x)(sum(is.na(x))<=maxNAs))
    }

    else{
        ## splitting data based on conditions
        matNA <- split(as.data.frame(t(matNA)), annotS[, infoCondition])
        matNA <- lapply(matNA, function(x){
            apply(t(x), 1, function(y){
                sum(is.na(y))<=maxNAs
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
