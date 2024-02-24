globalVariables(names = "quantName")

#' @title Creating quantity matrices from a MS report.
#'
#' @description Extract the quantity of interest, e.g. peptide intensities, and
#' protein quantities from a MS report/ MS reports.
#'
#' @usage extractMSData(reportLiP, reportTrP=NULL, sampleName, quantName,
#' quantValue, protValue, LiPonly=FALSE)
#'
#' @param reportLiP MS report of LiP data. Has to include columns with
#' the name of the samples (\code{sampleName}), name of quantities which will be
#' the row names of the quantity matrices (e.g. peptide names,
#' \code{quantName}), quantities matching the \code{quantName} (e.g. peptide
#' quantities, \code{quantValue}) and the protein quantities (\code{protValue}).
#' @param reportTrP MS report of the TrP data. Has to include the same
#' columns as the \code{reportLiP}. Does not have to be provided if the LiPonly
#' mode is run. For this set (\code{LiPonly} = TRUE). If LiP and TrP data were
#' processed together, please separate them into two input files.
#' @param sampleName A character string or numeric giving the column name or
#' column number in which sample names are provided in the MS report.
#' @param quantName A character string or numeric giving the column name or
#' column number in which the row names of the quantity matrices created in this
#' function are provided (e.g. peptide, modified peptide, precursor).
#' @param quantValue A character string or numeric giving the column name or
#' column number with quantities to infer structural alterations from (e.g.
#' peptide quantities, modified peptide quantities, precursors). These have to
#' match to the ID provided in \code{quantName}.
#' @param protValue A character string or numeric giving the column name or
#' column number were protein quantities are provided.
#' @param LiPonly A boolean value, set to 'TRUE' if you want to run the LiPonly
#' version o the package and not use trypsin-only data. Default is set to
#' 'FALSE'.
#'
#' @return
#' If run in default mode, a list of three matrices is returned:
#' \itemize{
#'   \item 'LiPPep': LiP peptide quantities (or alternatively the modified
#'   peptide/ precursor/other quantities, dependent on \code{quantName} and
#'   \code{quantValue})
#'   \item 'TrPPep': TrP peptide quantities (or alternatively the modified
#'   peptide/precursor/other quantities, dependent on \code{quantName} and
#'   \code{quantValue})
#'   \item TrP protein quantities
#'   }
#'
#' If run in \code{LiPonly} mode, a list of two matrices is returned:
#' \itemize{
#'   \item 'LiPPep': LiP peptide quantities (or alternatively the modified
#'   peptide/precursor/other quantities, dependent on \code{quantName} and
#'   \code{quantValue})
#'   \item 'LiPProt': LiP protein quantities
#' }
#'
#' Rows represent features and columns samples.
#'
#'
#' @export

extractMSData <- function(reportLiP, reportTrP=NULL, sampleName, quantName,
                          quantValue, protValue, LiPonly=FALSE){

    ## checking if 'reportTrP' is provided in case that LiPonly is FALSE
    if(!LiPonly & is.null(reportTrP)){
        stop("Please provide 'reportTrP' data or set LiPonly to 'TRUE'.")
    }

    ## extracting feature matrix for LiP quantities of interest
    message("Creating LiP matrix with '", quantName, "'.")
    LiPPep <- convert2Matrix(reportLiP, quantValue, quantName, sampleName)

    ## extracting feature matrix of the LiP proteins if LiPonly mode is run
    if(LiPonly){
        message("Creating LiP protein matrix with '", quantName, "'.")
        LiPProt <- convert2Matrix(reportLiP, protValue, quantName, sampleName)
        out <- list(LiPPep=LiPPep, LiPProt=LiPProt)
    }

    ## extracting feature matrix for TrP quantities
    else{
        message("Creating TrP matrix with '", quantName, "'.")
        TrPPep <- convert2Matrix(reportTrP, quantValue, quantName, sampleName)
        message("Creating TrP protein matrix with '", protValue, "'.")
        TrPProt <- convert2Matrix(reportTrP, protValue, quantName, sampleName)
        out <- list(LiPPep=LiPPep, TrPPep=TrPPep, TrPProt=TrPProt)
    }
    return(out)
}

#' @title Creating quantity matrices from a Spectronaut report.
#'
#' @description Extract the quantity of interest, e.g. peptide intensities,
#' and protein quantities from a Spectronaut report/Spectronaut reports.
#' Data must have been exported using the Spectronaut schema
#' 'SpectroSchema_LiPAnalyzerOut', else, please use the function
#' /code{extractMSData}.
#'
#' @usage extractSpectroData(spectroLiP, spectroTrP=NULL, analysisLvl="Peptide",
#' sampleName="R.FileName", LiPonly=FALSE)
#'
#' @param spectroLiP Spectronaut report of LiP data. Spectronaut report has to
#' be exported using the Spectronaut schema SpectroSchema_LiPAnalyzerOut'.
#' @param spectroTrP Spectronaut report of TrP data. Spectronaut report has to
#' be exported using the Spectronaut schema 'SpectroSchema_LiPAnalyzerOut'.
#' If LiP and TrP data were processed together in Spectronaut, please separate
#' them into two input files.
#' @param analysisLvl A character string defining the level on which the
#' peptide/protein quantities should be exported. Is set to 'Peptide' by
#' default, can alternatively be set to 'ModifiedPeptide' or 'Precursor'.This
#' will define the column which will be converted into row.names of the
#' matrices ('Peptide'='PEP.Quantity', ModifiedPeptide'='FG.Quantity',
#' 'Precursor'='FG.Quantity'). If you are want to choose a different column,
#'  else, please use the function /code{extractMSData}.
#'@param sampleName A character string or numeric giving the column name or
#' column number in which sample names are provided in the Spectronaut report.
#' Set to R.FileName' by default.
#' @param LiPonly A boolean value, set to 'TRUE' if you want to run the LiPonly
#' version o the package and not use trypsin-only data. Default is set to
#' 'FALSE'.
#'
#' @return
#' If run in default, a list of three matrices is returned ('quantityList'):
#' \itemize{
#'   \item 'LiPPep': LiP peptide quantities (or alternatively the modified
#'   peptide/ precursor/other quantities, dependent on \code{quantName} and
#'   \code{quantValue})
#'   \item 'TrPPep': TrP peptide quantities (or alternatively the modified
#'   peptide/precursor/other quantities, dependent on \code{quantName} and
#'   \code{quantValue})
#'   \item TrP protein quantities
#'   }
#'
#' If run in \code{LiPonly} mode, a list of two matrices is returned:
#' \itemize{
#'   \item 'LiPPep': LiP peptide quantities (or alternatively the modified
#'   peptide/precursor/other quantities, dependent on \code{quantName} and
#'   \code{quantValue})
#'   \item 'LiPProt': LiP protein quantities
#' }
#'
#' Rows represent features and columns samples.
#'
#' @export

extractSpectroData <- function(spectroLiP, spectroTrP=NULL,
                               analysisLvl="Peptide", sampleName="R.FileName",
                               LiPonly=FALSE){

    ## checking if value of 'AnalysisLvl' is set correctly
    if(!tolower(analysisLvl) %in% c("peptide", "modifiedpeptide", "precursor")){
        stop(analysisLvl, " is an unknown input for 'analysisLvl'. Please choose
'Peptide' or 'modifiedPeptide' or 'Precursor'.")
    }

    ## defining names feature and quantity columns based on 'AnalysisLvl'

    if(tolower(analysisLvl) == "peptide"){
        quantName <- "PEP.StrippedSequence"
        quantValue <- "PEP.Quantity"
    }
    else if(tolower(analysisLvl) == "modifiedpeptide"){
        quantName <- "EG.ModifiedSequence"
        quantValue <- "FG.Quantity"
    }
    else if(tolower(analysisLvl) == "precursor"){
        quantName <- "EG.PrecursorId"
        quantValue <- "FG.Quantity"
    }
    protValue <- "PG.Quantity"

    ## extracting feature matrices from reports of LiP (and TrP) data
    quantityList <- extractMSData(spectroLiP, spectroTrP, sampleName, quantName,
                                  quantValue, protValue, LiPonly)

    return(quantityList)
}


#' @title Extracting a single peptide/protein matrix from a MS report.
#'
#' @description Export quantities from the MS report, writing them into a
#' matrix were rows represent features and columns samples.
#'
#' @usage convert2Matrix(reportOut, quantValue, rows, cols)
#'
#' @param reportOut MS report exported including sample and quantity
#' information, e.g. Spectronaut report.
#' @param quantValue A character string or numeric variable referring to the
#' column name or column number from the MS report in which quantities or
#' interest are provided. If there are multiple quantities referring to the same
#' \code{rows}, the mean of these will be estimated and returned.
#' @param rows A character string or numeric giving the column name or column
#' number in the  MS report where quantity IDs are provided.
#' @param cols A character string or numeric giving the column name or column
#' number in the MS report where sample names are provided.
#'
#' @return Returns a matrix with quantities, \code{rows} represent features
#' and \code{cols} refer to the samples.
#'
#' @export

convert2Matrix <- function(reportOut, quantValue, rows, cols){
    message(quantValue, " used as intensities.")
    myMat <- data.frame(rows=gsub("_", "", reportOut[, rows]),
                        cols=reportOut[, cols],
                        values=reportOut[, quantValue])
    myMat <- reshape2::dcast(data=myMat,
                             formula=rows ~ cols,
                             fun.aggregate=mean,
                             fill=0, # adding 0 if no quantity given
                             value.var="values")

    ## setting all not measured values to NA
    myMat[myMat == 0] <- NA
    myMat[myMat == "NaN"] <- NA

    row.names(myMat) <- gsub(" ", "", myMat[,1])
    myMat <- myMat[,-1]
    return(myMat)
}


#' @title Creating an annotation data.frame with important peptide and protein
#' features.
#'
#' @description Create an annotation data.frame for the quantities from a MS
#' report/MS reports. This annotation file is needed for further downstream
#' analysis and provides a structured overview about peptide and protein
#' features.
#'
#' @usage getPepProtAnnot(reportOut, reportOut2=NULL, quantName, pepName,
#' protName, isTryptic, startPosition, allProtName=NULL, NMissedCleavages=NULL,
#' isProteotypic=NULL)
#'
#' @param reportOut MS report exported including sample and quantity
#' information, e.g. Spectronaut report.
#' The MS report has to include the following columns:
#' \itemize{
#'   \item names/identifiers of the quantity to be analyzes (e.g. peptide,
#'   modified peptide, precursor)
#'   \item the peptide names (i.e., AS sequence). This may be identical to the
#'   names/identifiers of the quantity to be analyzed.
#'   \item the matching protein names
#'   \item the information the peptide is full- or half-tryptic
#'   \item the start position of the peptide in the protein sequence
#' }
#' Optional columns to include are:
#' \itemize{
#'   \item all proteins matching the sequence
#'   \item number of missed cleavages
#'   \item proteotypic state
#'   }
#' @param reportOut2 Optional additional MS report including sample and
#' quantity information, e.g. Spectronaut report. \code{reportOut2} should
#' contain the same columns as \code{reportOut}. If \code{reportOut2} is added,
#' the peptide & protein annotation file will also include quantities uniquely
#' measured in \code{reportOut}. Default is 'NULL', assuming only one MS report
#' being presented.
#' @param quantName A character string or numeric giving the column name or
#' column number in which the row names of value used in the quantity matrices
#' are provided. This  These variable should be identical to the
#' \code{quantName} provided in \code{extractMSData} and therefore match the
#' row.names of the LiP and/or TrP data matrices.
#' @param pepName A character string or numeric giving the column name or
#' column number in which peptide names (i.e., AA sequences) are provided.
#' @param protName A character string or numeric giving the column name or
#' column number in which protein names are provided.
#' @param isTryptic A character string or numeric giving the column name or
#' column number in which annotation of the digest type are provided.
#' @param startPosition A character string or numeric giving the column name or
#' column number in which start position from each peptide in its protein
#' sequence are provided. Different start positions for the same peptide should
#' be separated by ',' for peptides mapping to the same protein multiple times
#' and ';' for peptides mapping to multiple proteins.
#' @param allProtName A character string or numeric giving the column name or
#' column number in which all protein names mapping to a peptide sequence are
#' provided.
#' @param NMissedCleavages A character string or numeric giving the column name
#' or column number in which number of missed cleavages in each peptide are
#' provided.
#' @param isProteotypic A character string or numeric giving the column name or
#' column number in which annotation of the peptide is proteotypic are provided.
#'
#' @return Returns a data.frame including all necessary annotation on peptides
#' (or precursors) and the matching proteins.
#'
#' @export

getPepProtAnnot <- function(reportOut,
                            reportOut2=NULL,
                            quantName,
                            pepName,
                            protName,
                            isTryptic,
                            startPosition,
                            allProtName=NULL,
                            NMissedCleavages=NULL,
                            isProteotypic=NULL){

    ## join Spectronaut reports if two are provided
    if(!is.null(reportOut2)){
        reportOut <-rbind(reportOut, reportOut2)
    }

    ## creating basic peptide & protein annotation file
    annotPP <- data.frame(quantID=gsub(" ", "", reportOut[, quantName]),
                          Peptide=gsub(" ", "", reportOut[, pepName]),
                          Protein=gsub(" ", "", reportOut[, protName]),
                          isTryptic=reportOut[, isTryptic])

    if(!is.null(allProtName)){
        annotPP$AllProtein <- reportOut[, allProtName]
        annotPP <- annotPP[, c(1,2,3,5,4)]
    }
    if(!is.null(NMissedCleavages)){
        annotPP$NMissedCleavages <- reportOut[, NMissedCleavages]
    }
    if(!is.null(isProteotypic)){
        annotPP$isProteotypic <- reportOut[, isProteotypic]
    }

    annotPP$startPosition <- reportOut[, startPosition]
    annotPP <- unique(annotPP) ## removing all duplicated rows from data.frame

    ## join protein names if quantity of interest was matched to different PG groups in
    ## different MS reports
    if(!is.null(reportOut2)){
        annotPP <- joinPG(annotPP, "quantID")
    }
    ## get end position of peptides
    annotPP$endPosition <- getEndPositionOfPep(annotPP)

    row.names(annotPP) <- annotPP[, "quantID"]
    return(annotPP)
}

#' @title Creating a data.frame with peptide and protein annotation.
#'
#' @description Create an annotation data.frame for the quantities from a
#' Spectronaut report/Spectronaut reports. The here created annotation file is
#' provides a structured overview about peptide and protein features and is
#' required in further downstream analysis.
#' Data must have been exported using the Spectronaut schema
#' 'SpectroSchema_LiPAnalyzerOut', else, please use the function
#' \code{extractMSData}.
#'
#' @usage getPepProtAnnotSpectro(spectroOut, spectroOut2=NULL,
#' analysisLvl="Peptide")
#'
#' @param spectroOut MS report exported including sample and quantity information,
#' e.g. Spectronaut report.
#' @param spectroOut2 An optional additional MS report including sample and
#' quantity information, e.g. Spectronaut report. If included theh peptide &
#' protein annotation file will also include quantities uniquely measured in
#' \code{reportOut}. Default is 'NULL', assuming only one MS report being
#' presented.
#' @param analysisLvl A character string defining the level on which the
#' peptide/protein quantities should be exported. Peptide' by default, can
#' alternatively be set to 'ModifiedPeptide' or Precursor'. Should be identical
#' to the \code{analysisLvl} defined in \code{extractSpectroData}.
#'
#' @return Returns a data.frame including all necessary annotation on peptides
#' and proteins
#'
#' @export

getPepProtAnnotSpectro <- function(spectroOut,
                                   spectroOut2=NULL,
                                   analysisLvl="Peptide"){

    ## checking if value of AnalysisLvl is set correctly
    if(!tolower(analysisLvl) %in% c("peptide", "modifiedpeptide", "precursor")){
        stop(analysisLvl, " is an unknown input for 'analysisLvl'. Please choose
'peptide' or 'modifiedpeptide' or 'precursor'.")
    }
    ## assigning quantName accordingly to choosen analysis level
    quantName <- switch(tolower(analysisLvl),
                        peptide = "PEP.StrippedSequence",
                        modifiedpeptide = "EG.ModifiedSequence",
                        precursor = "EG.PrecursorId")

    ## call getPepProtAnnot with specified Spectronaut names
    annotPP <- getPepProtAnnot(spectroOut, spectroOut2, quantName,
                               "PEP.StrippedSequence", "PG.ProteinGroups",
                               "PEP.DigestType....Trypsin.P.",
                               "PEP.PeptidePosition",
                               "PEP.AllOccurringProteinAccessions",
                               "PEP.NrOfMissedCleavages","PEP.IsProteotypic")
    return(annotPP)
}

#' @title Estimating end position of peptide sequences
#'
#' @description Getting end position of peptides sequences from the starting
#' positions provided. Function takes into account that peptides can map to
#' multiple regions in the proteome, providing multiple end positions if
#' multiple start positions are given.
#'
#' @usage getEndPositionOfPep(annotPP, startPosition="startPosition",
#' pepCol="Peptide")
#'
#' @param annotPP A data.frame with peptide and protein annotation were rows are
#' features. The data.frame must include a column providing the peptide
#' sequences and one providing the start position the peptide has in protein
#' sequence. Different start positions for the same peptide should be separated
#' by ',' for peptides mapping to the same protein multiple times and ';' for
#' peptides mapping to multiple proteins.
#' @param startPosition A character string or numeric giving the column name or
#' column number defining column in which start positions of peptides are
#' provided.
#' Default is 'startPosition'.
#' @param pepCol A character string or numeric giving the column name or
#' column number defining column in which peptide sequences are provide. Peptide
#' sequences have to be provided without modifications or charges.
#' Default is 'Peptide'.
#'
#' @return Returns a character vector with the end positions of all peptides.
#' Different start positions for the same peptide should be separated by ',' for
#' peptides mapping to the same protein multiple times and ';' for peptides
#' mapping to multiple proteins.

getEndPositionOfPep <- function(annotPP, startPosition="startPosition",
                                pepCol="Peptide"){
    end <- sapply(seq(1:nrow(annotPP)), function(i){
        start <- as.character(annotPP[, startPosition][i])
        length <- nchar(annotPP[, pepCol][i])-1

        ## estimate end position, if only one start position is provided
        if(!grepl("[;,]", start)){
            end <- as.character(as.numeric(start)+length)
        }

        ## split start positions and estimate end position if multiple start
        ## positions of peptide are provided
        else{
            start <- unlist(strsplit(start, ";"))
            end <- sapply(start, function(x){
                x <- unlist(strsplit(x, ","))
                end <- paste(unname(sapply(x, function(y){as.character(
                    as.numeric(y)+length)})), collapse=",")
            })
            end <- paste(unname(end), collapse=";")
        }
        return(end)
    })

    return(end)
}

#' @title Combining protein names of identical peptides in the peptide and
#' protein annotation data.frame.
#'
#' @description Joining protein names together in the respective column in the
#' peptide and protein annotation data.frame in cases were a peptide was matched
#' to different PG groups in different MS reports. Only relevant if two MS
#' reports are used to create the annotation data.frame.
#'
#' @param annotPP data.frame with peptide and protein annotation were rows are
#' features. The data.frame must include a column with the identifiers for the
#' quantity of interest, which should be unique after this function is run, e.g.
#' peptide, modified peptide or precursor ID. It mus additionally include a
#' column providing the matching protein names. Protein names for the same
#' peptide should be separated by ',' if a peptide maps to the same protein
#' multiple times and ';' if a peptide maps to multiple proteins.
#' @param iCol A character string or numeric giving the column name or
#' column number containing the identifiers for the quantity of interest.
#' @param protCol A character string or numeric giving the column name or
#' column number in which protein names  are provide. Set to 'Protein' by
#' default.
#'
#' @return A data.frame were the chosen quantity of interest, e.g. peptide,
#' modified peptide or precursor ID, is unique.

joinPG <- function(annotPP, iCol, protCol="Protein"){
    annotPP <- split(annotPP, annotPP[, iCol])
    annotPP <- do.call(rbind, lapply(annotPP, function(x){
        if(nrow(x)>1){
            Proteins <- unlist(strsplit(paste(x[, protCol], collapse=";"), ";"))
            Proteins <- paste(Proteins[!duplicated(Proteins)], collapse=";")
            x <- x[1, , drop=FALSE]
            x [, protCol] <- Proteins
        }
        return(x)
    }))
    annotPP <- unique(annotPP)
    return(annotPP)
}

#' @title Creates data.frame with sample annotation from MS report
#'
#' @description Create an annotation file on samples including information on
#' the conditions of samples. This function only runs, if all the annotations
#' were added to the MS report, else please create the sample annotation file
#' yourself. The file should include all variables which should be used in the
#' LiPAnalyzeR models. Rows should be samples, with row names having the
#' matching row.names to the column names in the quantity matrices. Columns
#' should contain the necessary sample information, for example a 'Condition'
#' column providing the sample-specific condition or a 'Batch' column.
#'
#' @usage getSampleAnnot(reportOut, sampleName="R.FileName",
#' sampleCondition="R.Condition", typeCondition="factor",
#' contrastCoding="dummyCoding", baseLevel=NULL)
#'
#' @param reportOut Spectronaut report of MS data. Spectronaut report should
#' have been exported using the Spectronaut schema
#' 'SpectroSchema_LiPAnalyzerOut'.
#' @param sampleName A character string or numeric giving the column name or
#' column number providing the sample names in the MS report.
#' Default is 'R.FileName'.
#' @param sampleCondition A character string or numeric giving the column name
#' or column number providing sample conditions in the MS report.
#' Default is 'R.Condition'.
#' @param typeCondition A character string providing information if condition is
#' a 'factor' or 'continuous' variable. If condition is 'continuous', the column
#' containing the condition in the MS report has to be numeric variable. If
#' the condition is a 'factor', contrast coding for later models is defined
#' within this function.
#' Default is 'factor'.
#' @param contrastCoding A character string providing information which type of
#' contrast coding to use if \code{typeCondition} is 'factor'. Can be set to
#' 'dummyCoding', 'sumCoding' or 'wecCoding'. If you do not want to define
#' the contrast coding method within this function set to NULL and subsequently
#' define it before running the contrast model.
#' Default is 'dummyCoding'.
#' @param baseLevel A character strong giving name of reference level if
#' \code{typeCondition} is a factor. If not provided the function will choose
#' the first occurring condition in the MS report as the reference level.
#'
#' @return Returns a data.frame containing sample and condition annotation were
#' rows are samples.
#'
#' @export

getSampleAnnot <- function(reportOut,
                           sampleName="R.FileName",
                           sampleCondition="R.Condition",
                           typeCondition="factor",
                           contrastCoding="dummyCoding",
                           baseLevel=NULL){

    ## check if typeCondition input is as expected
    if(!tolower(typeCondition) %in% c("factor", "continuous")){
        stop("Please set 'typeCondition' to 'factor' or 'continuous.")
    }

    ## check if contrastCoding input is as expected

    if(!is.null(contrastCoding)){
        if(!tolower(contrastCoding) %in% c("dummycoding", "sumcoding",
                                           "weccoding")){
            stop("Please set 'contrastCoding' to 'dummycoding', 'sumcoding' or
'weccoding'. If you do not want to set the contrast method yet, set to 'NULL'.")
        }

    }

    ## Factorial condition settings
    if(tolower(typeCondition) == "factor"){
        message("Setting condition variable to a factor variable.")
        Condition <- as.factor(as.character(reportOut[, sampleCondition]))

        if(is.null(contrastCoding)){
            message("Not defining contrast coding. Please define contrast coding
method yourself before running the contrast model, else the function will
automatically choose dummy coding and the condition first in the alphabetical
order as the reference level.")
        }

        ## Applying coding methods to factorial condition
        else{
            if(is.null(baseLevel)){
                message("No 'baseLevel' provided, choosing the first occuring
condition in the MS report file as the reference level.")
                baseLevel <- as.character(Condition[1])
            }
            NbaseLevel <- which(levels(Condition) == baseLevel)

            if(tolower(contrastCoding) == "dummycoding"){
                message("Applying dummy coding to condition variable.")
                stats::contrasts(Condition) <- stats::contr.treatment(nlevels(
                    Condition), base=which(levels(Condition) == baseLevel))

            }
            else if(tolower(contrastCoding) == "sumcoding"){
                message("Applying sum coding to condition variable..")
                stats::contrasts(Condition) <- stats::contr.sum(Condition)
            }
            else if(tolower(contrastCoding) == "weccoding"){
                message("Applying wec coding to condition variable.")
                stats::contrasts(Condition) <- wec::contr.wec(Condition,
                                                              omitted=baseLevel)
            }
        }
    }

    ## Continious condition settings
    else if(tolower(typeCondition) == "continuous"){
        message("Setting condition variable to a continuous variable.")
        Condition <- as.numeric(reportOut[, sampleCondition])
    }

    annotS <- data.frame(SampleName=reportOut[, sampleName],
                        Condition= Condition)
    annotS <- unique(annotS)
    row.names(annotS) <- annotS$SampleName

    return(annotS)
}
