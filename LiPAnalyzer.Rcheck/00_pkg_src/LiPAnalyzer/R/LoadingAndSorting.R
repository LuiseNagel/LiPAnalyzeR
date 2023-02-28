
#' @title Creating  peptide and protein matrices from the Spectronaut report.
#'
#' @description Function which exports peptide and protein quantities from a Spectronaut report and writes them into a matrix.
#' Data is exported from both the LiP and Trp Spectronaut reports per default.
#' Alternatively, function can be run in \code{LiPonly} mode, only requiring LiP data.
#'
#' @usage ExtractDataFromSpectro(spectroLiP, spectroTrp = NULL, analysisLvl = "Peptide",
#' sampleName = "R.FileName", valuePep = NULL,
#' valueProt = "PG.Quantity", LiPonly = F)
#'
#' @param spectroLiP Spectronaut report of LiP data.
#' Spectronaut report should have been exported using the Spectronaut schema SpectroSchema_LiPAnalyzerOut.
#' @param spectroTrp Spectronaut report of Trp data
#' Spectronaut report should have been exportet using the Spectronaut schema SpectroSchema_LiPAnalyzerOut.
#' @param analysisLvl a character defining the level on which the peptide/protein quantities should be exported.
#' 'Peptide' by default, can alternatively be set to "ModifiedPeptide" or "Precursor".
#' @param sampleName a character defining the column of the Spectronaut report where samples names are provided.
#' 'R.Filename' by default.
#' @param valuePep a character defining the column from which peptide/precursor quantities should be taken from.
#' NULL by default, automatically setting \code{valuePep} based on the \code{AnalysisLvl}.
#' @param valueProt a character defining the column from which protein quantities should be taken from.
#' 'PG.Quantity' by default.
#' @param LiPonly a boolean value, set to 'TRUE' if you are running the
#' LiPonly version o the package and not providing trypsin-only data
#'
#' @return list of three matrices, containing the LiP peptide/precursor,
#' Trp peptide/precursor and Trp protein quantities.
#' Rows represent features and columns refer to the samples.
#'
#' @export

ExtractDataFromSpectro <- function(spectroLiP, spectroTrp = NULL, analysisLvl = "Peptide",
                                   sampleName = "R.FileName", valuePep = NULL,
                                   valueProt = "PG.Quantity", LiPonly = F){

    # checking if value of AnalysisLvl is set correctly
    if(!tolower(analysisLvl) %in% c("peptide", "modifiedpeptide", "precursor")){
        stop(analysisLvl, " is an unknown input for 'analysisLvl'. Please choose 'peptide' or 'modifiedpeptide' or 'precursor'.")
    }

    # checking if spectroTrp is [provided in case that LiPonly is False
    if(!LiPonly & is.null(spectroTrp)){
        stop("Please add spectroTrp data or set LiPonly to TRUE.")
    }

    # Setting names feature and quantity columns based on AnalysisLvl
    if(tolower(analysisLvl) == "peptide"){
        rows <- "PEP.StrippedSequence"
        if(is.null(valuePep)){
            valuePep <- "PEP.Quantity"
        }
    }
    else if(tolower(analysisLvl) == "modifiedpeptide"){
        rows <- "EG.ModifiedSequence"
        if(is.null(valuePep)){
            valuePep <- "FG.Quantity"
        }
    }
    else if(tolower(analysisLvl) == "precursor"){
        rows <- "EG.PrecursorId"
        if(is.null(valuePep)){
            valuePep <- "FG.Quantity"
        }
    }

    # Writing feature matrices for LiP and Trp data
    message("Creating LiP ", tolower(analysisLvl), " matrix.")
    LiPPep <- Convert2Matrix(spectroLiP, valuePep, rows, sampleName)
    if(LiPonly){
        message("Creating LiP protein matrix.")
        LiPProt <- Convert2Matrix(spectroLiP, valueProt, rows, sampleName)
        out <- list(LiPPep = LiPPep, LiPProt = LiPProt)
    }
    else{
        message("Creating Trp ", tolower(analysisLvl), " matrix.")
        TrpPep <- Convert2Matrix(spectroTrp, valuePep, rows, sampleName)
        message("Creating Trp protein matrix.")
        TrpProt <- Convert2Matrix(spectroTrp, valueProt, rows, sampleName)
        out <- list(LiPPep = LiPPep, TrpPep = TrpPep, TrpProt = TrpProt)
    }

    return(out)
}


#' @title Extracting a single matrix from Spectronaut report
#'
#' @description Function which exports quantities from the Spectronaut report and writes them into a matrix.
#' Rows are peptides/proteins and columns are samples.
#'
#' @usage Convert2Matrix(spectroOut, values, rows, cols)
#'
#' @param spectroOut Spectronaut report
#' @param values character defining column in Spectronaut report where quantities to export are located
#' @param rows character defining column in Spectronaut report where the names of the quantities to export are located.
#' If there are different quantities for the same \code{rows} the mean of these values is used.
#' @param cols character defining  in Spectronaut report where the sample or file name are located
#'
#' @return matrix with quantities, \code{rows} represent features and \code{cols} refer to the samples.
#'
#' @export

Convert2Matrix <- function(spectroOut, values, rows, cols){
    message(values, " used as intensities.")
    myMat <- data.frame(rows = gsub("_", "", spectroOut[, rows]),
                        cols = spectroOut[, cols],
                        values = spectroOut[, values])
    myMat <- reshape2::dcast(data = myMat,
                             formula = rows ~ cols,
                             fun.aggregate = mean,
                             fill = 0, # adding 0 if no quantity included in Spectronaut report
                             value.var = "values")

    # Setting all not measured values to NA in matrix
    myMat[myMat == 0] <- NA
    myMat[myMat == "NaN"] <- NA

    rownames(myMat) <- gsub(" ", "", myMat[,1])
    myMat <- myMat[,-1]
    return(myMat)
}


#' @title Creating data.frame with peptide and protein annotation
#'
#' @description Creates an annotation file on the peptides and proteins in the Spectronaut report.
#'
#' @usage GetPepProtAnnot(spectroOut, spectroOut2 = NULL,analysisLvl = "Peptide",
#' Precursor = "EG.PrecursorId", modPeptide = "EG.ModifiedSequence", Peptide = "PEP.StrippedSequence",
#' Protein = "PG.ProteinGroups", AllProtein = "PEP.AllOccurringProteinAccessions", NMissedCleavages = "PEP.NrOfMissedCleavages",
#' isProteotypic = "PEP.IsProteotypic", isTryptic = "PEP.DigestType....Trypsin.P.", start = "PEP.PeptidePosition")
#'
#' @param spectroOut Spectronaut report of MS data.
#' Spectronaut report should have been exported using the Spectronaut schema SpectroSchema_LiPAnalyzerOut.
#' @param spectroOut2 optional additional Spectronaut report, default is NULL.
#' Spectronaut report should have been exported using the Spectronaut schema SpectroSchema_LiPAnalyzerOut.
#' @param analysisLvl character defining the level on which the peptide/protein quantities should be exported.
#' 'Peptide' by default, can alternatively be set to 'ModifiedPeptide' or 'Precursor'.
#' Should be identical to how \code{analysisLvl} defined in ExtractDataFromSpectro.
#' @param Precursor character giving the column name in which precursor IDs can be found in the Spectronaut report.
#' 'EG.PrecursorId' by default.
#' Only used if \code{analysisLvl} is set to 'Precursor'.
#' @param modPeptide character giving the column in which modified peptide sequence can be found in the Spectronaut report.
#' 'EG.ModifiedSequence' by default.
#' Only used if \code{analysisLvl} is set to 'ModifiedPeptide'.
#' @param Peptide character giving the column in which peptide sequence can be found in the Spectronaut report.
#' 'PEP.StrippedSequence' by default.
#' @param Protein character giving the column in which protein names can be found in the Spectronaut report.
#' 'PG.ProteinGroups' by default.
#' @param AllProtein character giving the column in which all protein names mapping to a peptide sequence can be found in the Spectronaut report.
#' 'PEP.AllOccurringProteinAccessions' by default.
#' @param NMissedCleavages character giving the column in which number of missed cleavages in each peptide can be found in the Spectronaut report.
#' 'PEP.NrOfMissedCleavages' by default.
#' @param isProteotypic character giving the column in which annotation if the peptide is proteotypic can be found in the Spectronaut report.
#' 'PEP.IsProteotypic' by default.
#' @param isTryptic character giving the column in which annotation of the digest type can be found in the Spectronaut report.
#' 'PEP.DigestType....Trypsin.P.' by default.
#' @param start character giving the column in which start position in each peptide can be found in the Spectronaut report.
#' 'PEP.PeptidePosition' by default.
#'
#'  @return data.frame including all necessary annotation on peptides and proteins
#'
#' @export

GetPepProtAnnot <- function(spectroOut,
                            spectroOut2 = NULL,
                            analysisLvl = "Peptide",
                            Precursor = "EG.PrecursorId",
                            modPeptide = "EG.ModifiedSequence",
                            Peptide = "PEP.StrippedSequence",
                            Protein = "PG.ProteinGroups",
                            AllProtein = "PEP.AllOccurringProteinAccessions",
                            NMissedCleavages = "PEP.NrOfMissedCleavages",
                            isProteotypic = "PEP.IsProteotypic",
                            isTryptic = "PEP.DigestType....Trypsin.P.",
                            start = "PEP.PeptidePosition"){

    # checking if value of AnalysisLvl is set correctly
    if(!tolower(analysisLvl) %in% c("peptide", "modifiedpeptide", "precursor")){
        stop(analysisLvl, " is an unknown input for 'analysisLvl'. Please choose 'peptide' or 'modifiedpeptide' or 'precursor'.")
    }

    if(!is.null(spectroOut2)){
        spectroOut <-rbind(spectroOut, spectroOut2) # join spectronaut reports if two are provided
    }

    # creating basic peptide & protein annotation file
    annotPP <- data.frame(Peptide = spectroOut[, Peptide],
                          Protein = spectroOut[, Protein],
                          AllProtein = spectroOut[, AllProtein],
                          NMissedCleavages = spectroOut[, NMissedCleavages],
                          isProteotypic = spectroOut[, isProteotypic],
                          isTryptic  = spectroOut[, isTryptic],
                          startPosition = spectroOut[, start])

    # adding modified peptides or precursors if a different analysisLvl is chosen
    if(tolower(analysisLvl) == "modifiedpeptide"){
        annotPP$modPeptide <- gsub("[ _]", "", spectroOut[, modPeptide])
        annotPP <- annotPP[, c(ncol(annotPP), 1:ncol(annotPP)-1)]
    }
    if(tolower(analysisLvl) == "precursor"){
        annotPP$Precursor <- gsub("[ _]", "", spectroOut[, Precursor])
        annotPP <- annotPP[, c(ncol(annotPP), 1:ncol(annotPP)-1)]
    }

    annotPP <- unique(annotPP) # removing all duplicated rows from data.frame
    iRow <- switch(tolower(analysisLvl), "peptide" = "Peptide", "modifiedpeptide" = "modPeptide", "precursor" = "Precursor")

    # join protein names if peptide was matched to different PG groups in different spectronaut reports
    if(!is.null(spectroOut2)){
        annotPP <- JoinPG(annotPP, iRow)
    }
    annotPP$endPosition <- GetEndPositionOfPep(annotPP) # get end position of peptides

    row.names(annotPP) <- annotPP[, iRow]
    return(annotPP)
}

#' @title Estimating end position of peptide sequences
#'
#' @description Getting end position fo peptides sequences from the starting positions provided.
#' Function takes into account that peptides can map to multiple regions in the proteome.
#'
#' @usage GetEndPositionOfPep(annotPP, startCol = "startPosition", pepCol = "Peptide")
#'
#' @param annotPP a data.frame with peptide and protein annotation.
#' Rows are features. Data frame must include a column providing the start position of each peptide in the proteins
#' and one containing the peptide sequence peptides are provided.
#' Different start positions for the same peptide should be separated by
#' ',' for peptides mapping to the same protein multiple times and
#' ';' for peptides mapping to multiple proteins.
#' @param startCol a character defining column in which start positions of peptides are provided
#' 'startPosition' by default.
#' @param pepCol character defining column in which peptide sequences are provide.
#' Peptide sequences have to be provided without modifications or charges.
#' 'Peptide' by default.
#'
#' @return character vector of end positions of peptides
#' #' Different start positions for the same peptide should be separated by
#' ',' for peptides mapping to the same protein multiple times and
#' ';' for peptides mapping to multiple proteins.

GetEndPositionOfPep <- function(annotPP, startCol = "startPosition", pepCol = "Peptide"){
    end <- sapply(seq(1:nrow(annotPP)), function(i){
        start <- as.character(annotPP[, startCol][i])
        length <- nchar(annotPP[, pepCol][i])-1

        # estimate end position if  one start position of peptide is provided
        if(!grepl("[;,]", start)){
            end <- as.character(as.numeric(start)+length)
        }

        # split start positions and estimate end position if multiple start positions of peptide are provided
        else{
            start <- unlist(strsplit(start, ";"))
            end <- sapply(start, function(x){
                x <- unlist(strsplit(x, ","))
                end <- paste(unname(sapply(x, function(y){as.character(as.numeric(y)+length)})), collapse = ",")
            })
            end <- paste(unname(end), collapse = ";")
        }

        return(end)
    })

    return(end)
}

#' @title Combining protein names of identical peptides
#'
#' @description Joining protein names together if a peptide was matched to different PG groups in different Spectronaut reports.
#'
#' @param annotPP data.frame with peptide and protein annotation.
#' Rows are features, must include a column protein name(s)
#' @param iRow character describing row in which features are provided
#' Based on the chosen \code{analysisLvl}.
#' @param protCol character defining column in which protein name(s) are provide.
#' 'Protein' by default.
#'
#' @return data.frame including all necessary annotation on peptides and proteins with joined protein names

JoinPG <- function(annotPP, iRow, protCol = "Protein"){
    annotPP <- split(annotPP, annotPP[, iRow])
    annotPP <- do.call(rbind, lapply(annotPP, function(x){
        if(nrow(x)>1){
            Proteins <- unlist(strsplit(paste(x[, protCol], collapse = ";"), ";"))
            Proteins <- paste(Proteins[!duplicated(Proteins)], collapse = ";")
            x <- x[1, , drop = F]
            x [, protCol] <- Proteins
        }
        return(x)
    }))
    annotPP <- unique(annotPP)
    return(annotPP)
}

#' @title Creating data.frame with sample annotation from Spectronaut report
#'
#' @description Creates an annotation file on samples including information on the conditions of samples.
#' This function only runs, if all the annotations were added to the Spectronaut report, else
#' please create the sample annotation file yourself.
#'
#' @usage GetSampleAnnot(spectroOut, sampleName = "R.FileName", sampleCondition = "R.Condition",
#' typeCondition = "factor", contrastCoding = "dummyCoding", baseLevel = NULL)
#'
#' @param spectroOut Spectronaut report.
#' Spectronaut report should have been exported using the Spectronaut schema SpectroSchema_LiPAnalyzerOut.
#' @param sampleName character giving the column in which sample names can be found in the Spectronaut report.
#' 'R.FileName' by default.
#' Only used if \code{analysisLvl} is set to 'Precursor'.
#' @param sampleCondition character giving the column in which conditions which should be compared can be found in the Spectronaut report.
#' 'R.Condition' by default.
#' @param typeCondition character providing information if condition is a factor or continuous variable
#' If condition is continuous has to be provided as numeric.
#' 'factor' is default.
#' @param contrastCoding character providing information which type of contrast coding to use if condition is a factor
#' 'dummyCoding' is default, alternatives 'sumCoding' or 'wecCoding' can be chosen
#' Can be set to NULL to not define contrast coding automatically
#' @param baseLevel character giving name of reference level if condition is a factor
#' and dummy or wec coding is used
#' If not provided will choose the first condition in Spectronaut report
#'
#' @return data.frame containing sample and condition annotation.
#' Rows are samples.
#'
#' @export

GetSampleAnnot <- function(spectroOut,
                           sampleName = "R.FileName",
                           sampleCondition = "R.Condition",
                           typeCondition = "factor",
                           contrastCoding = "dummyCoding",
                           baseLevel = NULL){
    if(tolower(typeCondition) == "factor"){
        Condition <- as.factor(as.character(spectroOut[, sampleCondition]))
        if(is.null(baseLevel)){
            baseLevel <- as.character(Condition[1])
        }
        NbaseLevel <- which(levels(Condition) == baseLevel)
        if(is.null(contrastCoding)){
            message("Not defining contrast coding.")
        }
        else if(tolower(contrastCoding) == "dummycoding"){
            message("Using dummy coding for condition.")
            stats::contrasts(Condition) <- stats::contr.treatment(nlevels(Condition),
                                                                  base = which(levels(Condition) == baseLevel))

        }
        else if(tolower(contrastCoding) == "sumcoding"){
            message("Using sum coding for condition.")
            stats::contrasts(Condition) <- stats::contr.sum(Condition)
        }
        else if(tolower(contrastCoding) == "weccoding"){
            message("Using wec coding for condition.")
            stats::contrasts(Condition) <- wec::contr.wec(Condition,
                                                         omitted = baseLevel)
        }
    }
    else if(tolower(typeCondition) == "continuous"){
        Condition <- as.numeric(spectroOut[, sampleCondition])
    }
    annotS <- data.frame(SampleName = spectroOut[, sampleName],
                        Condition =  Condition)
    annotS <- unique(annotS)
    row.names(annotS) <- annotS$SampleName
    return(annotS)
}
