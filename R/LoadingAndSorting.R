
#' @title Creating peptide and protein matrices from the Spectronaut report.
#'
#' @description This function writes peptide and protein quantities from a
#' Spectronaut report into a matrix.Data is exported from both the LiP and Trp
#' Spectronaut reports per default. Alternatively, function can be run in
#' \code{LiPonly} mode, only requiring LiP data.
#'
#' @usage extractDataFromSpectro(spectroLiP, spectroTrp=NULL,
#' analysisLvl="Peptide", rowName=NULL, sampleName="R.FileName", valuePep=NULL,
#' valueProt="PG.Quantity", LiPonly=FALSE)
#'
#' @param spectroLiP Spectronaut report of LiP data. Spectronaut report should
#' have been exported using the Spectronaut schema
#' 'SpectroSchema_LiPAnalyzerOut'.
#' @param spectroTrp Spectronaut report of Trp data. Spectronaut report should
#' have been exported using the Spectronaut schema
#' 'SpectroSchema_LiPAnalyzerOut'.
#' @param analysisLvl A character string defining the level on which the
#' peptide/protein quantities should be exported. Is set to 'Peptide' by
#' default, can alternatively be set to 'ModifiedPeptide' or 'Precursor'.This
#' will define the column which will be converted into row.names of the
#' matrices ('Peptide'='PEP.Quantity', ModifiedPeptide'='FG.Quantity',
#' 'Precursor'='FG.Quantity'). If you are using an alternative input format
#' and/or want to choose a different column, please set this variable to
#' 'individual' and specify the column in the variable \code{rowName}. You then
#' also have to define the \code{'valuePep'}.
#' @param rowName A character string or numeric giving the column name or column
#' number in which the row names of the quantity matrixes are defined. This
#' parameter will only be evaluated if \code{analysisLvl} is set to
#' individual'.
#' @param sampleName A character string or numeric giving the column name or
#' column number in which sample names are defined in the Spectronaut report.
#' Set to R.FileName' by default.
#' @param valuePep A character string or numeric giving the column name or
#' column number from which peptide/precursor quantities should be taken from.
#' Defined as 'NULL' by default, resulting in the function automatically setting
#' \code{valuePep} based on the \code{AnalysisLvl} ('Peptide'='PEP.Quantity',
#' 'ModifiedPeptide'='FG.Quantity', 'Precursor'='FG.Quantity').
#' @param valueProt A character string or numeric giving the column name or
#' column number from which protein quantities should be taken from. Set to
#' 'PG.Quantity' by default.
#' @param LiPonly A boolean value, set to 'TRUE' if you want to run the LiPonly
#' version o the package and not providing trypsin-only data.Default is set to
#' 'FALSE'.
#'
#' @return Per default a list of three matrices is returned, containing the LiP
#' peptide/precursor, Trp peptide/precursor and Trp protein quantities. If run
#' in \code{LiPonly} mode will return a list of two matrcies with the LiP
#' peptide/precursor and Trp protein quantities. Rows represent features and
#' columns refer to the samples.
#'
#' @export

extractDataFromSpectro <- function(spectroLiP, spectroTrp=NULL,
                                   analysisLvl="Peptide", rowName=NULL,
                                   sampleName="R.FileName", valuePep=NULL,
                                   valueProt="PG.Quantity", LiPonly=FALSE){

    ## checking if value of 'AnalysisLvl' is set correctly
    if(!tolower(analysisLvl) %in% c("peptide", "modifiedpeptide", "precursor",
                                    "individual")){
        stop(analysisLvl, " is an unknown input for 'analysisLvl'. Please choose
             'peptide' or 'modifiedpeptide' or 'precursor' or 'individual'.")
    }

    ## checking if 'spectroTrp' is provided in case that LiPonly is FALSE
    if(!LiPonly & is.null(spectroTrp)){
        stop("Please add spectroTrp data or set LiPonly to TRUE.")
    }

    ## Checking AnalysisLvl
    if(is.null(rowName) & tolower(analysisLvl) == "individual"){
        stop("'rowName' is defined but 'analysisLvl' is not set to
             'individual'. Either set rowName to NULL or set 'analysisLvl' to
             'individual'.")
        if(is.null(valuePep)){
            stop("Since you are running the analysis level 'individual',
                 please define 'valuePep'.")
        }
    }

    ## defining names feature and quantity columns based on 'AnalysisLvl'
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
    else if(tolower(analysisLvl) == "individual"){
        rows <- rowName
    }

    ## extracting feature matrices for LiP and Trp data
    message("Creating LiP ", tolower(analysisLvl), " matrix.")
    LiPPep <- convert2Matrix(spectroLiP, valuePep, rows, sampleName)

    if(LiPonly){
        message("Creating LiP protein matrix.")
        LiPProt <- convert2Matrix(spectroLiP, valueProt, rows, sampleName)
        out <- list(LiPPep=LiPPep, LiPProt=LiPProt)
    }

    else{
        message("Creating Trp ", tolower(analysisLvl), " matrix.")
        TrpPep <- convert2Matrix(spectroTrp, valuePep, rows, sampleName)
        message("Creating Trp protein matrix.")
        TrpProt <- convert2Matrix(spectroTrp, valueProt, rows, sampleName)
        out <- list(LiPPep=LiPPep, TrpPep=TrpPep, TrpProt=TrpProt)
    }

    return(out)
}


#' @title Extracting a single peptide/protein matrix from Spectronaut report
#'
#' @description This function exports quantities from the Spectronaut report
#' and writes them into a matrix. Rows are peptides/proteins and columns are
#' samples.
#'
#' @usage convert2Matrix(spectroOut, values, rows, cols)
#'
#' @param spectroOut Spectronaut report exported using the Spectronaut schema
#' 'SpectroSchema_LiPAnalyzerOut'.
#' @param values A character string or numeric giving the column name or
#' column number from the Spectronaut report in which quantities that should be
#' writen into matrix are located.If there are different quantities for the
#' same \code{rows} the mean of these values is used.
#' @param rows A character string or numeric giving the column name or
#' column number from the  Spectronaut report where the names of the quantities
#' writen into matrix are located.
#' @param cols A character string or numeric giving the column name or
#' column number from the  Spectronaut report where the sample or file names are
#' located.
#'
#' @return Returns a matrix with quantities, \code{rows} represent features
#' and \code{cols} refer to the samples.
#'
#' @export

convert2Matrix <- function(spectroOut, values, rows, cols){
    message(values, " used as intensities.")
    myMat <- data.frame(rows=gsub("_", "", spectroOut[, rows]),
                        cols=spectroOut[, cols],
                        values=spectroOut[, values])
    myMat <- reshape2::dcast(data=myMat,
                             formula=rows ~ cols,
                             fun.aggregate=mean,
                             fill=0, # adding 0 if no quantity given
                             value.var="values")

    ## setting all not measured values to NA
    myMat[myMat == 0] <- NA
    myMat[myMat == "NaN"] <- NA

    rownames(myMat) <- gsub(" ", "", myMat[,1])
    myMat <- myMat[,-1]
    return(myMat)
}


#' @title Creating a data.frame with peptide and protein annotation.
#'
#' @description Creates an annotation file on the peptides and proteins in the
#' Spectronaut report.
#'
#' @usage getPepProtAnnot(spectroOut, spectroOut2=NULL,analysisLvl="Peptide",
#' Precursor="EG.PrecursorId", modPeptide="EG.ModifiedSequence",
#' Peptide="PEP.StrippedSequence", Protein="PG.ProteinGroups",
#' AllProtein="PEP.AllOccurringProteinAccessions",
#' NMissedCleavages="PEP.NrOfMissedCleavages",
#' isProteotypic="PEP.IsProteotypic", isTryptic="PEP.DigestType....Trypsin.P.",
#' start="PEP.PeptidePosition")
#'
#' @param spectroOut Spectronaut report of MS data. Spectronaut report should
#' have been exported using the Spectronaut schema
#' 'SpectroSchema_LiPAnalyzerOut'.
#' @param spectroOut2 An optional additional Spectronaut report, default is
#' 'NULL', assuming only one Spectronaut report is present. Spectronaut
#' report(s) should have been exported using the Spectronaut schema
#' 'SpectroSchema_LiPAnalyzerOut'.
#' @param analysisLvl A character string defining the level on which the
#' peptide/protein quantities should be exported. Peptide' by default, can
#' alternatively be set to 'ModifiedPeptide' oR Precursor'. Should be identical
#' to the \code{analysisLvl} defined in \code{ExtractDataFromSpectro}.
#' @param Precursor  A character string giving the column name in which
#' precursor IDs can be found in the Spectronaut report. Defined as
#' 'EG.PrecursorId' by default. Only used if \code{analysisLvl} is set to
#' 'Precursor'.
#' @param modPeptide A character string or numeric giving the column name or
#' column number in which modified peptide sequence can be found in the
#' Spectronaut report. Defined as EG.ModifiedSequence' by default. Only used if
#' \code{analysisLvl} is set to 'ModifiedPeptide'.
#' @param Peptide A character string or numeric giving the column name or
#' column number in which peptide sequence can be found in the Spectronaut
#' report. Defined as PEP.StrippedSequence' by default.
#' @param Protein A character string or numeric giving the column name or
#' column number in which protein names can be found in the Spectronaut report.
#' Defined as 'PG.ProteinGroups' by default.
#' @param AllProtein A character string or numeric giving the column name or
#' column number in which all protein names mapping to a peptide sequence can be
#' found in the Spectronaut report. Defined as
#' 'PEP.AllOccurringProteinAccessions' by default.
#' @param NMissedCleavages A character string or numeric giving the column name or
#' column number in which number of missed cleavages in each peptide can be
#' found in the Spectronaut report. Defined as 'PEP.NrOfMissedCleavages' by
#' default.
#' @param isProteotypic A character string or numeric giving the column name or
#' column number in which annotation of the peptide is proteotypic can be found
#' in the Spectronaut report. Defined as 'PEP.IsProteotypic' by default.
#' @param isTryptic A character string or numeric giving the column name or
#' column number in which annotation of the digest type can be found in the
#' Spectronaut report. Defined as PEP.DigestType....Trypsin.P.' by default.
#' @param start A character string or numeric giving the column name or
#' column number in which start position of each peptide can be found in the
#' Spectronaut report. Defined as PEP.PeptidePosition' by default.
#'
#'  @return Returns a data.frame including all necessary annotation on peptides
#'  and proteins
#'
#' @export

getPepProtAnnot <- function(spectroOut,
                            spectroOut2=NULL,
                            analysisLvl="Peptide",
                            Precursor="EG.PrecursorId",
                            modPeptide="EG.ModifiedSequence",
                            Peptide="PEP.StrippedSequence",
                            Protein="PG.ProteinGroups",
                            AllProtein="PEP.AllOccurringProteinAccessions",
                            NMissedCleavages="PEP.NrOfMissedCleavages",
                            isProteotypic="PEP.IsProteotypic",
                            isTryptic="PEP.DigestType....Trypsin.P.",
                            start="PEP.PeptidePosition"){

    ## checking if value of AnalysisLvl is set correctly
    if(!tolower(analysisLvl) %in% c("peptide", "modifiedpeptide", "precursor")){
        stop(analysisLvl, " is an unknown input for 'analysisLvl'. Please choose
             'peptide' or 'modifiedpeptide' or 'precursor'.")
    }

    ## join Spectronaut reports if two are provided
    if(!is.null(spectroOut2)){
        spectroOut <-rbind(spectroOut, spectroOut2)
    }

    ## creating basic peptide & protein annotation file
    annotPP <- data.frame(Peptide=spectroOut[, Peptide],
                          Protein=spectroOut[, Protein],
                          AllProtein=spectroOut[, AllProtein],
                          NMissedCleavages=spectroOut[, NMissedCleavages],
                          isProteotypic=spectroOut[, isProteotypic],
                          isTryptic =spectroOut[, isTryptic],
                          startPosition=spectroOut[, start])

    ## adding modified peptides/precursors based on 'analysisLvl'
    if(tolower(analysisLvl) == "modifiedpeptide"){
        annotPP$modPeptide <- gsub("[ _]", "", spectroOut[, modPeptide])
        annotPP <- annotPP[, c(ncol(annotPP), 1:ncol(annotPP)-1)]
    }
    if(tolower(analysisLvl) == "precursor"){
        annotPP$Precursor <- gsub("[ _]", "", spectroOut[, Precursor])
        annotPP <- annotPP[, c(ncol(annotPP), 1:ncol(annotPP)-1)]
    }

    annotPP <- unique(annotPP) ## removing all duplicated rows from data.frame
    iRow <- switch(tolower(analysisLvl), "peptide"="Peptide",
                   "modifiedpeptide"="modPeptide", "precursor"="Precursor")

    ## join protein names if peptide was matched to different PG groups in
    ## different Spectronaut reports
    if(!is.null(spectroOut2)){
        annotPP <- joinPG(annotPP, iRow)
    }
    ## get end position of peptides
    annotPP$endPosition <- getEndPositionOfPep(annotPP)

    row.names(annotPP) <- annotPP[, iRow]
    return(annotPP)
}

#' @title Estimating end position of peptide sequences
#'
#' @description Getting end position of peptides sequences from the starting
#' positions provided. Function takes into account that peptides can map to
#' multiple regions in the proteome, providing multiple end positions if
#' multiple start positions are given.
#'
#' @usage getEndPositionOfPep(annotPP, startCol="startPosition",
#' pepCol="Peptide")
#'
#' @param annotPP A data.frame with peptide and protein annotation. Rows are
#' features. The Data.frame must include a column providing the start position
#' of each peptide in the proteins and one containing the peptide sequences.
#' Different start positions for the same peptide should be separated by ',' for
#' peptides mapping to the same protein multiple times and ';' for peptides
#' mapping to multiple proteins.
#' @param startCol A character string or numeric giving the column name or
#' column number defining column in which start positions of peptides are
#' provided. Defined as 'startPosition' by default.
#' @param pepCol A character string or numeric giving the column name or
#' column number defining column in which peptide sequences are provide. Peptide
#' sequences have to be provided without modifications or charges. Set to
#' 'Peptide' by default.
#'
#' @return Returns a character vector with the end positions of all peptides.
#' Different start positions for the same peptide should be separated by ',' for
#' peptides mapping to the same protein multiple times and ';' for peptides
#' mapping to multiple proteins.

getEndPositionOfPep <- function(annotPP, startCol="startPosition",
                                pepCol="Peptide"){
    end <- sapply(seq(1:nrow(annotPP)), function(i){
        start <- as.character(annotPP[, startCol][i])
        length <- nchar(annotPP[, pepCol][i])-1

        ## estimate end position if  one start position of peptide is provided
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

#' @title Combining protein names of identical peptides
#'
#' @description Joining protein names together if a peptide was matched to
#' different PG groups in different Spectronaut reports.
#'
#' @param annotPP A data.frame with peptide and protein annotation.
#' Rows are feature and must include a column protein name(s)
#' @param iCol A character string or numeric giving the column name or
#' column number in which features are provided in the Spectronaut report
#' based on the chosen \code{analysisLvl}.
#' @param protCol A character string or numeric giving the column name or
#' column number in which protein name(s)  are provide. Set to 'Protein' by
#' default.
#'
#' @return A data.frame including all necessary annotation on peptides and
#' proteins with joined protein names.

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

#' @title Creating data.frame with sample annotation from Spectronaut report
#'
#' @description Creates an annotation file on samples including information on
#' the conditions of samples. This function only runs, if all the annotations
#' were added to the Spectronaut report, else please create the sample
#' annotation file yourself.
#'
#' @usage getSampleAnnot(spectroOut, sampleName="R.FileName",
#' sampleCondition="R.Condition", typeCondition="factor",
#' contrastCoding="dummyCoding", baseLevel=NULL)
#'
#' @param spectroOut Spectronaut report of MS data. Spectronaut report should
#' have been exported using the Spectronaut schema
#' 'SpectroSchema_LiPAnalyzerOut'.
#' @param sampleName A character string or numeric giving the column name or
#' column number in which sample names are defined in the Spectronaut report.
#' Set to R.FileName' by default.
#' @param sampleCondition A character string or numeric giving the column name
#' or  column number in which conditions can be found in the Spectronaut report.
#' Set to 'R.Condition' by
#' default.
#' @param typeCondition A character string providing information if condition is
#' a factor or continuous variable If condition is continuous it has to be
#' provided as numeric in the Spectronaut report. Defined as ' 'factor' is
#' default, alternatively set to 'continuous'/
#' @param contrastCoding A character string providing information which type of
#' contrast coding to use if condition is a factor. Defined as 'dummyCoding' by
#' default. Alternatively, 'sumCoding' or 'wecCoding' can be chosen. Can be set
#' to NULL to not define contrast coding automatically in case you want to
#' defined it yourself.
#' @param baseLevel A character strong giving name of reference level if
#' \code{typeCondition} is a factor. If not provided the function will choose
#' the first condition in Spectronaut report as the reference level.
#'
#' @return Returns a data.frame containing sample and condition annotation.
#' Rows are samples.
#'
#' @export

getSampleAnnot <- function(spectroOut,
                           sampleName="R.FileName",
                           sampleCondition="R.Condition",
                           typeCondition="factor",
                           contrastCoding="dummyCoding",
                           baseLevel=NULL){
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
            stats::contrasts(Condition) <- stats::contr.treatment(nlevels(
                Condition), base=which(levels(Condition) == baseLevel))

        }
        else if(tolower(contrastCoding) == "sumcoding"){
            message("Using sum coding for condition.")
            stats::contrasts(Condition) <- stats::contr.sum(Condition)
        }
        else if(tolower(contrastCoding) == "weccoding"){
            message("Using wec coding for condition.")
            stats::contrasts(Condition) <- wec::contr.wec(Condition,
                                                         omitted=baseLevel)
        }
    }
    else if(tolower(typeCondition) == "continuous"){
        Condition <- as.numeric(spectroOut[, sampleCondition])
    }
    annotS <- data.frame(SampleName=spectroOut[, sampleName],
                        Condition= Condition)
    annotS <- unique(annotS)
    row.names(annotS) <- annotS$SampleName
    return(annotS)
}
