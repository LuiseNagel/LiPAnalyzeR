globalVariables(names=c("XPep", "XProt", "Y"))

#' @title Fitting model to quantify  accessibility changes in LiP quantities
#'
#' @description Function to build linear regression models for fitting MS data
#' and retrieving structural variation between different conditions.
#'
#' @usage  analyzeLiPPepData(quantityList, annotS, infoCondition="Condition",
#' formulaRUV="Y~XPep+XProt", formulaContrast, lowRUV=c(-1e9,0,0),
#' upRUV=c(Inf, Inf, Inf), addRUVbounds=FALSE, mode="default")
#'
#' @param quantityList A list of preprocessed matrices, containing quantities of
#' interest(e.g. peptide, modified peptide, precursor) and protein abundances.
#' Rows represent features and columns samples and should match between the
#' different matrices contained in the list.
#' Output from \code{preprocessQuantityMatrix} is in the correct format.
#' The matrices of \code{quantityList} may have the following names
#' \itemize{
#'   \item c('LiPPep', 'TrPPep', 'TrPProt') for
#'   \code{mode = c("default", "HTonly")}
#'   \item c('LiPPep', 'TrPProt') for \code{mode = "FTHTjoin"}
#'   \item c('LiPPep', 'LiPProt') for \code{mode = "LiPonly"}
#'   }
#' @param annotS A data.frame containing sample annotation. Rows are samples and
#' must match to columns of the matrices in \code{quantityList}. Must include
#' columns of any further variables used in \code{formulaRUV} and
#' \code{formulaContrast}, including \code{infoCondition}.
#' @param infoCondition A character string providing column name of
#' \code{annotS} were condition to fit in the contrast model is provided.
#' Default is 'Condition'.
#' @param formulaRUV A character string or formula defining the RUV models.
#' Default is 'Y~XPep+XProt'.
#' @param formulaContrast A character string or formula defining the contrast
#' models.
#' Default is 'NULL', causing the function to set \code{formulaContrast} to
#' 'Y~\code{infoCondition}.
#' @param lowRUV A numeric vector defining lower boundaries of the coefficients
#' of the RUV models. Elements refer to definition of \code{formulaRUV}.
#' Default is 'c(-Inf, 0, 0)'.
#' @param upRUV A numeric vector defining upper boundaries of the coefficients
#' of the RUV models. Elements refer to definition of \code{formulaRUV}.
#' Default is 'c(Inf, Inf, Inf)'.
#' @param addRUVbounds A boolean value, if set to 'TRUE' as many bounds as
#' additionally needed based on \code{formulaRUV} in each RUV model are added to
#' \code{lowRUV} and \code{upRUV}. Added boundaries are automatically set to
#' \code{lowRUV = -Inf} and \code{upRUV = Inf}. Important to set to 'TRUE', if
#' you have categories with multiple levels in the RUV model and did not adjust
#' the RUV boundaries based in the number of levels. This might be the case if
#' you have more than two batches you aim to account for in the RUV model.
#' Default is 'FALSE'.
#' @param mode A character variable defining mode in which function is run. Can
#' be set to c('default', 'HTonly', 'FTHTjoin' or 'LiPonly').
#' \itemize{
#'   \item 'default': Correcting full-tryptic LiP quantities (e.g. peptides) for
#'   TrP peptide and protein quantities.
#'   \item 'HTonly': Correcting half-tryptic LiP quantities (e.g. peptides) for
#'   best matching TrP peptide and the corresponding protein quantities. Please
#'   run \code{preprocessQuantityMatrix} with \code{mode = 'HTonly'} to
#'   preprocess the \code{quantityList} prior to running the models.
#'   \item 'FTHTjoin': Correcting full-tryptic and half-tryptic LiP quantities
#'   (e.g. peptides) for TrP protein quantities. Please run
#'   \code{preprocessQuantityMatrix} with \code{mode = 'FTHTjoin'} to
#'   preprocess the \code{quantityList} prior to running the models.
#'   \item 'LiPonly': Correcting LiP quantities (e.g. peptides) for LiP protein
#'   quantities.
#'   }
#' When aiming to run different corrections (e.g. correct LiP peptide for TrP
#' peptide quantities) please use \code{runModel} instead. This function allows
#' for more concrete settings.
#'
#' @return If running RUV and contrast model a list containing the following
#' matrices is returned:
#' \itemize{
#'   \item Matrix containing coefficients estimated in both models.
#'   \item P-values estimated for coefficients of contrast model.
#'   \item Matrix with the residuals resulting from the RUV model, rows are
#'   features and columns are samples.
#' }
#' If only a contrast model is run a list containing the following matrices is
#' returned:
#' \itemize{
#'   \item Matrix containing coefficients estimated in the contrast model.
#'   \item P-values estimated for coefficients of contrast model.
#' }
#'
#' @export

analyzeLiPPepData <- function(quantityList, annotS, infoCondition="Condition",
                              formulaRUV="Y~XPep+XProt", formulaContrast=NULL,
                              lowRUV=c(-1e9, 0, 0), upRUV=c(Inf, Inf, Inf),
                              addRUVbounds=FALSE, mode="default"){

    ## Setting (and checking) contrast if necessary
    if(is.null(formulaContrast)){
        message("'formulaContrast' is set to 'NULL', function is setting
formulaContrast to Y ~ 'infoCondition'. If you wish to not run a contrast model
but only a RUV model please use the function 'runModel'.")
        if(!infoCondition %in% colnames(annotS)){
            stop("'infoCondition' is not a column provided in 'annotS'.")
        }
        formulaContrast <- paste0("Y~", infoCondition)
    }

    ## Checking input based on mode setting
    if(tolower(mode) == "default"|tolower(mode) == "htonly"){
        if(paste(names(quantityList), collapse="") != c("LiPPepTrPPepTrPProt")){
            stop(paste0("Running mode ='", mode, "', 'forumlaRUV' does not
include 'XPep' AND 'XProt. Please adjust 'formulaRUV'." ))
        }
        if(!(grepl("XPep", as.character(formulaRUV))&
             grepl("XProt", as.character(formulaRUV)))){
            stop(paste0("Running mode ='", mode, "'. Names of 'quantityList' do
not meet expectation to be 'LiPPep', 'TrpPep 'TrPProt'."))
        }
    }

    if(tolower(mode) == "fthtjoin"){
        if(paste(names(quantityList), collapse="") != c("LiPPepTrPProt")){
            stop("Running mode = 'FTHTjoin'. Names of 'quantityList' do not meet
expectation to be 'LiPPep' and 'TrPProt'.")
        }
        if(grepl("XPep", as.character(formulaRUV))){
            stop("Function is run in 'FTHTjoin' mode, but 'formulaRUV' includes
XPep. Please adjust 'formulaRUV'.")
        }
    }

    if(tolower(mode) == "liponly"){
        if(paste(names(quantityList), collapse="") != c("LiPPepLiPProt")){
            stop("Running mode = 'LiPonly'. Names of 'quantityList' do not meet
expectation to be 'LiPPep' and 'LiPProt'.")
        }
        if(grepl("XPep", as.character(formulaRUV))){
            stop("Function is run in 'LiPonly' mode, but 'formulaRUV' includes
XPep. Please adjust 'formulaRUV'.")
        }
    }

        LiPOut <- runModel(quantityList=quantityList,
                           annotS=annotS,
                           formulaRUV=formulaRUV,
                           formulaContrast=formulaContrast,
                           lowRUV=lowRUV,
                           upRUV=upRUV,
                           addRUVbounds=addRUVbounds)
    return(LiPOut)
}

#' @title Fitting model to quantify PK-independent changes in TrP peptides

#' @description Function to build linear regression models for fitting MS data
#' and retrieving PK-independent peptide variations between different
#' conditions.
#'
#' @usage analyzeTrPPepData(quantityList, annotS, infoCondition="Condition",
#' formulaRUV="Y~XProt", formulaContrast=NULL, lowRUV=c(-1e9, 0),
#' upRUV=c(Inf, Inf), addRUVbounds=FALSE)
#'
#' @param quantityList A list of preprocessed matrices, containing quantities of
#' interest(e.g. peptide, modified peptide, precursor) and protein abundances.
#' Rows represent features and columns samples and should match between the
#' different matrices contained in the list.
#' Output from \code{preprocessQuantityMatrix} is in the correct format.
#' The matrices of \code{quantityList} may have the following names
#' \itemize{
#'   \item c('LiPPep', 'TrPPep', 'TrPProt') - in this case, 'LiPPep' will not
#'   be taken into account.
#'   \item c('TrPPep', 'TrPProt)
#'   }
#' @param annotS A data.frame containing sample annotation. Rows are samples and
#' must match to columns of the matrices in \code{quantityList}. Must include
#' columns of any further variables used in \code{formulaRUV} and
#' \code{formulaContrast}, including \code{infoCondition}.
#' @param infoCondition A character string providing column name of
#' \code{annotS} were condition to fit in the contrast model is provided.
#' Default is 'Condition'.
#' @param formulaRUV A character string or formula defining the RUV models.
#' Default is defined as 'Y~XProt'.
#' @param formulaContrast A character string or formula defining the contrast
#' models.
#' Default is 'NULL', causing the function to set \code{formulaContrast} to
#' 'Y~\code{infoCondition}.
#' @param lowRUV A numeric vector defining lower boundaries of the coefficients
#' of the RUV models. Elements refer to definition of \code{formulaRUV}.
#' Default is 'c(-Inf, 0)'.
#'@param upRUV A numeric vector defining upper boundaries of the coefficients
#'of the RUV models. Elements refer to definition of \code{formulaRUV}.
#' Default is 'c(Inf, Inf)'.
#' @param addRUVbounds A boolean value, if set to 'TRUE' as many bounds as
#' additionally needed based on \code{formulaRUV} in each RUV model are added to
#' \code{lowRUV} and \code{upRUV}. Added boundaries are automatically set to
#' \code{lowRUV = -Inf} and \code{upRUV = Inf}. Important to set to 'TRUE', if
#' you have categories with multiple levels in the RUV model and did not adjust
#' the RUV boundaries based in the number of levels. This might be the case if
#' you have more than two batches you aim to account for in the RUV model.
#' Default is 'FALSE'.
#'
#' @return If running RUV and contrast model a list containing the following
#' matrices is returned:
#' \itemize{
#'   \item Matrix containing coefficients estimated in both models.
#'   \item P-values estimated for coefficients of contrast model.
#'   \item Matrix with the residuals resulting from the RUV model, rows are
#'   features and columns are samples.
#' }
#' If only a contrast model is run a list containing the following matrices is
#' returned:
#' \itemize{
#'   \item Matrix containing coefficients estimated in the contrast model.
#'   \item P-values estimated for coefficients of contrast model.
#' }
#'
#' @export

analyzeTrPPepData <- function(quantityList, annotS, infoCondition="Condition",
                              formulaRUV="Y~XProt", formulaContrast=NULL,
                              lowRUV=c(-1e9, 0), upRUV=c(Inf, Inf),
                              addRUVbounds=FALSE){

    ## Setting (and checking) contrast if necessary
    if(is.null(formulaContrast)){
        message("'formulaContrast' is set to 'NULL', function is setting
formulaContrast to Y ~ 'infoCondition'. If you wish to not run a contrast model
but only a RUV model please use the function 'runModel'.")
        if(!infoCondition %in% colnames(annotS)){
            stop("'infoCondition' is not a column provided in 'annotS'.")
        }
        formulaContrast <- paste0("Y~", infoCondition)
    }

    ## Checking names of quantityList
    if(paste(names(quantityList), collapse="") != c("LiPPepTrPPepTrPProt") &
       paste(names(quantityList), collapse="") != c("TrPPepTrPProt")){
        stop("Aiming to analyze 'TrPPep' data. Names of 'quantityList' do not
meet expectation to be 'LiPPep', 'TrpPep', 'TrpProt' OR 'TrpPep', 'TrPProt'.")
    }

    ## Running models
    quantityList <- list(Y=quantityList$TrPPep, XProt=quantityList$TrPProt)
    TrPOut <- runModel(quantityList=quantityList,
                       annotS=annotS,
                       formulaRUV=formulaRUV,
                       formulaContrast=formulaContrast,
                       lowRUV=lowRUV,
                       upRUV=upRUV,
                       addRUVbounds=addRUVbounds)

    return(TrPOut)
}

#' @title Fitting model to quantify protein abundance changes

#' @description Function to build linear regression models for fitting MS data
#' and retrieving protein abundance variation between different conditions.
#'
#' @usage analyzeProtData(quantityList, annotS, annotPP,
#' infoCondition="Condition", nameProtQuant="Protein", formulaRUV=NULL,
#' formulaContrast=NULL, lowRUV=NULL, upRUV=NULL, addRUVbounds=FALSE)
#'
#' @param quantityList A list of matrices, containing peptide/protein
#' quantities. Rows represent features and columns refer to the samples. Names
#' of the list items should be set to "LiPPep", "TrPPep" and "TrPProt" (if the
#' LiP only version is run, names should refer to "LiPPep" and "LiPProt").
#' @param annotS A data.frame containing sample annotation. Rows are samples and
#' must match to columns of the matrices in \code{quantityList}. Must include
#' columns of any further variables used in \code{formulaRUV} and
#' \code{formulaContrast}, including \code{infoCondition}.
#' @param annotPP A data.frame with peptides (/modified peptides/precursors) and
#' protein annotation. Rows are features and the row names of the
#' \code{quantityList} matrices must be found here. \code{annotPP} must include
#' a column named \code{nameProtQuant} providing protein (group) names.
#' The output from \code{getPepProtAnnot} can be given here.
#' @param infoCondition A character string providing column name of
#' \code{annotS} were condition to fit in the contrast model is provided.
#' Default is 'Condition'.
#' @param nameProtQuant A character string giving column of \code{annotPP} were
#' protein names are provided.
#' Default is 'Protein'.
#' @param formulaRUV A character string or formula defining the RUV models.
#' Default is defined as 'NULL'.
#' @param formulaContrast A character string or formula defining the contrast
#' models.
#' Default is 'NULL', causing the function to set \code{formulaContrast} to
#' 'Y~\code{infoCondition}.
#' @param lowRUV A numeric vector defining lower boundaries of the coefficients
#' of the RUV models. Elements refer to definition of \code{formulaRUV}.
#' Default is 'NULL'.
#' @param upRUV A numeric vector defining upper boundaries of the coefficients
#' of the RUV models. Elements refer to definition of \code{formulaRUV}.
#' Default is 'NULL'.
#' @param addRUVbounds A boolean value, if set to 'TRUE' as many bounds as
#' additionally needed based on \code{formulaRUV} in each RUV model are added to
#' \code{lowRUV} and \code{upRUV}. Added boundaries are automatically set to
#' \code{lowRUV = -Inf} and \code{upRUV = Inf}. Important to set to 'TRUE', if
#' you have categories with multiple levels in the RUV model and did not adjust
#' the RUV boundaries based in the number of levels. This might be the case if
#' you have more than two batches you aim to account for in the RUV model.
#' Default is 'FALSE'.
#'
#' @return If running RUV and contrast model a list containing the following
#' matrices is returned:
#' \itemize{
#'   \item Matrix containing coefficients estimated in both models.
#'   \item P-values estimated for coefficients of contrast model.
#'   \item Matrix with the residuals resulting from the RUV model, rows are
#'   features and columns are samples.
#' }
#' If only a contrast model is run a list containing the following matrices is
#' returned:
#' \itemize{
#'   \item Matrix containing coefficients estimated in the contrast model.
#'   \item P-values estimated for coefficients of contrast model.
#' }
#'
#' @export

analyzeProtData <- function(quantityList, annotS, annotPP,
                            infoCondition="Condition",
                            nameProtQuant="Protein", formulaRUV=NULL,
                            formulaContrast=NULL, lowRUV=NULL, upRUV=NULL,
                            addRUVbounds=FALSE){

    ## Setting (and checking) contrast if necessary
    if(is.null(formulaContrast)){
        message("'formulaContrast' is set to 'NULL', function is setting
formulaContrast to Y ~ 'infoCondition'. If you wish to not run a contrast model
but only a RUV model please use the function 'runModel'.")
        if(!infoCondition %in% colnames(annotS)){
            stop("'infoCondition' is not a column provided in 'annotS'.")
        }
        formulaContrast <- paste0("Y~", infoCondition)
    }

    ## check quantityList input
    if(sum(names(quantityList) %in% c("LiPProt", "TrPProt"))>1){
        stop("Providing more than one protein matrix in 'quantityList', please
provide only one matrix (either 'TrPProt' or 'LiPProt'.")
    }

    if(sum("TrPProt" %in% names(quantityList)) == 1){
        message("Analyzing TrPProt data.")
        protData <- quantityList$TrPProt
    }

    else if(sum("LiPProt" %in% names(quantityList)) == 1){
        message("Analyzing LiPProt data.")
        protData <- quantityList$LiPProt
    }
    else{
        stop("No protein matrix provided in  'quantityList Please provide a
matrix in list names either 'TrPProt' or 'LiPProt'.")
    }

    ## map peptides to proteins
    Peps2Prots <- split(row.names(protData), annotPP[row.names(protData),
                                                     nameProtQuant])
    Peps2Prots <- unlist(lapply(Peps2Prots, \(x) (x[1])))
    quantityList <- list(Y=protData[Peps2Prots,])
    row.names(quantityList$Y) <- names(Peps2Prots)

    ## running models on protein data
    ProtOut <- runModel(quantityList=quantityList,
                        annotS=annotS,
                        formulaRUV=formulaRUV,
                        formulaContrast=formulaContrast,
                        lowRUV=lowRUV,
                        upRUV=upRUV,
                        addRUVbounds=addRUVbounds)
    return(ProtOut)
}


#' @title Fitting linear regression models to MS data

#' @description Function to build linear regression models to
#' \enumerate{
#'    \item remove unwanted variation from MS data (RUV step)
#'    \item estimate effect sizes and corresponding p-values for a condition of
#'    interest (contrast modeling step)
#'    }
#'
#' @usage runModel(quantityList, annotS=NULL, formulaRUV="Y~XPep+XProt",
#' formulaContrast="Y~Condition", lowRUV=c(-1e9, 0, 0), upRUV=c(Inf, Inf, Inf),
#' addRUVbounds=FALSE, returnRUVmodels=FALSE, returnContrastmodels=FALSE)
#'
#' @param quantityList A list of preprocessed matrices, containing quantities of
#' interest(e.g. peptide, modified peptide, precursor) and protein abundances.
#' Rows represent features and columns samples and should match between the
#' different matrices contained in the list.
#' Output from \code{preprocessQuantityMatrix} is in the correct format.
#' The matrices of \code{quantityList} may have the following names
#' \itemize{
#'   \item c('LiPPep', 'TrPPep', 'TrPProt')
#'   \item c('LiPPep', 'TrPPep')
#'   \item c('LiPPep', 'TrPProt')
#'   \item c('LiPPep', 'LiPProt')
#'   or use the variable naming 'Y', 'XPep' and 'XProt'.
#'   }
#' @param annotS A data.frame containing sample annotation. Rows are samples and
#' must match to columns of the matrices in \code{quantityList}. Must include
#' columns of any further variables used in \code{formulaRUV} and
#' \code{formulaContrast}.
#' @param formulaRUV A character string or formula defining the RUV models
#' performing bounded variable least square regression. Use 'Y' for defining the
#' quantity matrix with values to predict. If additional quantity matrices
#' should be used as variables in the models, please refer to these as 'XPep'
#' and/or 'XProt'. All other variables in the formula should refer to columns in
#' \code{annotS}. If RUV should not be run, set to 'NULL'.
#' Default is 'Y~XPep+XProt'.
#' @param formulaContrast A character string or formula defining the contrast
#' models performing ordinary variable least square models. Use 'Y' for defining
#' the quantity matrix with values to predict. ´If additional quantity matrices
#' should be used as variables in the models, please refer to these as 'XPep'
#' and/or 'XProt'. All other variables in the formula should refer to columns in
#' \code{annotS}.If contrast models should not be run, set to 'NULL'. If set to
#' 'NULL', the residuals of the RUV models will additionally be returned.
#' Default is Y~Condition'.
#' @param lowRUV A numeric vector defining lower boundaries of the coefficients
#' of the RUV models. Elements refer to definition of \code{formulaRUV}.
#' Default is defined as 'c(-Inf, 0, 0)'.
#' @param upRUV A numeric vector defining upper boundaries of the coefficients
#' of the RUV models. Elements refer to definition of \code{formulaRUV}.
#' Default is defined as 'c(Inf, Inf, Inf)'.
#' @param addRUVbounds A boolean value, if set to 'TRUE' as many bounds as
#' additionally needed based on \code{formulaRUV} in each RUV model are added to
#' \code{lowRUV} and \code{upRUV}. Added boundaries are automatically set to
#' \code{lowRUV = -Inf} and \code{upRUV = Inf}. Important to set to 'TRUE', if
#' you have categories with multiple levels in the RUV model and did not adjust
#' the RUV boundaries based in the number of levels. This might be the case if
#' you have more than two batches you aim to account for in the RUV model.
#' Default is 'FALSE'.
#' @param returnRUVmodels A boolean value, set to 'TRUE' the function will
#' additionally return all RUV models.
#' Default is 'FALSE'.
#' @param returnContrastmodels A boolean value, set to 'TRUE' the function will
#' additionally return all contrast models
#' Default is 'FALSE'.
#'
#' @return If running RUV and contrast model a list containing the following
#' matrices is returned:
#' \itemize{
#'   \item Matrix containing coefficients estimated in both models.
#'   \item P-values estimated for coefficients of contrast model.
#'   \item Matrix with the residuals resulting from the RUV model, rows are
#'   features and columns are samples.
#' }
#' If only a RUV model is run a list containing the following matrices is
#' returned:
#' \itemize{
#'   \item Matrix containing coefficients estimated in the RUV model.
#'   \item Matrix with the residuals resulting from the RUV model, rows are
#'   features and columns are samples.
#' }
#' If only a contrast model is run a list containing the following matrices is
#' returned:
#' \itemize{
#'   \item Matrix containing coefficients estimated in the contrast model.
#'   \item P-values estimated for coefficients of contrast model.
#' }
#' If \code{returnContrastmodels} and/or \code{returnRUVmodels} is set to
#' 'TRUE', the RUV and/or contrast models for each peptide will be exported as
#' additional list element(s).
#'
#' @export

runModel <- function(quantityList, annotS=NULL, formulaRUV="Y~XPep+XProt",
                     formulaContrast="Y~Condition", lowRUV=c(-1e9, 0, 0),
                     upRUV=c(Inf, Inf, Inf), addRUVbounds=FALSE,
                     returnRUVmodels=FALSE, returnContrastmodels=FALSE){

    ## if necessary transforming formulas into character
    ## remove spaces in formula
    if(!is.null(formulaRUV)){
        if(inherits(formulaRUV,"formula")){
            formulaRUV <- Reduce(paste, deparse(formulaRUV))
        }
        formulaRUV <- gsub(" ", "", formulaRUV)
    }
    if(!is.null(formulaContrast)){
        if(inherits(formulaContrast,"formula")){
            formulaContrast <- Reduce(paste, deparse(formulaContrast))
        }
        formulaContrast <- gsub(" ", "", formulaContrast)
    }

    ## check input format of formulas
    if(!(is.null(formulaRUV)|is.character(formulaRUV))&
       (is.null(formulaContrast)|is.character(formulaContrast))){
        stop("Please provide 'formulaRUV' and 'formulaConstrast' in the correct
class. Use 'character' or 'formula'.")
    }
    ## stop if no formulas are provided
    if(is.null(formulaRUV) & is.null(formulaContrast)){
        stop("No RUV or contrast formula provided. Please provide at least one
of them for running models.")
    }
    resAll <- NULL

    ## checking if quantiyList is in data.frame/matrix format
    if(!all(unlist(lapply(quantityList, \(x)
                          inherits(x, c("matrix","data.frame")))))){
        stop("Elements of 'quantityList' have to be data.frames or
matrices.")
    }
    quantityList <- lapply(quantityList, as.data.frame)


    ## checking and potentially changing names of quantityList matrices
    if(sum(!names(quantityList) %in% c("LiPPep", "TrPPep", "TrPProt", "LiPProt",
                                       "Y", "XPep", "XProt"))>0){
        stop("Names of matrices in 'quantityList' not permitted. Please change
accordingly.")
    }

    if(paste(names(quantityList), collapse="") == c("LiPPepTrPPepTrPProt")){
        names(quantityList) <- c("Y", "XPep", "XProt")
    }

    if(paste(names(quantityList), collapse="") == c("LiPPepTrPPep")){
        names(quantityList) <- c("Y", "XPep")
    }

    else if(paste(names(quantityList), collapse="") == c("LiPPepLiPProt")|
            paste(names(quantityList), collapse="") == c("LiPPepTrPProt")){
        names(quantityList) <- c("Y", "XProt")
    }

    ## checking that max 1 Y, XPep and XProt are provided in quantityList
    if(sum(duplicated(names(quantityList)))>0){
        if(sum(grepl("Y", names(quantityList)))>1){
            stop("Multiple 'LiPPep'/'Y' matrices provided in the quantityList.
Please only povide one.")
        }
        if(sum(grepl("XPep", names(quantityList)))>1){
            stop("Multiple 'TrPPep'/'XPep' matrices provided in the
quantityList. Please only povide one.")
        }
        if(sum(grepl("XPep", names(quantityList)))>1){
            stop("Multiple 'TrPProt'/'LiPProt'/'XProt' matrices provided in the
quantityList. Please only povide one.")
        }
    }

    ###quantityList <- lapply(quantityList, as.data.frame)

    ## assuring row.names and colnames fit over all input data
    feat <- Reduce(intersect, lapply(quantityList, row.names))
    samples <- Reduce(intersect, lapply(quantityList, colnames))

    if(!is.null(annotS)){
        samples <- intersect(row.names(annotS), samples)
    }
    if(length(feat) == 0|length(samples) == 0){
        stop(length(feat), " peptides/proteins in ", length(samples), " samples
detected.\nPlease check your input as well as column and row names of all
provided data." )
    }
    message(length(feat), " quantities and ", length(samples),
            " samples used in models.")
    quantityList <- lapply(quantityList, function(x){
        x[feat, samples]
    })
    if(!is.null(annotS)){
        annotS <- annotS[samples, ]
    }

    ## running RUV models
    if(!is.null(formulaRUV)){
        message("Running RUV models.")
        modelMat <- createModelMatrix(quantityList, formulaRUV, annotS, samples)
        modelRUV <- runRUV(formulaRUV, modelMat, lowRUV, upRUV,
                             addRUVbounds)
        resRUV <- extractRUV(modelRUV, samples)
        if(is.null(formulaContrast)){
            resAll <- resRUV
            message("Returning RUV results including residuals and estimated
coefficients.")
        }
        else{
            dfRUV <- ncol(resRUV[[2]])-1
            quantityList[["Y"]] <- resRUV[[1]]
        }
    }
    else{
        dfRUV <- NULL
    }

    ## running contrast models
    if(!is.null(formulaContrast)){
        message("Running contrast models.")
        modelMat <- createModelMatrix(quantityList, formulaContrast, annotS,
                                      samples)
        modelContrast <- runContrast(modelMat)
        resContrast <- extractContrast(modelContrast, formulaContrast, dfRUV)
        if(!is.null(formulaRUV)){
            resContrast$modelCoeff <- cbind(resRUV$modelCoeff,
                                            resContrast$modelCoeff)
        }
        resAll <- resContrast
    }

    ## adding feature names to output
    resAll <- lapply(resAll, function(x){
        rownames(x) <- feat
        return(as.data.frame(x))
    })

    ## add message if there TrPPep or TrPProt coefficients are very high
    if(any(grepl("XPep|XProt", colnames(resAll[[1]])))){
        coeffPepProt <- unlist(c(resAll[[1]][, grepl("XPep",
                                                     colnames(resAll[[1]]))],
                                 resAll[[1]][, grepl("XProt",
                                                     colnames(resAll[[1]]))]))
        if(max(stats::na.omit(coeffPepProt))>5){
            message("At least one peptide/protein coefficient estiamted is
higher than 5, please manually check the results of the respective peptides for
plausibility.")
        }
    }

    ## add residuals to output, if not already included
    if(!is.null(formulaRUV)&!(is.null(formulaContrast))){
        resAll[[length(resAll)+1]] <- quantityList[["Y"]]
        names(resAll)[length(resAll)] <- "modelResid"
    }
    ## adding models to results if set in function input
    if(returnRUVmodels){
        names(modelRUV) <- row.names(resAll$modelCoeff)
        resAll[[length(resAll)+1]] <- modelRUV
        names(resAll)[length(resAll)] <- "modelRUV"
    }

    if(returnContrastmodels){
        names(modelContrast) <- row.names(resAll$modelCoeff)
        resAll[[length(resAll)+1]] <- modelContrast
        names(resAll)[length(resAll)] <- "modelContrast"
    }

    return(resAll)
}


#' @title Creating model matrices to run RUV or contrast models
#'
#' @description Creates one model matrices per peptide/protein which are
#' used in the RUV or contrast model function afterwards
#'
#' @usage createModelMatrix(quantityList, formula, annotS, samples)
#'
#' @param quantityList A list of preprocessed matrices, containing quantities of
#' interest(e.g. peptide, modified peptide, precursor) and protein abundances.
#' Rows represent features and columns samples and should match between the
#' different matrices contained in the list.
#' Output from \code{preprocessQuantityMatrix} have to use the variable naming
#' 'Y', 'XPep' and 'XProt'.
#' @param formula A formula providing structure of model matrices created in
#' this function
#' @param annotS A data.frame containing sample annotation. Must contain all
#' columns included in the RUV and contrast models. Rows are samples and must
#' match to columns of the matrices in \code{quantityList}. Must include
#' columns of any further variables used in \code{formulaRUV}.
#' @param samples A character vector providing sample names for the model.
#'
#' @return A list with model matrices for running the RUV or contrast models.

createModelMatrix <- function(quantityList, formula, annotS, samples){

    list2env(quantityList, envir=environment()) ## write matrices to environment

    ## create model matrices for each quantity
    modelMat <- lapply(as.list(seq(1, nrow(Y))), function(i){
        Y <- as.numeric(Y[i, ])
        formula <- stats::as.formula(formula)
        formulaVars <- formula.tools::get.vars(formula)[
            formula.tools::get.vars(formula)!="Y"]
        e <- environment()
        sapply(formulaVars, function(x){
            if(x == "XPep"){
                assign(x, as.numeric(XPep[i,]), envir=e)
            }
            else if(x == "XProt"){
                assign(x, as.numeric(XProt[i,]), envir=e)
            }
            else{
                assign(x, annotS[,x], envir=e)
            }
            return(NULL)
        })

        X <- stats::model.matrix(formula)
        attr(X, "samples") <- samples[as.numeric(row.names(X))]
        return(list(Y=Y[as.numeric(row.names(X))], X=X))
    })

    return(modelMat)
}


#' @title Running RUV models for every quantity provided
#'
#' @description Function to run RUV on a list of model matrices, one element in
#' the list should represent one peptide or protein.
#'
#' @usage runRUV(formula, modelMat, lowRUV, upRUV, addRUVbounds)
#'
#' @param formula A formula used to create \code{modelMat} for RUV models.
#' @param modelMat A list of model matrices to perform RUV using bounded
#' variable least square regression.
#' @param lowRUV A numeric vector defining lower boundaries of the coefficients
#' of the RUV models. Elements refer to definition of \code{formulaRUV}.
#' Default is 'c(-Inf, 0, 0)'.
#' @param upRUV A numeric vector defining upper boundaries of the coefficients
#' of the RUV models. Elements refer to definition of \code{formulaRUV}.
#' Default is 'c(Inf, Inf, Inf)'.
#' @param addRUVbounds A boolean value, if set to 'TRUE' as many bounds as
#' additionally needed based on \code{formulaRUV} in each RUV model are added to
#' \code{lowRUV} and \code{upRUV}. Added boundaries are automatically set to
#' \code{lowRUV = -Inf} and \code{upRUV = Inf}. Important to set to 'TRUE', if
#' you have categories with multiple levels in the RUV model and did not adjust
#' the RUV boundaries based in the number of levels. This might be the case if
#' you have more than two batches you aim to account for in the RUV model.
#' Default is 'FALSE'.
#'
#' @return Returns a list of all RUV models

runRUV <- function(formula, modelMat, lowRUV, upRUV, addRUVbounds){

    ## change -Inf to high negative number instead, since bvls() does not take
    ## -Inf as an input
    lowRUV[lowRUV == -Inf] <- -1e9

    modelRUV <- lapply(modelMat, function(data){

        Y <- data$Y
        X <- data$X

        ## adding as many as necessary -Inf and Inf to the boundaries for bvls()
        ## if 'addRUVbounds' is set to TRUE
        if(addRUVbounds){
            nX <- ncol(X)
            lowRep <- nX-length(lowRUV)
            lowRUV <- c(lowRUV, rep(-1e9, lowRep))
            upRep <- nX-length(upRUV)
            upRUV <- c(upRUV, rep(Inf, upRep))
        }

        ## running RUV
        RUV <- bvls::bvls(A=X,
                           b=Y,
                           bl= lowRUV,
                           bu=upRUV)

        ## add variable and sample annotation to RUV output
        varsRUV <- paste0(dimnames(X)[[2]], "_RUV")
        varsRUV[1] <- "Intercept_RUV"
        attr(RUV, "variables") <- varsRUV
        attr(RUV, "samples") <- attributes(X)$samples
        return(RUV)
    })
    return(modelRUV)
}

#' @title Extracting information from RUV models
#'
#' @description Function to extract the coefficients and residuals from each
#' RUV model
#'
#' @usage extractRUV(mRUV, samples)
#'
#' @param mRUV A list of RUV models estimated using bounded variable least
#' square regression
#' @param samples A character vector providing sample names for the model.
#'
#' @return A list were the first element is a data.frame with residuals from
#; the RUV models and the second element is a data.frame with coefficients
#' from the RUV models.

extractRUV <- function(mRUV, samples){

    ## extract residuals and coefficients from RUV models
    modelResid <- do.call(plyr::rbind.fill, lapply(mRUV, function(x){
        y <- as.data.frame(t(x$residuals))
        colnames(y) <- attributes(x)$samples
        return(y)
    }))

    modelCoeff <- do.call(plyr::rbind.fill, lapply(mRUV, function(x){
        y <- as.data.frame(t(x$x))
        colnames(y) <- attributes(x)$variables
        return(y)
    }))

    ## adding NA columns in residuals if peptide in sample could not be modeled
    startDf <- as.data.frame(matrix(NA, nrow=1, ncol=length(samples)))
    colnames(startDf) <- samples
    modelResid <- plyr::rbind.fill(startDf, modelResid)[-1,]

    return(list(modelResid=modelResid, modelCoeff=modelCoeff))
}


#' @title Running contrast models for every quantity provided
#'
#' @description Function to run contrast on a list of model matrices, one
#' element in the list should represent one peptide or protein.
#'
#' @usage runContrast(modelMat)
#'
#' @param modelMat A list of model matrices to perform contrast modeling on
#' facilitating ordinary least square regression.
#'
#' @return Returns a list of all contrast models

runContrast <- function(modelMat){
    modelRes <- lapply(modelMat, function(data){
        Y <- data$Y
        X <- data$X[, -1, drop=FALSE]
        stats::lm(Y ~ X)
    })
    return(modelRes)
}


#' @title Extracting information from contrast models
#'
#' @description Function to extract the coefficients and residuals from each
#' RUV model. If RUV was run before the contrast model, the p-values estimation
#' taks into account the degrees of freedom already used by the RUV step.
#'
#' @usage extractContrast(mContrast, formulaContrast, dfRUV, coeffPval="Pr(>|t|)",
#' coeffTval="t value")
#'
#' @param mContrast A list of contrast models
#' @param formulaContrast A formula used to create \code{modelMat} for
#' contrast models.
#' @param dfRUV A numberic value providing the number of degrees of freedom used
#' in the RUV models.
#' @param coeffPval A character variable giving the column name were p-values
#' are provided if \code{summary.lm(mRUV)}.
#' Default is "Pr(>|t|)".
#' @param coeffTval A character variable giving the column name were t-values
#' are provided if \code{summary.lm(mRUV)}.
#' Default is "t value".
#'
#' @return A list were the first element is a data.frame with residuals from
#; the contrast models and the second element is a data.frame with coefficients
#' from the contrast models.

extractContrast <- function(mContrast, formulaContrast, dfRUV,
                            coeffPval="Pr(>|t|)", coeffTval="t value"){

    ## extract coefficients from contrast models
    modelCoef <- do.call(plyr::rbind.fill, lapply(mContrast, function(x){
        as.data.frame(t(stats::coef(x)))
    }))

    ## if RUV was not run before, extract p-values directly from contrast models
    if(is.null(dfRUV)){
        modelPv <- do.call(plyr::rbind.fill, lapply(mContrast, function(x){
            as.data.frame(t(summary(x)$coefficients[, coeffPval]))
        }))
    }
    ## if RUV was run before, calculate p-values taking degrees of freedom
    ## used in RUV models into account
    else{
        message("Estimating p-values while removing degrees of freedom
previously used in the RUV models.")
        modelPv <- calcualtePvalAfterRUV(mContrast, coeffTval, dfRUV)
    }

    ## adjusting column names of coefficient and p-value data.frame
    if(ncol(modelCoef) == 2){
        formulaContrast <- stats::as.formula(formulaContrast)
        formulaVars <- formula.tools::get.vars(formulaContrast)[
            formula.tools::get.vars(formulaContrast)!="Y"]
        cols <- paste0(c("Intercept", formulaVars), "_Contrast")
    }
    else{
        cols <- colnames(modelCoef)[-1]
        cols <- substring(cols, 2)
        cols <- paste0(c("Intercept", cols), "_Contrast")
    }
    colnames(modelCoef) <- cols
    colnames(modelPv) <- cols

    return(list(modelCoeff=modelCoef, modelPv=modelPv))
}


#' @title Estimate p-values from contrast model taking degrees of freedom used
#' in RUV step into account
#'
#' @usage calcualtePvalAfterRUV(LM, coeffTval, dfRUV)
#'
#' @param LM A single linear regression model originating from contrast step
#' @param coeffTval A character variable giving the column name were t-values
#' are provided if \code{summary.lm(mRUV)}.
#' @param dfRUV A numberic value providing the number of degrees of freedom used
#' in the RUV models.
#'
#' @return Returns a data.frame with p-values.

calcualtePvalAfterRUV <- function(LM, coeffTval, dfRUV){
    modelPv <- do.call(plyr::rbind.fill.matrix, lapply(LM, function(x){
        x <- stats::summary.lm(x)
        df <- x$df[2] - dfRUV

        if(df<1){
            warning("Not enough degrees of freedom to estimate p-values, model
is not reliable! Returning NAs for affected quantity.")
            pv <- stats::setNames(rep(NA, nrow(x$coefficients)),
                                  row.names(x$coefficients))
            return(t(pv))
        }

        ## Estimating p-values from t-values of contrast models taking degrees
        ## of freedom used in RUV into account
        else{
            pv <- sapply(x$coefficients[, coeffTval], function(y){
                2*stats::pt(abs(y), df, lower.tail=FALSE)
            })
            return(t(pv))
        }
    }))

    return(modelPv)
}
