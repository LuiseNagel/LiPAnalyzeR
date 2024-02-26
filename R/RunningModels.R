##matrixStats
globalVariables(names=c("Condition", "XPep", "XProt",  "Y"))

#' @title Fitting model for accessibility changes in LiP peptides.
#'
#' @description Function to build linear regression models for fitting MS data
#' and retrieving structural variation between different conditions.
#'
#' @usage  analyzeLiPPepData(quantityList, annotS, infoCondition="Condition",
#' formulaRUV="Y~XPep+XProt", formulaContrast=NULL, lowRUV=c(-1e9, 0, 0),
#' upRUV=c(Inf, Inf, Inf), addRUVbounds=FALSE, LiPonly=FALSE, withHT=FALSE)
#'
#' @param quantityList A list of matrices, containing peptide/protein quantities.
#' Rows represent features and columns refer to the samples. Names of the list
#' items should be set to "LiPPep", "TrpPep" and "TrpProt" (if the LiP only
#' version is run, names should refer to "LiPPep" and "LiPProt").
#' @param annotS A data.frame containing sample annotation. Must contain all
#' columns required for the RUV and contrast models. Rows are samples and
#' must match to columns of the matrices in \code{quantityList}.
#' @param infoCondition A character string providing column name of
#' \code{annotS} in which condition is provided. Default is 'Condition'.
#' @param formulaRUV A character string or formula defining the RUV models.
#' Default is defined as 'Y~XPep+XProt'.
#' @param formulaContrast A character string or formula defining the contrast
#' models. If 'NULL' will be set to 'Y~\code{infoCondition}.
#' @param lowRUV A numeric vector defining lower boundaries of the coefficients
#' of the RUV models. Elements refer to definition of \code{formulaRUV}.
#' Default is defined as 'c(-Inf, 0, 0)'.
#'@param upRUV A numeric vector defining upper boundaries of the coefficients
#' of the RUV models. Elements refer to definition of \code{formulaRUV}.
#' Default is defined as 'c(Inf, Inf, Inf)'.
#' @param addRUVbounds A boolean value, if set to 'TRUE' as many bounds as
#' additionally needed in each RUV model are added to \code{lowRUV} and
#' \code{upRUV}. Added boundaries are automatically set to -Inf for
#' \code{lowRUV} and Inf for \code{upRUV}. Important to set to 'TRUE', if you
#' are for example also running batch correction in the RUV model.
#' @param LiPonly A boolean value to set to 'TRUE' if you are running the
#' LiPonly version of the package and not providing trypsin-only data. If set
#' to TRUE' RUV boundaries will be adjusted, \code{lowRUV} will be set to
#' c(-1e9, 0) and \code{RUV} will be set to c(Inf, Inf). Additionally,
#' \code{formulaRUV} will be adjusted. In case you want to add further
#' variables to the models, please use the \code{runModel} function. Default is
#' 'FALSE'.
#' @param withHT A boolean value to set to 'TRUE' if LiPPep should only be
#' corrected for TrpProt only. If set to TRUE' RUV boundaries will be adjusted,
#' \code{lowRUV} will be set to c(-1e9, 0) and \code{RUV} will be set to
#' c(Inf, Inf). Additionally, \code{formulaRUV} will be adjusted. In case you
#' want to add further variables to the models, please use the \code{runModel}
#' function. Default is 'FALSE'.
#'
#'
#' @export
analyzeLiPPepData <- function(quantityList, annotS, infoCondition="Condition",
                              formulaRUV="Y~XPep+XProt", formulaContrast=NULL,
                              lowRUV=c(-1e9, 0, 0), upRUV=c(Inf, Inf, Inf),
                              addRUVbounds=FALSE, LiPonly=FALSE, withHT=FALSE){
    if(is.null(formulaContrast)){
        formulaContrast <- paste0("Y~", infoCondition)
    }
    if(LiPonly|withHT){
        if(LiPonly){
            message("Running 'LiPonly' mode, performing RUV only with LiPProt.
If you want to add further variables to the RUV model please use the 'runModel'
function.")
        }
        else if(withHT){
            message("Running 'withHT' mode, performing RUV only with TrpProt.
If you want to add further variables to the RUV model please use the 'runModel'
function.")
        }
        formulaRUV <- "Y~XProt"
        lowRUV <- c(-1e9, 0)
        upRUV <- c(1e9, 1e9)
    }

    LiPOut <- runModel(quantityList=quantityList,
                       annotS=annotS,
                       formulaRUV=formulaRUV,
                       formulaContrast=formulaContrast,
                       lowRUV=lowRUV,
                       upRUV=upRUV,
                       addRUVbounds=addRUVbounds,
                       LiPonly=LiPonly,
                       withHT=withHT)
    return(LiPOut)
}

#' @title Fitting model for PK-independent changes in Trp peptides

#' @description Function to build linear regression models for fitting MS data
#' and retrieving PK-independent peptide variation between different conditions.
#'
#' @usage analyzeTrpPepData(quantityList, annotS, infoCondition="Condition",
#' formulaRUV="Y~XProt", formulaContrast=NULL, lowRUV=c(-1e9, 0),
#' upRUV=c(Inf, Inf), addRUVbounds=FALSE)
#'
#' @param quantityList A list of matrices, containing peptide/protein quantities.
#' Rows represent features and columns refer to the samples. Names of the list
#' items should be set to "LiPPep", "TrpPep" and "TrpProt".
#' @param annotS A data.frame containing sample annotation. Rows are samples and
#' must match to columns of \code{quantityList}.
#' @param infoCondition A character string providing column name of
#' \code{annotS} in which condition is provided. Default is 'Condition'.
#' @param formulaRUV A character string or formula defining the RUV models.
#' Default is defined as 'Y~XProt'.
#' @param formulaContrast A character string or formula defining the contrast
#' models. If 'NULL' will be set to 'Y~\code{infoCondition}.
#' @param lowRUV A numeric vector defining lower boundaries of the coefficients
#' of the RUV models. Elements refer to definition of \code{formulaRUV}.
#' Default is defined as 'c(-Inf, 0)'.
#'@param upRUV A numeric vector defining upper boundaries of the coefficients
#'of the RUV models. Elements refer to definition of \code{formulaRUV}.
#' Default is defined as 'c(Inf, Inf)'.
#' @param addRUVbounds A boolean value, if set to 'TRUE' as many bounds as
#' additionally needed in each RUV model are added to \code{lowRUV} and
#' \code{upRUV}. Added boundaries are automatically set to -Inf for
#' \code{lowRUV} and Inf for \code{upRUV}. Important to set to 'TRUE', if you
#' are for example also running batch correction in the RUV model.
#'
#' @export
analyzeTrpPepData <- function(quantityList, annotS, infoCondition="Condition",
                              formulaRUV="Y~XProt", formulaContrast=NULL,
                              lowRUV=c(-1e9, 0), upRUV=c(Inf, Inf),
                              addRUVbounds=FALSE){
    if(is.null(formulaContrast)){
        formulaContrast <- paste0("Y~", infoCondition)
    }
    quantityList <- list(Y=quantityList$TrpPep, XProt=quantityList$TrpProt)
    TrpOut <- runModel(quantityList=quantityList,
                       annotS=annotS,
                       formulaRUV=formulaRUV,
                       formulaContrast=formulaContrast,
                       lowRUV=lowRUV,
                       upRUV=upRUV,
                       addRUVbounds=addRUVbounds)
    colnames(TrpOut$modelCoeff)[2] <- c("TrpProt_RUV")
    return(TrpOut)
}

#' @title Fitting model for protein abundance changes

#' @description Function to build linear regression models for fitting MS data
#' and retrieving protein abundance variation between different conditions.
#'
#' @usage analyzeTrpProtData(quantityList, annotS, annotPP,
#' infoCondition="Condition", infoProtName="Protein", formulaRUV=NULL,
#' formulaContrast=NULL, lowRUV=NULL, upRUV=NULL, addRUVbounds=FALSE,
#' LiPonly=FALSE)
#'
#' @param quantityList A list of matrices, containing peptide/protein quantities.
#' Rows represent features and columns refer to the samples. Names of the list
#' items should be set to "LiPPep", "TrpPep" and "TrpProt" (if the LiP only
#' version is run, names should refer to "LiPPep" and "LiPProt").
#' @param annotS A data.frame containing sample annotation. Rows are samples and
#' must match to columns of \code{quantityList}.
#' @param annotPP A data.frame with peptide and protein annotatioon. Rows are
#' features and must match to the row.names of quantityList/QuantityMatrix.
#' @param infoCondition A character string providing column name of
#' \code{annotS} in which condition is provided. Default is 'Condition'.
#' @param infoProtName A character string providing column name of
#' \code{annotPP} in which the protein names are provided. Default is 'Protein'.
#' @param formulaRUV A character string or formula defining the RUV models.
#' Default is defined as 'NULL'.
#' @param formulaContrast A character string or formula defining the contrast
#' models. If 'NULL' will be set to 'Y~\code{infoCondition}
#' @param lowRUV A numeric vector defining lower boundaries of the coefficients
#' of the RUV models. Elements refer to definition of \code{formulaRUV}.
#' Default is defined as 'NULL'.
#' @param upRUV A numeric vector defining upper boundaries of the coefficients
#' of the RUV models. Elements refer to definition of \code{formulaRUV}.
#' Default is defined as 'NULL'.
#' @param addRUVbounds A boolean value, if set to 'TRUE' as many bounds as
#' additionally needed in each RUV model are added to \code{lowRUV} and
#' \code{upRUV}. Added boundaries are automatically set to -Inf for
#' \code{lowRUV} and Inf for \code{upRUV}. Important to set to 'TRUE', if you
#' are for example also running batch correction in the RUV model.
#' @param LiPonly A boolean value to set to 'TRUE' if you are running the
#' LiPonly version of the package and not providing trypsin-only data.

#'
#' @export
analyzeTrpProtData <- function(quantityList, annotS, annotPP,
                               infoCondition="Condition",
                               infoProtName="Protein", formulaRUV=NULL,
                               formulaContrast=NULL, lowRUV=NULL, upRUV=NULL,
                               addRUVbounds=FALSE, LiPonly=FALSE){
    if(is.null(formulaContrast)){
        formulaContrast <- paste0("Y~", infoCondition)
    }

    ## map peptides to proteins
    if(!LiPonly){
        protData <- quantityList$TrpProt
    }
    else{
        protData <- quantityList$LiPProt
    }
    Peps2Prots <- split(row.names(protData), annotPP[row.names(protData),
                                                     infoProtName])
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
                        addRUVbounds=addRUVbounds,
                        LiPonly=LiPonly)
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
#' addRUVbounds=FALSE, returnRUVmodels=FALSE, returnContrastmodels=FALSE,
#' LiPonly=FALSE, withHT=FALSE)
#'
#' @param quantityList A list of preprocessed matrices, containing quantities of
#' interest(e.g. peptide, modified peptide, precursor) and protein abundances.
#' Rows represent features and columns samples and should match between the
#' different matrices contained in the list.
#' Output from \code{preprocessQuantityMatrix} is in the correct format.
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
#' @param annotS A data.frame containing sample annotation. Must contain all
#' columns required for the RUV and contrast models. Rows are samples and
#' must match to columns of the matrices in \code{quantityList}.
#' @param formulaRUV A character string or formula defining the bounded
#' variable least square models. Use 'Y' for defining the quantity matrix with
#' values to predict. ´ If additional quantity matrices should be used as
#' variables in the models, please refer to these as 'XPep' and/or 'XProt'. All
#' other variables in the formula should refer to columns in \code{annotS}. If
#' RUV should not be run, set to 'NULL'. Default is defined as 'Y~XPep+XProt'.
#' @param formulaContrast A character string or formula defining the ordinary
#' variable least square models. Use 'Y' for defining the quantity matrix with
#' values to predict. ´If additional quantity matrices should be used as
#' variables in the models, please refer to these as 'XPep' and/or 'XProt'.
#' All other variables in the formula should refer to columns in \code{annotS}.
#' If contrast models should not be run, set to 'NULL'. Default is defined as
#' 'Y~Condition'. If set to 'NULL', the residuals of the RUV models will
#' additionally be returned.
#' @param lowRUV A numeric vector defining lower boundaries of the coefficients
#' of the RUV models. Elements refer to definition of \code{formulaRUV}.
#' Default is defined as 'c(-Inf, 0, 0)'
#' @param upRUV A numeric vector defining upper boundaries of the coefficients
#' of the RUV models. Elements refer to definition of \code{formulaRUV}.
#' Default is defined as 'c(Inf, Inf, Inf)'
#' @param addRUVbounds A boolean value, if set to 'TRUE' as many bounds as
#' additionally needed in each RUV model are added to \code{lowRUV} and
#' \code{upRUV}. Added boundaries are automatically set to -Inf for
#' \code{lowRUV} and Inf for \code{upRUV}. Important to set to 'TRUE', if you
#' are for example also running batch correction in the RUV model.
#' @param returnRUVmodels A boolean value, set to 'TRUE' if you want the
#' function to additionally return all RUV models.
#' @param returnContrastmodels A boolean value, set to 'TRUE' if you want the
#' function to additionally return all contrast models
#' @param LiPonly A boolean value to set to 'TRUE' if you are running the
#' LiPonly version of the package and not providing trypsin-only data. Default
#' is set to 'FALSE'.
#' @param withHT A boolean value to set to 'TRUE' if LiPPep should only be
#' corrected for TrpProt only. Default is set to 'FALSE'.
#'
#' @return If RUV & contrast model or only the contrast model is run, a list of
#' two data.frames will be returned. The first contains the model coefficients
#' of all models, the second one provides the p-values estimated in the contrast
#' model. If only RUV is run, it will return a list of two data.frame, the first
#' containing the residuals of the RUV model, the second one provides the
#' coefficients from the model. If \code{returnContrastmodels} and/or
#' \code{returnRUVmodels} is set to 'TRUE', least square models for each
#' peptide will be exported as additional list elements.
#'
#' @export
runModel <- function(quantityList, annotS=NULL, formulaRUV="Y~XPep+XProt",
                     formulaContrast="Y~Condition", lowRUV=c(-1e9, 0, 0),
                     upRUV=c(Inf, Inf, Inf), addRUVbounds=FALSE,
                     returnRUVmodels=FALSE, returnContrastmodels=FALSE,
                     LiPonly=FALSE, withHT=FALSE){

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
        stop("Please provide 'formula' in the correct class. Use 'character' or
'formula'.")
    }

    ## checking if 'formulaRUV' contains unexpected variables
    if(LiPonly){
        if(grepl("XPep", as.character(formulaRUV))){
            stop("Function is run in 'LiPonly' mode, but 'formulaRUV' includes
XPep. Please adjust 'formulaRUV'.")
        }
    }

    if(withHT){
        if(grepl("XPep", as.character(formulaRUV))){
            stop("Function is run in 'withHT' mode, but 'formulaRUV' includes
XPep. Please adjust 'formulaRUV'.")
        }
    }

    ## adjusting names of quantityList matrices
    if(paste(names(quantityList), collapse="") == c("LiPPepTrpPepTrpProt")){
        names(quantityList) <- c("Y", "XPep", "XProt")
    }

    else if(paste(names(quantityList), collapse="") == c("LiPPepLiPProt")|
            paste(names(quantityList), collapse="") == c("LiPPepTrpProt")){
        names(quantityList) <- c("Y", "XProt")
    }

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
    message("Using ", length(feat), " peptides/proteins and ", length(samples),
            " samples.")
    quantityList <- lapply(quantityList, function(x){
        x[feat, samples]
    })
    if(!is.null(annotS)){
        annotS <- annotS[samples, ]
    }

    ## stop if no formulas are provided
    if(is.null(formulaRUV) & is.null(formulaContrast)){
        stop("No RUV or contrast formula provided. Please provide at least one
of them to define the model(s) you want to run.")
    }
    resAll <- NULL

    ## running RUV models
    if(!is.null(formulaRUV)){
        message("Running bounded variable least square models.")
        modelMat <- createModelMatrix(quantityList, formulaRUV, annotS, samples)
        modelRUV <- runRUV(formulaRUV, modelMat, lowRUV, upRUV,
                             addRUVbounds)
        resRUV <- extractRUV(modelRUV, samples)
        if(is.null(formulaContrast)){
            resAll <- resRUV
            message("Returning RUV results including residuals and estimated
                    coefficients.")
        }
        dfRUV <- ncol(resRUV[[2]])-1
        quantityList[["Y"]] <- resRUV[[1]]
    }
    else{
        dfRUV <- NULL
    }

    ## running contrast models
    if(!is.null(formulaContrast)){
        message("Running ordinary least square models.")
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

    ## add message if there TrpPep or TrpProt coefficients are very high
    if(any(grepl("XPep|XProt", colnames(resAll[[1]])))){
        coeffPepProt <- unlist(c(resAll[[1]][, grepl("XPep",
                                                     colnames(resAll[[1]]))],
                                 resAll[[1]][, grepl("XProt",
                                                     colnames(resAll[[1]]))]))
        if(max(stats::na.omit(coeffPepProt))>5){
            message("At least one peptide/protein coefficient is higher than 5,
                     please check the results of the effected peptides for
                     plausibility.")
        }
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


## @title Creating model matrices to run RUV or contrast models
##
## @description Creates one model matrices per peptide/protein which are given
## into the RUV or contrast model function afterwards
##
## @return list with matrices for linear modelling
createModelMatrix <- function(quantityList, formula, annotS, samples){

    list2env(quantityList, envir=environment()) ## write matrices to environment
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

## @title Running RUV models
##
## @description Function to run RUV on a list of model matrices, one element in
## the list should represent one peptide or protein.
##
## @return Returns complete RUV model
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

## @title Extracting information from RUV models
##
## @description Function to extract the coefficients and residuals for each
## peptide/protein from the RUV models
##
## @return List with the first element being a data.frame with residuals from
## the RUV models and the second element being a data.frame with coefficients
## from the RUV models
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

## @title Running contrast models
##
## @description Function to run contrast models on a list of model matrices,
## one element in the list should represent one peptide or protein.
##
## @return Returns complete contrast model
runContrast <- function(modelMat){
    modelRes <- lapply(modelMat, function(data){
        Y <- data$Y
        X <- data$X[, -1, drop=FALSE]
        stats::lm(Y ~ X)
    })
    return(modelRes)
}

## @title Extracting information from Contrast models
##
## @description Function to extract the coefficients and p-values for each
## peptide/protein from the Contrast models. If RUV was run before the contrast
## model, the p-values are estimated taken the degrees of freedom already used
## by the RUV into account.
##
## @return List with the first element being a data.frame with coefficients from
## the contrast models and the second element being a data.frame with p- values from
## the contrast models
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
                consumed by RUV.")
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

## @title Estimate p-values taking degrees of freedom into account
##
## @description Function estimates p-values for all coefficients in the contrast
## model
## model, but takes the degrees of freedom used when running the RUV into
## account.
##
## @return Returns a data.frame with p-values
calcualtePvalAfterRUV <- function(LM,coeffTval, dfRUV){
    modelPv <- do.call(plyr::rbind.fill.matrix, lapply(LM, function(x){
        x <- stats::summary.lm(x)
        df <- x$df[2] - dfRUV

        if(df<1){
            warning("Not enough degrees of freedom to estimate p-values,
                    model is not reliable! Returning NAs for affected peptides/
                    proteins.")
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
