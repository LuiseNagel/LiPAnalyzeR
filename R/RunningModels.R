##matrixStats
globalVariables(names=c("Condition", "XPep", "XProt",  "Y"))

#' @title Fitting model for accessibility changes in LiP peptides.
#'
#' @description Function to build linear regression models for fitting MS data
#' and retrieving structural variation between different conditions.
#'
#' @usage  analyzeLiPPepData(spectroList, annotS, infoCondition="Condition",
#' formulaBVLS="Y~XPep+XProt", formulaOLS=NULL, lowBVLS=c(-1e9, 0, 0),
#' upBVLS=c(Inf, Inf, Inf), addBVLSbounds=FALSE, LiPonly=FALSE, withHT=FALSE)
#'
#' @param spectroList A list of matrices, containing peptide/protein quantities.
#' Rows represent features and columns refer to the samples. Names of the list
#' items should be set to "LiPPep", "TrpPep" and "TrpProt" (if the LiP only
#' version is run, names should refer to "LiPPep" and "LiPProt").
#' @param annotS A data.frame containing sample annotation. Rows are samples and
#' must match to columns of \code{spectroList}.
#' @param infoCondition A character string providing column name of
#' \code{annotS} in which condition is provided. Default is 'Condition'.
#' @param formulaBVLS A character string or formula defining the BVLS models.
#' Default is defined as 'Y~XPep+XProt'.
#' @param formulaOLS A character string or formula defining the OLS models.
#' If 'NULL' will be set to 'Y~\code{infoCondition}.
#' @param lowBVLS A numeric vector defining lower boundaries of the coefficients
#' of the BVLS models. Elements refer to definition of \code{formulaBVLS}.
#' Default is defined as 'c(-Inf, 0, 0)'.
#'@param upBVLS A numeric vector defining upper boundaries of the coefficients
#' of the BVLS models. Elements refer to definition of \code{formulaBVLS}.
#' Default is defined as 'c(Inf, Inf, Inf)'.
#' @param addBVLSbounds A boolean value, if set to 'TRUE' as many bounds as
#' additionally needed in each BVLS model are added to \code{lowBVLS} and
#' \code{upBVLS}. Added boundaries are automatically set to -Inf for
#' \code{lowBVLS} and Inf for \code{upBVLS}. Important to set to 'TRUE', if you
#' are for example also running batch correction in the BVLS model.
#' @param LiPonly A boolean value to set to 'TRUE' if you are running the
#' LiPonly version of the package and not providing trypsin-only data. If set
#' to TRUE' BVLS boundaries will be adjusted, \code{lowBVLS} will be set to
#' c(-1e9, 0) and \code{BVLS} will be set to c(Inf, Inf). Additionally,
#' \code{formulaBVLS} will be adjusted. In case you want to add further
#' variables to the models, please use the \code{runModel} function. Default is
#' 'FALSE'.
#' @param withHT A boolean value to set to 'TRUE' if LiPPep should only be
#' corrected for TrpProt only. If set to TRUE' BVLS boundaries will be adjusted,
#' \code{lowBVLS} will be set to c(-1e9, 0) and \code{BVLS} will be set to
#' c(Inf, Inf). Additionally, \code{formulaBVLS} will be adjusted. In case you
#' want to add further variables to the models, please use the \code{runModel}
#' function. Default is 'FALSE'.
#'
#'
#' @export
analyzeLiPPepData <- function(spectroList, annotS, infoCondition="Condition",
                              formulaBVLS="Y~XPep+XProt", formulaOLS=NULL,
                              lowBVLS=c(-1e9, 0, 0), upBVLS=c(Inf, Inf, Inf),
                              addBVLSbounds=FALSE, LiPonly=FALSE, withHT=FALSE){
    if(is.null(formulaOLS)){
        formulaOLS <- paste0("Y~", infoCondition)
    }
    if(LiPonly|withHT){
        if(LiPonly){
            message("Running 'LiPonly' mode and only regressing out LiPProt
            quantities from LiPPeps in BVLS model. If you want to add further
            variables to the BVLS please use the 'runModel' function.")
        }
        else if(withHT){
            message("Running 'withHT' mode and only regressing out LiPProt
            quantities from LiPPeps in BVLS model. If you want to add further
            variables to the BVLS please use the 'runModel' function.")
        }

        formulaBVLS <- "Y~XProt"
        lowBVLS <- c(-1e9, 0)
        upBVLS <- c(1e9, 1e9)
    }
    LiPOut <- runModel(spectroList=spectroList,
                       annotS=annotS,
                       formulaBVLS=formulaBVLS,
                       formulaOLS=formulaOLS,
                       lowBVLS=lowBVLS,
                       upBVLS=upBVLS,
                       addBVLSbounds=addBVLSbounds,
                       LiPonly=LiPonly,
                       withHT=withHT)
    return(LiPOut)
}

#' @title Fitting model for PK-independent changes in Trp peptides

#' @description Function to build linear regression models for fitting MS data
#' and retrieving PK-independent peptide variation between different conditions.
#'
#' @usage analyzeTrpPepData(spectroList, annotS, infoCondition="Condition",
#' formulaBVLS="Y~XProt", formulaOLS=NULL, lowBVLS=c(-1e9, 0),
#' upBVLS=c(Inf, Inf), addBVLSbounds=FALSE)
#'
#' @param spectroList A list of matrices, containing peptide/protein quantities.
#' Rows represent features and columns refer to the samples. Names of the list
#' items should be set to "LiPPep", "TrpPep" and "TrpProt".
#' @param annotS A data.frame containing sample annotation. Rows are samples and
#' must match to columns of \code{spectroList}.
#' @param infoCondition A character string providing column name of
#' \code{annotS} in which condition is provided. Default is 'Condition'.
#' @param formulaBVLS A character string or formula defining the BVLS models.
#' Default is defined as 'Y~XProt'.
#' @param formulaOLS A character string or formula defining the OLS models.
#' If 'NULL' will be set to 'Y~\code{infoCondition}.
#' @param lowBVLS A numeric vector defining lower boundaries of the coefficients
#' of the BVLS models. Elements refer to definition of \code{formulaBVLS}.
#' Default is defined as 'c(-Inf, 0)'.
#'@param upBVLS A numeric vector defining upper boundaries of the coefficients
#'of the BVLS models. Elements refer to definition of \code{formulaBVLS}.
#' Default is defined as 'c(Inf, Inf)'.
#' @param addBVLSbounds A boolean value, if set to 'TRUE' as many bounds as
#' additionally needed in each BVLS model are added to \code{lowBVLS} and
#' \code{upBVLS}. Added boundaries are automatically set to -Inf for
#' \code{lowBVLS} and Inf for \code{upBVLS}. Important to set to 'TRUE', if you
#' are for example also running batch correction in the BVLS model.
#'
#' @export
analyzeTrpPepData <- function(spectroList, annotS, infoCondition="Condition",
                              formulaBVLS="Y~XProt", formulaOLS=NULL,
                              lowBVLS=c(-1e9, 0), upBVLS=c(Inf, Inf),
                              addBVLSbounds=FALSE){
    if(is.null(formulaOLS)){
        formulaOLS <- paste0("Y~", infoCondition)
    }
    spectroList <- list(Y=spectroList$TrpPep, XProt=spectroList$TrpProt)
    TrpOut <- runModel(spectroList=spectroList,
                       annotS=annotS,
                       formulaBVLS=formulaBVLS,
                       formulaOLS=formulaOLS,
                       lowBVLS=lowBVLS,
                       upBVLS=upBVLS,
                       addBVLSbounds=addBVLSbounds)
    colnames(TrpOut$modelCoeff)[2] <- c("TrpProt_BVLS")
    return(TrpOut)
}

#' @title Fitting model for protein abundance changes

#' @description Function to build linear regression models for fitting MS data
#' and retrieving protein abundance variation between different conditions.
#'
#' @usage analyzeTrpProtData(spectroList, annotS, annotPP,
#' infoCondition="Condition", infoProtName="Protein", formulaBVLS=NULL,
#' formulaOLS=NULL, lowBVLS=NULL, upBVLS=NULL, addBVLSbounds=FALSE,
#' LiPonly=FALSE)
#'
#' @param spectroList A list of matrices, containing peptide/protein quantities.
#' Rows represent features and columns refer to the samples. Names of the list
#' items should be set to "LiPPep", "TrpPep" and "TrpProt" (if the LiP only
#' version is run, names should refer to "LiPPep" and "LiPProt").
#' @param annotS A data.frame containing sample annotation. Rows are samples and
#' must match to columns of \code{spectroList}.
#' @param annotPP A data.frame with peptide and protein annotatioon. Rows are
#' features and must match to the row.names of SpectroList/QuantityMatrix.
#' @param infoCondition A character string providing column name of
#' \code{annotS} in which condition is provided. Default is 'Condition'.
#' @param infoProtName A character string providing column name of
#' \code{annotPP} in which the protein names are provided. Default is 'Protein'.
#' @param formulaBVLS A character string or formula defining the BVLS models.
#' Default is defined as 'NULL'.
#' @param formulaOLS A character string or formula defining the OLS models.
#' If 'NULL' will be set to 'Y~\code{infoCondition}
#' @param lowBVLS A numeric vector defining lower boundaries of the coefficients
#' of the BVLS models. Elements refer to definition of \code{formulaBVLS}.
#' Default is defined as 'NULL'.
#' @param upBVLS A numeric vector defining upper boundaries of the coefficients
#' of the BVLS models. Elements refer to definition of \code{formulaBVLS}.
#' Default is defined as 'NULL'.
#' @param addBVLSbounds A boolean value, if set to 'TRUE' as many bounds as
#' additionally needed in each BVLS model are added to \code{lowBVLS} and
#' \code{upBVLS}. Added boundaries are automatically set to -Inf for
#' \code{lowBVLS} and Inf for \code{upBVLS}. Important to set to 'TRUE', if you
#' are for example also running batch correction in the BVLS model.
#' @param LiPonly A boolean value to set to 'TRUE' if you are running the
#' LiPonly version of the package and not providing trypsin-only data.

#'
#' @export
analyzeTrpProtData <- function(spectroList, annotS, annotPP,
                               infoCondition="Condition",
                               infoProtName="Protein", formulaBVLS=NULL,
                               formulaOLS=NULL, lowBVLS=NULL, upBVLS=NULL,
                               addBVLSbounds=FALSE, LiPonly=FALSE){
    if(is.null(formulaOLS)){
        formulaOLS <- paste0("Y~", infoCondition)
    }

    ## map peptides to proteins
    if(!LiPonly){
        protData <- spectroList$TrpProt
    }
    else{
        protData <- spectroList$LiPProt
    }
    Peps2Prots <- split(row.names(protData), annotPP[row.names(protData),
                                                     infoProtName])
    Peps2Prots <- unlist(lapply(Peps2Prots, \(x) (x[1])))
    spectroList <- list(Y=protData[Peps2Prots,])
    row.names(spectroList$Y) <- names(Peps2Prots)

    ## running models on protein data
    ProtOut <- runModel(spectroList=spectroList,
                        annotS=annotS,
                        formulaBVLS=formulaBVLS,
                        formulaOLS=formulaOLS,
                        lowBVLS=lowBVLS,
                        upBVLS=upBVLS,
                        addBVLSbounds=addBVLSbounds,
                        LiPonly=LiPonly)
    return(ProtOut)
}


#' @title Fitting models to MS data

#' @description Function to build linear regression models for fitting MS data.
#'
#' @usage runModel(spectroList, annotS=NULL, formulaBVLS="Y~XPep+XProt",
#' formulaOLS="Y~Condition", lowBVLS=c(-1e9, 0, 0), upBVLS=c(Inf, Inf, Inf),
#' addBVLSbounds=FALSE, returnBVLSmodels=FALSE, returnOLSmodels=FALSE,
#' LiPonly=FALSE, withHT=FALSE)
#'
#' @param spectroList A list of matrices, containing peptide/protein quantities.
#' Rows represent features and columns refer to the samples. Names of the list
#' items should be set to "LiPPep", "TrpPep" and "TrpProt" (if the LiP only
#' version is run, names should refer to "LiPPep" and "LiPProt").
#' @param annotS A data.frame containing sample annotation. Rows are samples and
#' must match to columns of \code{spectroList}.
#' @param formulaBVLS A character string or formula defining the bounded
#' variable least square models. Use 'Y' for defining the quantity matrix with
#' values to predict. ´ If additional quantity matrices should be used as
#' variables in the models, please refer to these as 'XPep' and/or 'XProt'. All
#' other variables in the formula should refer to columns in \code{annotS}. If
#' BVLS should not be run, set to 'NULL'. Default is defined as 'Y~XPep+XProt'.
#' @param formulaOLS A character string or formula defining the ordinary
#' variable least square models. Use 'Y' for defining the quantity matrix with
#' values to predict. ´If additional quantity matrices should be used as
#' variables in the models, please refer to these as 'XPep' and/or 'XProt'.
#' All other variables in the formula should refer to columns in \code{annotS}.
#' If OLS should not be run, set to 'NULL'. Default is defined as 'Y~Condition'.
#' If set to 'NULL', the residuals of the BVLS models will additionally be
#' returned.
#' @param lowBVLS A numeric vector defining lower boundaries of the coefficients
#' of the BVLS models. Elements refer to definition of \code{formulaBVLS}.
#' Default is defined as 'c(-Inf, 0, 0)'
#' @param upBVLS A numeric vector defining upper boundaries of the coefficients
#' of the BVLS models. Elements refer to definition of \code{formulaBVLS}.
#' Default is defined as 'c(Inf, Inf, Inf)'
#' @param addBVLSbounds A boolean value, if set to 'TRUE' as many bounds as
#' additionally needed in each BVLS model are added to \code{lowBVLS} and
#' \code{upBVLS}. Added boundaries are automatically set to -Inf for
#' \code{lowBVLS} and Inf for \code{upBVLS}. Important to set to 'TRUE', if you
#' are for example also running batch correction in the BVLS model.
#' @param returnBVLSmodels A boolean value, set to 'TRUE' if you want the
#' function to additionally return all BVLS models.
#' @param returnOLSmodels A boolean value, set to 'TRUE' if you want the
#' function to additionally return all OLS models
#' @param LiPonly A boolean value to set to 'TRUE' if you are running the
#' LiPonly version of the package and not providing trypsin-only data. Default
#' is set to 'FALSE'.
#' @param withHT A boolean value to set to 'TRUE' if LiPPep should only be
#' corrected for TrpProt only. Default is set to 'FALSE'.
#'
#' @return If run with BVLS & OLS or only OLS model it will return a list of two
#' data.frames, the first contains the model coefficients of the BVLS & OLS
#' model, the second one provides the p-values estimated in the OLS model. If
#' only BVLS is run, it will return a list of two data.frame, the first
#' containing the residuals of the BVLS model, the second one provides the
#' coefficients from the model. If \code{returnOLSmodels} and/or
#' \code{returnBVLSmodels} is set to 'TRUE', least square models for each
#' peptide will be exported as additional list elements.
#'
#' @export
runModel <- function(spectroList, annotS=NULL, formulaBVLS="Y~XPep+XProt",
                     formulaOLS="Y~Condition", lowBVLS=c(-1e9, 0, 0),
                     upBVLS=c(Inf, Inf, Inf), addBVLSbounds=FALSE,
                     returnBVLSmodels=FALSE, returnOLSmodels=FALSE,
                     LiPonly=FALSE, withHT=FALSE){

    ## if necessary transforming formulas into character
    ## remove spaces in formula
    if(!is.null(formulaBVLS)){
        if(inherits(formulaBVLS,"formula")){
            formulaBVLS <- Reduce(paste, deparse(formulaBVLS))
        }
        formulaBVLS <- gsub(" ", "", formulaBVLS)
    }
    if(!is.null(formulaOLS)){
        if(inherits(formulaOLS,"formula")){
            formulaOLS <- Reduce(paste, deparse(formulaOLS))
        }
        formulaOLS <- gsub(" ", "", formulaOLS)
    }

    ## check input format of formulas
    if(!(is.null(formulaBVLS)|is.character(formulaBVLS))&
       (is.null(formulaOLS)|is.character(formulaOLS))){
        stop("Please provide 'formula' in the correct class. Use 'character' or
             'formula'.")
    }

    ## checking if 'formulaBVLS' contains unexpected variables
    if(LiPonly){
        if(grepl("XPep", as.character(formulaBVLS))){
            stop("Function is run in 'LiPonly' mode, but 'formulaBVLS' includes
                 XPep. Please adjust 'formulaBVLS'.")
        }
    }

    if(withHT){
        if(grepl("XPep", as.character(formulaBVLS))){
            stop("Function is run in 'withHT' mode, but 'formulaBVLS' includes
                 XPep. Please adjust 'formulaBVLS'.")
        }
    }

    ## adjusting names of SpectroList matrices
    if(paste(names(spectroList), collapse="") == c("LiPPepTrpPepTrpProt")){
        names(spectroList) <- c("Y", "XPep", "XProt")
    }

    else if(paste(names(spectroList), collapse="") == c("LiPPepLiPProt")|
            paste(names(spectroList), collapse="") == c("LiPPepTrpProt")){
        names(spectroList) <- c("Y", "XProt")
    }

    ## assuring row.names and colnames fit over all input data
    feat <- Reduce(intersect, lapply(spectroList, row.names))
    samples <- Reduce(intersect, lapply(spectroList, colnames))

    if(!is.null(annotS)){
        samples <- intersect(row.names(annotS), samples)
    }
    if(length(feat) == 0|length(samples) == 0){
        stop(length(feat), " peptides/proteins in ", length(samples), " samples
             detected.\nPlease check your input as well as column and row names
             of all provided data." )
    }
    message("Using ", length(feat), " peptides/proteins and ", length(samples),
            " samples.")
    spectroList <- lapply(spectroList, function(x){
        x[feat, samples]
    })
    if(!is.null(annotS)){
        annotS <- annotS[samples, ]
    }

    ## stop if no formulas are provided
    if(is.null(formulaBVLS) & is.null(formulaOLS)){
        stop("No BVLS or OLS formula provided. Please provide at least one of
             them to define the model(s) you want to run.")
    }
    resAll <- NULL

    ## running BVLS models
    if(!is.null(formulaBVLS)){
        message("Running bounded variable least square models.")
        modelMat <- createModelMatrix(spectroList, formulaBVLS, annotS, samples)
        modelBVLS <- runBVLS(formulaBVLS, modelMat, lowBVLS, upBVLS,
                             addBVLSbounds)
        resBVLS <- extractBVLS(modelBVLS, samples)
        if(is.null(formulaOLS)){
            resAll <- resBVLS
            message("Returning BVLS results including residuals and estimated
                    coefficients.")
        }
        dfBVLS <- ncol(resBVLS[[2]])-1
        spectroList[["Y"]] <- resBVLS[[1]]
    }
    else{
        dfBVLS <- NULL
    }

    ## running OLS models
    if(!is.null(formulaOLS)){
        message("Running ordinary least square models.")
        modelMat <- createModelMatrix(spectroList, formulaOLS, annotS, samples)
        modelOLS <- runOLS(modelMat)
        resOLS <- extractOLS(modelOLS, formulaOLS, dfBVLS)
        if(!is.null(formulaBVLS)){
            resOLS$modelCoeff <- cbind(resBVLS$modelCoeff,resOLS$modelCoeff)
        }
        resAll <- resOLS
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
    if(returnBVLSmodels){
        resAll[[length(resAll)+1]] <- list(modelBVLS)
        names(resAll)[length(resAll)] <- "modelBVLS"
    }

    if(returnOLSmodels){
        resAll[[length(resAll)+1]] <- list(modelOLS)
        names(resAll)[length(resAll)] <- "modelOLS"
    }

    return(resAll)
}


## @title Creating model matrices to run BVLS or OLS on
##
## @description Creates one model matrices per peptide/protein which are given
## into the BVLS or OLS functions afterwards
##
## @return list with matrices for linear modelling
createModelMatrix <- function(spectroList, formula, annotS, samples){

    list2env(spectroList, envir=environment()) ## write matrices to environment
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

## @title Running BVLS models
##
## @description Function to run BVLS on a list of model matrices, one element in
## the list should represent one peptide or protein.
##
## @return Returns complete BVLS model
runBVLS <- function(formula, modelMat, lowBVLS, upBVLS, addBVLSbounds){

    ## change -Inf to high negative number instead, since bvls() does not take
    ## -Inf as an input
    lowBVLS[lowBVLS == -Inf] <- -1e9

    modelRes <- lapply(modelMat, function(data){

        Y <- data$Y
        X <- data$X

        ## adding as many as necessary -Inf and Inf to the boundaries for bvls()
        ## if 'addBVLSbounds' is set to TRUE
        if(addBVLSbounds){
            nX <- ncol(X)
            lowRep <- nX-length(lowBVLS)
            lowBVLS <- c(lowBVLS, rep(-1e9, lowRep))
            upRep <- nX-length(upBVLS)
            upBVLS <- c(upBVLS, rep(Inf, upRep))
        }

        ## running BVLS
        BVLS <- bvls::bvls(A=X,
                           b=Y,
                           bl= lowBVLS,
                           bu=upBVLS)

        ## add variable and sample annotation to BVLS output
        varsBVLS <- paste0(dimnames(X)[[2]], "_BVLS")
        varsBVLS[1] <- "Intercept_BVLS"
        attr(BVLS, "variables") <- varsBVLS
        attr(BVLS, "samples") <- attributes(X)$samples
        return(BVLS)
    })
    return(modelRes)
}

## @title Extracting information from BVLS models
##
## @description Function to extract the coefficients and residuals for each
## peptide/protein from the BVLS models
##
## @return List with the first element being a data.frame with residuals from
## the BVLS models and the second element being a data.frame with coefficients
## from the BVLS models
extractBVLS <- function(BVLS, samples){

    ## extract residuals and coefficients from BVLS models
    modelResid <- do.call(plyr::rbind.fill, lapply(BVLS, function(x){
        y <- as.data.frame(t(x$residuals))
        colnames(y) <- attributes(x)$samples
        return(y)
    }))

    modelCoeff <- do.call(plyr::rbind.fill, lapply(BVLS, function(x){
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

## @title Running OLS models
##
## @description Function to run OLS on a list of model matrices, one element in
## the list should represent one peptide or protein.
##
## @return Returns complete OLS model
runOLS <- function(modelMat){
    modelRes <- lapply(modelMat, function(data){
        Y <- data$Y
        X <- data$X[, -1, drop=FALSE]
        stats::lm(Y ~ X)
    })
    return(modelRes)
}

## @title Extracting information from OLS models
##
## @description Function to extract the coefficients and p-values for each
## peptide/protein from the OLS models. If BVLS was run before the OLS, the
## p-values are estimated taken the degrees of freedom already used by the BVLS
## into account.
##
## @return List with the first element being a data.frame with coefficients from
## the OLS models and the second element being a data.frame with p- values from
## the OLS models
extractOLS <- function(OLS, formulaOLS, dfBVLS, coeffPval="Pr(>|t|)",
                       coeffTval="t value"){

    ## extract coefficients from OLS
    modelCoef <- do.call(plyr::rbind.fill, lapply(OLS, function(x){
        as.data.frame(t(stats::coef(x)))
    }))

    ## if BVLS was not run before, extract p-values directly from OLS
    if(is.null(dfBVLS)){
        modelPv <- do.call(plyr::rbind.fill, lapply(OLS, function(x){
            as.data.frame(t(summary(x)$coefficients[, coeffPval]))
        }))
    }
    ## if BVLS was run before, calculate p-values taking degrees of freedom
    ## used in BVLS models into account
    else{
        message("Estimating p-values while removing degrees of freedom
                consumed by BVLS.")
        modelPv <- calcualtePvalAfterBVLS(OLS, coeffTval, dfBVLS)
    }

    ## adjusting column names of coefficient and p-value data.frame
    if(ncol(modelCoef) == 2){
        formulaOLS <- stats::as.formula(formulaOLS)
        formulaVars <- formula.tools::get.vars(formulaOLS)[
            formula.tools::get.vars(formulaOLS)!="Y"]
        cols <- paste0(c("Intercept", formulaVars), "_OLS")
    }
    else{
        cols <- colnames(modelCoef)[-1]
        cols <- substring(cols, 2)
        cols <- paste0(c("Intercept", cols), "_OLS")
    }
    colnames(modelCoef) <- cols
    colnames(modelPv) <- cols

    return(list(modelCoeff=modelCoef, modelPv=modelPv))
}

## @title Estimate p-values taking degrees of freedom into account
##
## @description Function estimates p-values for all coefficients in the OLS
## model, but takes the degrees of freedom used when running the BVLS into
## account.
##
## @return Returns a data.frame with p-values
calcualtePvalAfterBVLS <- function(LM,coeffTval, dfBVLS){
    modelPv <- do.call(plyr::rbind.fill.matrix, lapply(LM, function(x){
        x <- stats::summary.lm(x)
        df <- x$df[2] - dfBVLS

        if(df<1){
            warning("Not enough degrees of freedom to estimate p-values,
                    model is not reliable! Returning NAs for affected peptides/
                    proteins.")
            pv <- stats::setNames(rep(NA, nrow(x$coefficients)),
                                  row.names(x$coefficients))
            return(t(pv))
        }

        ## Estimating p-values from t-values of OLS taking degrees of freedom
        ## used in BVLS into account
        else{
            pv <- sapply(x$coefficients[, coeffTval], function(y){
                2*stats::pt(abs(y), df, lower.tail=FALSE)
            })
            return(t(pv))
        }
    }))

    return(modelPv)
}
