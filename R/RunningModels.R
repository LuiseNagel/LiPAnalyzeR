##matrixStats
globalVariables(names=c("Condition", "XPep", "XProt",  "Y"))

#' @title Fitting model for accessibility changes in LiP peptides.
#'
#' @description Function to build linear regression models for fitting MS data
#' and retrieving structural variation between different conditions.
#'
#' @usage  AnalyzeLiPPepData(spectroList, annotS, infoCondition="Condition",
#' formulaBVLS="Y~XPep+XProt", formulaOLS=NULL, lowBVLS=c(-1e9, 0, 0),
#' upBVLS=c(Inf, Inf, Inf), addBVLSbounds=FALSE, LiPonly=FALSE)
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
#' LiPonly version of the package and not providing trypsin-only data.
#'
#' @export
AnalyzeLiPPepData <- function(spectroList, annotS, infoCondition="Condition",
                              formulaBVLS="Y~XPep+XProt", formulaOLS=NULL,
                              lowBVLS=c(-1e9, 0, 0), upBVLS=c(Inf, Inf, Inf),
                              addBVLSbounds=FALSE, LiPonly=FALSE){
    if(is.null(formulaOLS)){
        formulaOLS <- paste0("Y~", infoCondition)
    }
    LiPOut <- RunModel(spectroList=spectroList,
                       annotS=annotS,
                       formulaBVLS=formulaBVLS,
                       formulaOLS=formulaOLS,
                       lowBVLS=lowBVLS,
                       upBVLS=upBVLS,
                       addBVLSbounds=addBVLSbounds,
                       LiPonly=LiPonly)
    return(LiPOut)
}

#' @title Fitting model for PK-independent changes in Trp peptides

#' @description Function to build linear regression models for fitting MS data
#' and retrieving PK-independent peptide variation between different conditions.
#'
#' @usage AnalyzeTrpPepData(spectroList, annotS, infoCondition="Condition",
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
AnalyzeTrpPepData <- function(spectroList, annotS, infoCondition="Condition",
                              formulaBVLS="Y~XProt", formulaOLS=NULL,
                              lowBVLS=c(-1e9, 0), upBVLS=c(Inf, Inf),
                              addBVLSbounds=FALSE){
    if(is.null(formulaOLS)){
        formulaOLS <- paste0("Y~", infoCondition)
    }
    spectroList <- list(Y=spectroList$TrpPep, XProt=spectroList$TrpProt)
    TrpOut <- RunModel(spectroList=spectroList,
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
#' @usage AnalyzeTrpProtData(spectroList, annotS, annotPP,
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
AnalyzeTrpProtData <- function(spectroList, annotS, annotPP,
                               infoCondition="Condition",
                               infoProtName="Protein", formulaBVLS=NULL,
                               formulaOLS=NULL, lowBVLS=NULL, upBVLS=NULL,
                               addBVLSbounds=FALSE, LiPonly=FALSE){
    if(is.null(formulaOLS)){
        formulaOLS <- paste0("Y~", infoCondition)
    }

    # map peptides to proteins
    if(!LiPonly){
        protData <- spectroList$TrpProt
    }
    else{
        protData <- spectroList$LiPProt
    }
    Peps2Prots <- split(row.names(protData), annotPP[row.names(protData),
                                                     infoProtName])
    table(unlist(lapply(Peps2Prots, length)))
    Peps2Prots <- unlist(lapply(Peps2Prots, \(x) (x[1])))
    spectroList <- list(Y=protData[Peps2Prots,])
    row.names(spectroList$Y) <- names(Peps2Prots)

    ProtOut <- RunModel(spectroList=spectroList,
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
#' @usage RunModel(spectroList, annotS=NULL, formulaBVLS="Y~XPep+XProt",
#' formulaOLS="Y~Condition", lowBVLS=c(-1e9, 0, 0), upBVLS=c(Inf, Inf, Inf),
#' addBVLSbounds=FALSE, returnBVLSmodels=FALSE, returnOLSmodels=FALSE,
#' LiPonly=FALSE)
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
#' LiPonly version of the package and not providing trypsin-only data.
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
RunModel <- function(spectroList, annotS=NULL, formulaBVLS="Y~XPep+XProt",
                     formulaOLS="Y~Condition", lowBVLS=c(-1e9, 0, 0),
                     upBVLS=c(Inf, Inf, Inf), addBVLSbounds=FALSE,
                     returnBVLSmodels=FALSE, returnOLSmodels=FALSE,
                     LiPonly=FALSE){

    # if necessary transforming formulas into character
    # remove spaces in formula
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

    # check input format of formulas
    if(!(is.null(formulaBVLS)|is.character(formulaBVLS))&
       (is.null(formulaOLS)|is.character(formulaOLS))){
        stop("Please provide 'formula' in the corrrect class. Use 'character' or
             'formula'.")
    }

    if(LiPonly){
        if(as.character(formulaBVLS) == "Y~XPep+XProt"){
            formulaBVLS <- "Y~XProt"
            lowBVLS <- c(-1e9, 0)
            upBVLS <- c(Inf, Inf)
        }
    }

    if(paste(names(spectroList), collapse="_") == c("LiPPep_TrpPep_TrpProt")){
        names(spectroList) <- c("Y", "XPep", "XProt")
    }

    else if(paste(names(spectroList), collapse="_") == c("LiPPep_LiPProt")){
        names(spectroList) <- c("Y", "XProt")
    }

    # Assuring row.names and colnames fit over all input data
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

    # returning warning message if no formulas are provided
    if(!is.null(formulaBVLS) & !is.null(formulaOLS)){
        stop("No BVLS or OLS formula provided. Please provide at least one of
             them to define the model(s) you want to run.")
    }
    resAll <- NULL

    # running BVLS
    if(!is.null(formulaBVLS)){
        message("Running bounded variable least square models.")
        modelMat <- CreateModelMatrix(spectroList, formulaBVLS, annotS, samples)
        modelBVLS <- RunBVLS(formulaBVLS, modelMat, lowBVLS, upBVLS,
                             addBVLSbounds)
        resBVLS <- ExtractBVLS(modelBVLS, samples)
        if(is.null(formulaOLS)){
            resALL <- resBVLS
            message("Returning BVLS results including residuals and estimated
                    coefficients.")
        }
        dfBVLS <- ncol(resBVLS[[2]])-1
        spectroList[["Y"]] <- resBVLS[[1]]
    }
    else{
        dfBVLS <- NULL
    }

    # running OLS
    if(!is.null(formulaOLS)){
        message("Running ordinary least square models.")
        modelMat <- CreateModelMatrix(spectroList, formulaOLS, annotS, samples)
        modelOLS <- RunOLS(modelMat)
        resOLS <- ExtractOLS(modelOLS, formulaOLS, dfBVLS)
        if(!is.null(formulaBVLS)){
            resOLS$modelCoeff <- cbind(resBVLS$modelCoeff,resOLS$modelCoeff)
        }
        resALL <- resOLS
    }

    # Adding feature names to output)
    resALL <- lapply(resALL, function(x){
        rownames(x) <- feat
        return(as.data.frame(x))
    })

    # Add message if there TrpPep or TrpProt coefficients are unrealistic
    if(any(grepl("XPep|XProt", colnames(resAll[[1]])))){
        if(max(stats::na.omit(c(resAll[[1]][, grepl("XPep",
                                                    colnames(resAll[[1]]))],
                                resAll[[1]][, grepl("XProt",
                                                    colnames(resAll[[1]]))])))
           >5){
            message("At least one peptide/protein coefficient is higher than 5,
                    please check the results of the effected peptides for
                    plausability.")
        }
    }

    # Adding models to results if required
    if(returnBVLSmodels){
        resALL[[length(resALL)+1]] <- list(modelBVLS)
        names(resALL)[length(resALL)] <- "modelBVLS"
    }

    if(returnOLSmodels){
        resALL[[length(resALL)+1]] <- list(modelOLS)
        names(resALL)[length(resALL)] <- "modelOLS"
    }

    # Potentially adding models
    return(resALL)
}


# Function for getting list of Y and model matrix for all peptides
CreateModelMatrix <- function(spectroList, formula, annotS, samples){
    #formula <- gsub(".*~", "~", formula)
    list2env(spectroList, envir=environment()) # write matrixes into environment
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


# Running BVLS, returning residuals and coeff
RunBVLS <- function(formula, modelMat, lowBVLS, upBVLS, addBVLSbounds){

    # bvls does not take -Inf as an input
    # change -Inf to high negative number instead
    lowBVLS[lowBVLS == -Inf] <- -1e9

    # run BVLS models
    modelRes <- lapply(modelMat, function(data){

        ### Add progress bar
        Y <- data$Y
        X <- data$X

        # if set to true, adding as many as necessary -Inf and Inf to the
        # boundaries of the BVLS
        if(addBVLSbounds){
            nX <- ncol(X)
            lowRep <- nX-length(lowBVLS)
            lowBVLS <- c(lowBVLS, rep(-1e9, lowRep))
            upRep <- nX-length(upBVLS)
            upBVLS <- c(upBVLS, rep(Inf, upRep))
        }

        # running BVLS
        BVLS <- bvls::bvls(A=X,
                           b=Y,
                           bl= lowBVLS,
                           bu=upBVLS)

        ## Add variable and sample annotation to BVLS output
        varsBVLS <- paste0(dimnames(X)[[2]], "_BVLS")
        varsBVLS[1] <- "Intercept_BVLS"
        attr(BVLS, "variables") <- varsBVLS
        attr(BVLS, "samples") <- attributes(X)$samples
        return(BVLS)
    })
    return(modelRes)
}

# Getting Ceoff from BVLS
ExtractBVLS <- function(BVLS, samples){

    # extract residuals and coefficients from BVLS models
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

    # adding NA columns in resid
    startDf <- as.data.frame(matrix(NA, nrow=1, ncol=length(samples)))
    colnames(startDf) <- samples
    modelResid <- plyr::rbind.fill(startDf, modelResid)[-1,]

    return(list(modelResid=modelResid, modelCoeff=modelCoeff))
}

# Running OLS
## Add warning for degrees of freedom
RunOLS <- function(modelMat){
    modelRes <- lapply(modelMat, function(data){
        ### Add progress bar
        Y <- data$Y
        X <- data$X[, -1, drop=FALSE]
        stats::lm(Y ~ X)
    })
    return(modelRes)
}


# Getting Ceoff and P-values from OLS
ExtractOLS <- function(OLS, formulaOLS, dfBVLS, coeffPval="Pr(>|t|)",
                       coeffTval="t value"){
    modelCoef <- do.call(plyr::rbind.fill, lapply(OLS, function(x){
        as.data.frame(t(stats::coef(x)))
    }))

    if(is.null(dfBVLS)){
        modelPv <- do.call(plyr::rbind.fill, lapply(OLS, function(x){
            as.data.frame(t(summary(x)$coefficients[, coeffPval]))
        }))
    }

    else{
        message("Estimating p-values while removing degrees of freedom
                consumed by BVLS.")
        modelPv <- CalcualtePvalAfterBVLS(OLS, coeffTval, dfBVLS)
    }
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

# Calculate p-value from OLS model taking degrees of freedom already considered
#in BVLS into account
CalcualtePvalAfterBVLS <- function(LM,coeffTval, dfBVLS){
    modelPv <- do.call(plyr::rbind.fill.matrix, lapply(LM, function(x){
        x <- stats::summary.lm(x)
        df <- x$df[2] - dfBVLS

        if(df<1){
            warning("Not enough degrees of freedom to estimate p-values,
                    model is not reliable!")
            pv <- stats::setNames(rep(NA, nrow(x$coefficients)),
                                  row.names(x$coefficients))
            return(t(pv))
        }
        else{
            pv <- sapply(x$coefficients[, coeffTval], function(y){
                2*stats::pt(abs(y), df, lower.tail=FALSE)
            })
            return(t(pv))
        }
    }))
    return(modelPv)
}
