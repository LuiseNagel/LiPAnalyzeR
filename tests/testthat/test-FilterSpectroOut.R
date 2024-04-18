library(testthat)

# Test for preprocessQuantityMatrix function
test_that("preprocessQuantityMatrix returns a list of matrices", {

    # Mock input data
    set.seed(42)
    LiPPep <- matrix(rnorm(100, mean=14), nrow=10)
    TrPPep <- matrix(rnorm(100, mean=15), nrow=10)
    TrPProt <- matrix(rep(rnorm(50, mean=16),2), byrow = TRUE , nrow=10)
    quantityList <- list(LiPPep=LiPPep, TrPPep=TrPPep, TrPProt=TrPProt)
    quantityList <- lapply(quantityList, \(x){
        row.names(x) <- paste0("Peptide_", seq(1, 10))
        colnames(x) <- paste0("Sample_", seq(1, 10))
        return(x)
    })

    # Mock annotPP data
    annotPP <- data.frame(quantID=row.names(quantityList$LiPPep),
                          Protein=paste0("Protein", rep(c(1:5), 2)),
                          row.names=row.names(quantityList$LiPPep))

    # Call the function
    processedMatrix <- preprocessQuantityMatrix(quantityList=quantityList,
                                                annotPP=annotPP,
                                                logT = FALSE,
                                                filterNA="none",
                                                filterTryptic=FALSE)

    # Check the output
    expect_type(processedMatrix, "list")
    expect_equal(length(processedMatrix), length(quantityList))
    expect_equal(dim(processedMatrix$LiPPep), dim(LiPPep))
})

test_that("preprocessQuantityMatrix handles missing data correctly", {
    # Mock input data with missing values
    set.seed(42)
    LiPPep <- matrix(c(rnorm(50, mean=14), rep(NA, 50)), byrow=TRUE, nrow=10)
    TrPPep <- matrix(rnorm(100, mean=15), nrow=10)
    TrPProt <- matrix(rep(rnorm(50, mean=16),2), byrow = TRUE , nrow=10)
    quantityList <- list(LiPPep=LiPPep, TrPPep=TrPPep, TrPProt=TrPProt)
    quantityList <- lapply(quantityList, \(x){
        row.names(x) <- paste0("Peptide_", seq(1, 10))
        colnames(x) <- paste0("Sample_", seq(1, 10))
        return(x)
    })

    # Mock annotPP data
    annotPP <- data.frame(quantID=row.names(quantityList$LiPPep),
                          Protein=paste0("Protein", rep(c(1:5), 2)),
                          row.names=row.names(quantityList$LiPPep))

    # Call the function
    processedMatrix <- preprocessQuantityMatrix(quantityList=quantityList,
                                                annotPP=annotPP,
                                                logT = FALSE,
                                                filterNA="all",
                                                filterTryptic=FALSE)

    # Check the output
    expect_true(all(rowSums(is.na(processedMatrix$LiPPep)) == 0))
    expect_equal(unname(unlist(lapply(processedMatrix, dim))), rep(c(5, 10), 3))
})

test_that("preprocessQuantityMatrix calls errors if wrong input is given", {
    # Mock input data
    set.seed(42)
    LiPPep <- matrix(rnorm(100, mean=14), nrow=10)
    TrPPep <- matrix(rnorm(100, mean=15), nrow=10)
    TrPProt <- matrix(rep(rnorm(50, mean=16),2), byrow = TRUE , nrow=10)
    quantityList <- list(LiPPep=LiPPep, TrPPep=TrPPep, TrPProt=TrPProt)
    quantityList <- lapply(quantityList, \(x){
        row.names(x) <- paste0("Peptide_", seq(1, 10))
        colnames(x) <- paste0("Sample_", seq(1, 10))
        return(x)
    })

    # Mock annotPP data
    annotPP <- data.frame(quantID=row.names(quantityList$LiPPep),
                          Protein=paste0("Protein", rep(c(1:5), 2)),
                          row.names=row.names(quantityList$LiPPep))

    # Expecting errors when calling functions
    ## invalide mode
    expect_error(preprocessQuantityMatrix(quantityList=quantityList,
                                          annotPP=annotPP,
                                          mode="invalid_mode"))
    ## quantityList input has wrong format
    expect_error(preprocessQuantityMatrix(quantityList=LiPPep,
                                          annotPP=annotPP))
    ## mode requires annotPP to be provided
    expect_error(preprocessQuantityMatrix(quantityList=quantityList,
                                          mode = "FTHTjoin"))

})
