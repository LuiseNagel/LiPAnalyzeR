library(testthat)
library(LiPAnalyzeR)

# Test for extractMSData function
test_that("extractMSData function returns expected output", {

    # Create mock MS reports
    reportLiP <- data.frame(SampleName = c("S1", "S2"),
                            Peptide = c("Pep1", "Pep2"),
                            PepQuant = c(1, 2),
                            ProtQuant = c(10, 20))

    reportTrp <- data.frame(SampleName = c("S1", "S2"),
                            Peptide = c("Pep3", "Pep4"),
                            PepQuant = c(3, 4),
                            ProtQuant = c(30, 40))

    # Test with LiPonly = FALSE
    result <- extractMSData(reportLiP = reportLiP,
                            reportTrp = reportTrp,
                            sampleName = "SampleName",
                            quantName = "Peptide",
                            quantValue = "PepQuant",
                            protValue = "ProtQuant")
    expect_equal(length(result), 3)

    # Test with LiPonly = TRUE
    result <- extractMSData(reportLiP = reportLiP,
                            sampleName = "SampleName",
                            quantName = "Peptide",
                            quantValue = "PepQuant",
                            protValue = "ProtQuant",
                            LiPonly = TRUE)
    expect_equal(length(result), 2)
})

# Test for convert2Matrix function
test_that("convert2Matrix function returns expected output", {

    # Create mock MS report
    report <- data.frame(Peptide = c("Pep1", "Pep2", "Pep1", "Pep2"),
                         Sample = c("S1", "S1", "S2", "S2"),
                         quantValue = c(1, 2, 3, 4))

    # Test the function
    result <- convert2Matrix(reportOut = report,
                             quantValue = "quantValue",
                             rows = "Peptide",
                             cols = "Sample")
    expect_equal(dim(result), c(2, 2))
    expect_equal(result["Pep1","S1"], 1)
})


# Test for getPepProtAnnot function
test_that("getPepProtAnnot function returns expected output", {

    # Create mock MS report
    report <- data.frame(Peptide = c("ABCD", "EFGHIJKLM"),
                         Protein = c("ProtA", "ProtB"),
                         isTryptic = c("FT", "HT"),
                         startPosition = c(1, 10))

    # Test the function
    result <- getPepProtAnnot(reportOut = report,
                              reportOut2 = report,
                              quantName = "Peptide",
                              pepName = "Peptide",
                              protName = "Protein",
                              isTryptic = "isTryptic",
                              startPosition = "startPosition")

    expect_equal(dim(result), c(2,6))
    expect_equal(row.names(result), result[, "quantID"])
    expect_type(result[,"endPosition"], "character")
})


# Test for getSampleAnnot function
test_that("getSampleAnnot function returns expected output", {

    # Create mock MS report with categorical and continious condition
    report <- data.frame(SampleName = c("S1", "S2", "S3", "S4"),
                         ConditionFactor = c("CX", "CY", "CZ", "CY"),
                         ConditionContinuous = c(1, 4, 7, 3))

    # Test default settings (CY condition is reference, dummy coding)
    result <- getSampleAnnot(report,
                             sampleName = "SampleName",
                             sampleCondition = "ConditionFactor",
                             typeCondition = "factor",
                             baseLevel = "CY")
    lm <- lm(report$ConditionContinuous ~ result$Condition)

    # Check resulting structure and used reference level in lm
    expect_equal(dim(result), c(4,2))
    expect_false("resultCondition2" %in% names(lm$coefficients))

    # Test default settings (first condition is reference, wec coding)
    result <- getSampleAnnot(report,
                             sampleName = "SampleName",
                             sampleCondition = "ConditionFactor",
                             typeCondition = "factor",
                             contrastCoding = "weccoding")

    # Check that wec coding worked as expected
    expect_equal(as.vector(contrasts(result$Condition)[1,]), c(-2, -1))

    # Test continious variable
    result <- getSampleAnnot(report,
                             sampleName = "SampleName",
                             sampleCondition = "ConditionContinuous",
                             typeCondition = "continuous")

    # Check that variable is continuous
    expect_type(result$Condition, "numeric")
})
