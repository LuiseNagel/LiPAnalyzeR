# LiPAnalyzeR

LiPAnalyzeR is a R package for inferring structural changes from *Li*mited-*P*roteolysis *M*ass *S*pectrometry (LiP-MS) data by identitying peptide-specific changes in the PK accessibility of LiP peptides. We provide functions to import MS reports, e.g. from Spectronaut searches, preprocess and filter the data, remove of unwated variation (RUV) through constrianed regression from the LiP signal, inferring effect sizes for structural variation between conditions with corresponding p-values and plotting relevant results. 

LiPAnalyzer was developed in R 4.2.2 under Ubuntu 22.04.2 LTS.

## Installation

LiPAnalyzeR requires R version >=4.2.0. Please use the following R commands to install LiPAnalyzeR from github:
```
if(!require("devtools")) install.packages("devtools")
library(devtools)
if(!require("LiPAnalyzeR")) install_github("LuiseNagel/LiPAnalyzeR")
library(LiPAnalyzeR)
```

## Quickstart

### Import MS reports to R

The LiP-MS data have to be provided as MS report(s) (long-table format) for example from Spectronaut or XXX. 

The MS report has to include the following columns:
- sample names or identifiers
- quantity to analyze (e.g. peptide intensities)
- names/identifiers of the quantity to be analyzes (e.g. AS sequence of peptides)
- peptide names (i.e., AS sequence, may be identical to the names/identifiers of the quantity to be analyzed) 
- protein quantities
- protein names (matching the protein quantities)
- information the peptide is full- or half-tryptic
- start position of the peptide in the protein sequence

Optional columns to include are:
- all protein names that match the peptide sequence
- number of missed cleavages
- information if a peptide is proteotypic or not

If you are using Spectronaut, it is recommended to use the provided [Spectronaut schema](https://github.com/LuiseNagel/LiPAnalyzeR/blob/main/SpectroSchema_LipAnalyzerOut.rs) for exporting your data from Spectronaut. You can then use the Spectronaut-specific functions provided in LiPAnalyzeR for sorting the data.


Subsequently, read in your MS reports.

```
reportLiP <- read.csv("/reportLiP.csv")
reportTrP <- read.csv("/reprotTrP.csv")
```

If LiP & TrP are in a combined report, please seperate them into two data matrices after importing the data in R.


### Create quantity matrices

In a first step the MS quantity matrices have to be created using ```extractMSData``` (or ```extractSpectroData```). 

```
quantityList <- extractMSData(reportLiP = reportLiP, 
			      reportTrP = reportTrP,
			      sampleName =  _column in MS reports providing sample names_,
			      quantName =  _column in MS reports providing quantity identifiers_,
			      quantValue =  _column in MS reports providing quantities of interest_,
			      protValue =  _column in MS reports providing protein quantities_)
```


If you exported your MS data using the provided Spectronaut schema you can use the following functin instead:

```
quantityList <- extractSpectroData(reportLiP = reportLiP, 
			           reportTrP = reportTrP,
                                   analysisLvl = "Peptide")
```

```analysisLvl = "Peptide"``` will pick 'PEP.Quantity' as the column to extract the ```quantValue```, for analyzing peptide intensities. Alternatively ```analysisLvl``` can also be set to ```modifiedPeptide``` or ```Precursor```. Protein quantities will be extracted from the 'PG.Quantity' column of the Spectronaut report.



Subsequently, please ensure that the column names (sample annotations) of all list elements in the ```quantityList``` data are identical. This will be the case if the sample names provided in the LiP and TrP MS reports were already identical. If this is not the case, please adjust them to match. It is important, that every LiP sample has a (biological) matching TrP sample provided in the data. Samples without matching data will be removed in the following analysis steps.


### Create peptide & protein as well as sample annotation files

Create annotation files for the quantity of interest, providng information on matching peptides, proteins, start position of the peptides in the protein(s) and tryptic status of each peptide:

```
annotPepProt <- getPepProtAnnot(reportOut = reportLiP, 
				reportOut2 = reportTrP,
                            	quantName,
                            	pepName, #may be identical to quantName
                            	protName,
                            	isTryptic,
                            	startPosition)

```

Create an annotation file including sample information:

```
annotSample <- getSampleAnnot(spectroOut = reportTrP,
			      sampleName = "R.FileName",
			      sampleCondition = "R.Condition")
```

Creating the sample annotation with this function only works, if all necessary information is already added in the MS report. As this is commonly not the case, please create the sample annotation file by hand, making sure that the rows are samples, with the row names matching the column names of all ```quantityList``` matrices.


### Preprocessing and filtering LiP & TrP quantity matrices

Data should be log-transformed and filtered. Per default all peptides are log2 transformed and quantities below 10 are set to NA, this setting can be accessed wiht ```thresholdMinLogQuant```. Additionally, it is important to limit the number of NAs measured per peptide. Per default LiPAnalyzeR removes all peptides with at least one NA. Depending on the experimental set-up, this may not be appropriate for you data, you can alter these settings and thresholds in the function call. The function per default also filters for full-trytpic peptides, expecting these to be annotated as ```Specific``` in the ```annotPP``` data. If you want to run different modes than the ```default``` mode, such as ```HTonly```, ```FTHTjoin``` or ```LiPonly``` please check the specific descriptions in the function call. 
```
quantityList <- preprocessQuantityMatrix(SpectroList = quantityList,
                                         annotPP = annotPepProt,
                                         annotS = annotSample)
```
All quantity matrices in the preprocessed ```quantityList``` contain the same features and samples. Therefore, if sample naming was not adjusted to be identical between the TrP and LiP quantities all samples will be removed. 


### Running RUV & contrast models

LiPAnalyzeR allows to remove unwanted variation from the LiP peptide quantities, carving out the structural signal and provides the additional option to subsequently inferring differences in the structural proteome in between conditions. Additionally, the coefficients of the structural effect as well as the p-values can be extraced from the results, including the option for p-value correction.

```
modelStrucVar <- analyzeLiPPepData(spectroList = QuantityList, 
                                   annotS = annotSample)
resStrucVar <- summarizeModelResults(resModel = modelStrucVar, 
                                     evalCovariable = "Condition_RUV")
```
