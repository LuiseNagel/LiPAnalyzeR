# LiPAnalyzeR

Package for inferring peptide-specific changes in the PK accessibility from LiP-MS data. Provides funtions from reading in the data from Spectronaut reports over filtering to running different regression models to remove unwanted variation and carve out signal due to structural differences between conditions.

LiPAnalyzer was developed unter R 4.2.2 under Ubuntu 22.04.2 LTS.

## Installation

LiPAnalyzer requires R version 4.1.3. The following R commands allow you to install LiPAnalyzeR from github:
```
if (!require("devtools")) install.packages("devtools")
library(devtools)
if(!require("LiPAnalyzeR")) install_github("LuiseNagel/LiPAnalyzeR")
library(LiPAnalyzeR)
```

## Quickstart

### Getting Spectronaut data into R

To get started please export your LiP-MS data from Spectronaut to _.csv_ using the povided [Spectronaut schema](https://github.com/LuiseNagel/LiPAnalyzeR/blob/main/SpectroSchema_LipAnalyzerOut.rs).

Subsequently, read in your Spectronaut files.

```
LiPSpectro <- read.csv("/LiPSpectronautReprot.csv")
TrpSpectro <- read.csv("/LiPSpectronautReprot.csv")
```

### Creating peptide, protein and annotation matrices in R

First, the necessary peptide and protein matrices have to be create. This can be archieved with calling the function ```ExtractDataFromSpectro```. Per default, peptide(```PEP.Quantity```) and protein(```PG.Quantity```) quantities are extracted from the Spectronaut report. Alternatively, you can also choose to extract modified peptides or precursors instead of peptides by changing ```analysisLvl``` in the function.

```
QuantityList <- ExtractDataFromSpectro(spectroLiP = LiPSpectro,
                                       spectroTrp = TrpSpectro)
```

Subsequently, please ensure that the column names (sample annotations) are identical. This depends on if the LiP and Trp names where assigned the same names in Spectronaut itself. If this is not the case, please adjust them, making sure that they match.

Additionally, create annotation files of the peptides and proteins as well as the samples. 

```
annotPepProt <- GetPepProtAnnot(spectroOut = TrpSpectro, spectroOut2 = LiPSpectro)
annotSample <- GetSampleAnnot(spectroOut = TrpSpectro
```

Creating the sample annotation only works, if all necessary information is already added in Spectronaut. If this is not the case, please create the sample annotation file by hand, making sure that the rows are samples, with the row names matching the column names of all ```QuantityList``` matrices.


### Filtering and log-transforming peptide and protein quantities

Before models are fitted to the data should be log-transformed as well as filtered. Per default all peptide all log2 transformed quantities below 10 are set to NA, this setting can be accessed wiht ```thresholdMinLogQuant```. All peptides with NAs will be removed by the filtering. This is a very strict filter and can easily be adjusted in the function by defining the number of NAs allowed per condition ```maxNAperCondition```. If too many NAs are allowed the results will not be very reliable or models cannot be fitted due to a lack of degrees of freedom available. How strict the NA filter is set can therefore be very dataset dependent. 

```
QuantityList <- PreprocessQuantityMatrix(SpectroList = QuantityList,
                                         annotPP = annotPepProt,
                                         annotS = annotSample)
```

### Running regression models

LiPAnalyzeR allows to remove unwanted variation from the LiP peptide quantities, carving out the structural signal and provides the additional option to subsequently inferring differences in the structural proteome in between conditions. Additionally, the coefficients of the structural effect as well as the p-values can be extraced from the results, including the option for p-value correction.

```
modelStrucVar <- AnalyzeLiPPepData(spectroList = QuantityList, 
                                   annotS = annotSample)
resStrucVar <- SummarizeModelResults(resModel = modelStrucVar, 
                                     evalCovariable = "Condition_OLS")
```






