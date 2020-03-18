# RNCE: Network Integration with Reciprocal Neighbors Contextual Encoding for Multi-modal Drug Community Study


The following steps will reproduce the output files (figure, tables...) mentioned in the main manuscript.
The script will be using data files such as:

## Files to run the scripts 

- Drug-target (benchmark) from CHEMBL and CTRP
- Sensitivity measurements (drugs x cell lines) from CTRPv2 and NCI-60
- Transcriptomic data from the LINCS database http://lincs.hms.harvard.edu/ created using PharmacoGx package https://cran.r-project.org/web/packages/PharmacoGx/index.html

## Run the R scripts 

* main.R 
(this will execute the following R codes):

```
## The code assumes that the working directory is "RNCE" folder
rm(list=ls())

#Generate comparing networks results
source("RCode/generatecompare.R")

#Test paramerers
source("RCode/paratest-ctrpv-lincs.R")
source("RCode/paratest-nci60-lincs.R")

#Generate networks and generated ROC PRC plots
source("RCode/main-ctrpv-lincs.R")
source("RCode/main-nci60-lincs.R")

#Test ARI of community detection and test the number of communities
source("RCode/communityARI.R")

#Generate community plots and community targets
source("RCode/main-network-generation-ctrpv2.R")
source("RCode/main-network-generation-nci60.R")
source("RCode/communityTargets.R")
source("RCode/CommunityEnrichmentNovel.R")

#Genearte barplots of the correlation analysis
source("RCode/CorrelationAnalysis.R")

#Generate unit test plots
source("RCode/main-ctrpv-layertest.R")
source("RCode/main-nci60-layertest.R")


```

# Set up the software environment
The software is developed and tested in the environment:
```
R version 3.6.1 (2019-07-05) -- "Action of the Toes"
Copyright (C) 2019 The R Foundation for Statistical Computing
Platform: x86_64-w64-mingw32/x64 (64-bit)
```

Try to run the scrpit "Rcode/dependencyInstaller.R" to get the packages required to generate the results

