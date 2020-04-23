## The code assumes that the working directory is "drugSNF" folder
rm(list=ls())

source('RCode/preprocessInput.R')
source('RCode/sensitivityData.R')
source('RCode/perturbationData.R')
source('RCode/structureData.R')
source('RCode/constStructureLayer.R')
source('RCode/constSensitivityLayer.R')
source('RCode/constPerturbationLayer.R')
source('RCode/integrateStrctSensPert.R')
source('RCode/drugTargetBench.R')
source('RCode/ATCbench.R')
source('RCode/generateDrugPairs.R')
source('RCode/compConcordIndx.R')
source('RCode/generateRocPlot.R')
source('RCode/generatePRPlot.R')
source('RCode/predPerf.R')
source('RCode/communityGen.R')
source('RCode/cindexComp2.R')
source('RCode/utilities.R')
source('Rcode/RSNF.R')


library(PharmacoGx) 
library(apcluster)
library(rcdk)
library(fingerprint)
library(annotate)
library(org.Hs.eg.db)
library(SNFtool)
library(ROCR)
library(survcomp)
library(reshape2)
library(proxy)
library(OptimalCutpoints)
library(VennDiagram)
library(devEMF)


badchars <- "[\xb5]|[\n]|[,]|[;]|[:]|[-]|[+]|[*]|[%]|[$]|[#]|[{]|[}]|[[]|[]]|[|]|[\\^]|[/]|[\\]|[.]|[_]|[ ]"

## creating the output directory if not exists
if ( ! file.exists("Output")) {
  dir.create("Output")
}
#datestr <- date()
datestr <- Sys.Date()
mainDir <- getwd()
datestr <- gsub(badchars, "",  datestr)
subDir<-paste("Output/",datestr,sep = "")
options(someUniqueTag.mainDir = mainDir)
options(someUniqueTag.subDir = "subDir")

if(!file_test("-d", file.path(mainDir, subDir))){
  if(file_test("-f", file.path(mainDir, subDir))) {
    stop("Path can't be created because a file with that name already exists.")
  } else {
    dir.create(file.path(mainDir, subDir))
  }
}

cDrugs <- preprocessInput(dname="nci60", "lincs",venn_plot=FALSE)
dim(cDrugs$lincsboth)  ##  238 X 28
dim(cDrugs$nciboth)  ##  238 X 66

df.druginfo <- read.csv("Data/nci60fulldruginfo.csv")

drug.FDA <- unlist(df.druginfo[df.druginfo$FDA_APPROVED==TRUE,]$MOLECULE_NAME)

cDrugs$lincsboth <- cDrugs$lincsboth[cDrugs$lincsboth$pert_iname %in% drug.FDA,]
cDrugs$nciboth <- cDrugs$nciboth[cDrugs$nciboth$Drug.name %in% drug.FDA,]


## loading and cleaning data for the layers 
sensData <- sensitivityData("nci60", cDrugs$nciboth)
dim(sensData)  ##[1]  60 238
pertData <- perturbationData("lincs", cDrugs$lincsboth, "nci60")
dim(pertData)  ##[1]  978 238
strcData <- structureData("lincs", cDrugs$lincsboth)  ## a vector
length(strcData) ## 238

#Reduce All Matrices to lowest common set of drugs across all 3
# Get 237 drugs now in the reduced sets
commonDrugs <- Reduce(intersect,list(sort(names(strcData)),sort(colnames(sensData)),
                                     sort(colnames(pertData))))

##Reduce the data to the set of common drugs across the three layers
strcData<- strcData[commonDrugs]
sensData <- sensData[,commonDrugs]
pertData <- pertData[,commonDrugs]
## filter out cDrugs$lincsboth "dataframe" as well accordingly
cDrugs$lincsboth <- cDrugs$lincsboth[cDrugs$lincsboth$pert_iname %in% commonDrugs, ]
#

#Sanity Checks
if (ncol(sensData) != ncol(pertData)) stop(sprintf("error!"))
if (ncol(sensData) !=length(strcData)) stop(sprintf("error!"))
if (all(colnames(pertData) != colnames(sensData))) stop(sprintf("error!"))
if (all(colnames(pertData) != names(strcData))) stop(sprintf("error!"))


## network layer construction and integration by SNF
strcAffMat <- constStructureLayer(strcData)
sensAffMat <- constSensitivityLayer(sensData)
pertAffMat <- constPerturbationLayer(pertData)

strcAffMat <- normalize(strcAffMat)
sensAffMat <- normalize(sensAffMat)
pertAffMat <- normalize(pertAffMat)



integrtStrctSensPert <- integrateStrctSensPert_RNCE(sensAffMat,strcAffMat, pertAffMat,k1=24,k2=6,saverr=FALSE)
integrtDNF <- integrateStrctSensPert(sensAffMat,strcAffMat, pertAffMat)
integrtCISNF <- integrateStrctSensPert_cisnf(sensAffMat,strcAffMat, pertAffMat,k1=24,k2=6)



# Sanity Check - should all have the same dimensions: 237 X 237 --> 238...
dim(strcAffMat)
dim(sensAffMat)
dim(pertAffMat)
dim(integrtStrctSensPert)



## Benchmarking and validation
## 1- DRUG-TARGET 
dataBench1.5 <- drugTargetBench("uniprot", commonDrugs)
dim(dataBench1.5) ##[1] 86

pairs1 <- generateDrugPairs(dataBench1.5, strcAffMat, sensAffMat, pertAffMat, integrtStrctSensPert, integrtStrctSensPert, integrtStrctSensPert, NULL, integrtStrctSensPert,list(integrtDNF,integrtCISNF))


## ROC and PR plots
cutoffCHEMBL<-customRocPlot(pairs1, d1Name="nci60", d2Name="FDA", benchNam="drug-target(UNIPROT)",datestr=datestr,list("DNF","CISNF"))
customPRPlot(pairs1, d1Name="nci60", d2Name="FDA", benchNam="drug-target(UNIPROT)",datestr=datestr,list("DNF","CISNF"))


## 2- ATC
dataBench3 <- ATCBench("chembl-new", cDrugs$lincsboth,benchname="NCI60")
dim(dataBench3) ##[1] 72 72

pairs2 <- generateDrugPairs(dataBench3, strcAffMat, sensAffMat, pertAffMat, integrtStrctSensPert, integrtStrctSensPert, integrtStrctSensPert, integrtStrctSensPert, NULL,list(integrtDNF,integrtCISNF))

cutoffATC<-customRocPlot(pairs2, d1Name="nci60", d2Name="FDA", benchNam="ATC(CHEMBL)-Zscore",datestr=datestr,list("DNF","CISNF"))
customPRPlot(pairs2, d1Name="nci60", d2Name="FDA", benchNam="ATC(CHEMBL)-Zscore",datestr=datestr,list("DNF","CISNF"))
