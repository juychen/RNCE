
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

load("Data/pertSimilarityMat-nci60.RData")
load("Data/sensSimilarityMat-nci60.RData")
load("Data/strcSimilarityMat-nci60.RData")
load("Data/integrationSimilarityMat-nci60.RData")







# Sanity Check - should all have the same dimensions: 237 X 237 --> 238...
dim(strcAffMat)
dim(sensAffMat)
dim(pertAffMat)
dim(integrtStrctSensPert)



## Benchmarking and validation
## 1- DRUG-TARGET 

load("Data/averageIorioPGX-all.RData") ## drug similarity matrix calculated based on Iorio snd Iskar et al. score (see iorioDIPS_PGX.R)
load("Data/averageIskarFinal.RData") ##load "iskar" results here ... (see iskar.R)
load("Data/drugERankSimil-nci60-kendall.Rdata")
load("Output/dataBench-nci60.RData")


## Store the vairable to avoid the load to corrpute the file
storeintegrtStrctSensPert <- integrtStrctSensPert
integrtContexLoad <- load("Data/rerkSimilarityMat-nci60.Rdata")
integrtContex <- get(integrtContexLoad)
integrRSNFLoad <- load("Data/RSNF-nci60.RData")
integrRSNF <-get(integrRSNFLoad)
integrRNCILoad <- load("Data/RNCI-nci60.RData")
integrRNCI <-get(integrRNCILoad)


## Network in each step
name.initial.nci60.rnce <- "Data/initSimilarityMat-nci60.Rdata"
load(name.initial.nci60.rnce)
integrIni<- integration

integrtStrctSensPert <- storeintegrtStrctSensPert

pairs1 <- generateDrugPairs(dataBench1.5, strcAffMat, sensAffMat, pertAffMat, integrtStrctSensPert, average, finalIskarScore, NULL, drugERankSimil2,list(integrRSNF,integrtContex,integrRNCI))

cutoffCHEMBL<-customRocPlot(pairs1, d1Name="nci60", d2Name="layertest", benchNam="drug-target(UNIPROT)",datestr=datestr,list("KRNN-SNF","KRNN-context","Both"))
customPRPlot(pairs1, d1Name="nci60", d2Name="layertest", benchNam="drug-target(UNIPROT)",datestr=datestr,list("KRNN-SNF","KRNN-context","Both"))

pairssteps <- generateDrugPairs(dataBench1.5, strcAffMat, sensAffMat, pertAffMat, integrtStrctSensPert, average, finalIskarScore, NULL, drugERankSimil2,list(integrRSNF,integrtContex,integrRNCI))

cutoffCHEMBL<-customRocPlot(pairs1, d1Name="nci60", d2Name="layertest", benchNam="drug-target(UNIPROT)",datestr=datestr,list("KRNN-SNF","KRNN-context","Both"))
customPRPlot(pairs1, d1Name="nci60", d2Name="layertest", benchNam="drug-target(UNIPROT)",datestr=datestr,list("KRNN-SNF","KRNN-context","Both"))


## 2- ATC

##load "superPred" post-processed results here ... (see compareTo_suprPred.R)
load("Data/superPredSimil-nci60-kendall.Rdata")
load("Data/averageIorioPGX-all.RData") ## drug similarity matrix calculated based on Iorio snd Iskar et al. score (see iorioDIPS_PGX.R)
##load "iskar" results here ... (see iskar.R)
load("Data/averageIskarFinal.RData")
load("Output/ATCBench-nci60.RData")

pairs2 <- generateDrugPairs(dataBench3, strcAffMat, sensAffMat, pertAffMat, integrtStrctSensPert, average, finalIskarScore, superPredSimil2, NULL,list(integrRSNF,integrtContex,integrRNCI))

## ROC and PR plots
cutoffATC<-customRocPlot(pairs2, d1Name="nci60", d2Name="layertest", benchNam="ATC(CHEMBL)-Zscore",datestr=datestr,list("KRNN-SNF","KRNN-context","Both"))
customPRPlot(pairs2, d1Name="nci60", d2Name="layertest", benchNam="ATC(CHEMBL)-Zscore",datestr=datestr,list("KRNN-SNF","KRNN-context","Both"))



