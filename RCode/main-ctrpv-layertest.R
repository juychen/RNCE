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
library(rcdk)
library(fingerprint)
library(annotate)
library(org.Hs.eg.db)
library(SNFtool)
library(ROCR)
library(survcomp)
library(reshape2)
library(proxy)
library(PRROC)
library(apcluster)
library(OptimalCutpoints)
library(VennDiagram)
library(devEMF)


badchars <- "[\xb5]|[\n]|[,]|[;]|[:]|[-]|[+]|[*]|[%]|[$]|[#]|[{]|[}]|[[]|[]]|[|]|[\\^]|[/]|[\\]|[.]|[_]|[ ]"

## creating the output directory if not exists
#outputDir <- paste(getwd(), "/outputDir", sep="")
if ( ! file.exists("Output")) {
    dir.create("Output")
}
#datestr <- date()
datestr <- Sys.Date()
simestr <- Sys.Date()
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


load("Data/pertSimilarityMat-ctrpv2.RData")
load("Data/sensSimilarityMat-ctrpv2.RData")
load("Data/strcSimilarityMat-ctrpv2.RData")

load("Data/integrationSimilarityMat-ctrpv2.RData")



## Benchmarking and validation
## 1- DRUG-TARGET 
## loading and cleaning benchmark dataset
load("Output/dataBench-ctrpv2.RData")
load("Data/averageIorioPGX-all.RData") ## drug similarity matrix calculated based on Iorio snd Iskar et al. score (see iorioDIPS_PGX.R)
load("Data/averageIskarFinal.RData") ##load "iskar" results here ... (see iskar.R)
load("Data/drugERankSimil-CTRPV2-kendall.Rdata")

## Store the vairable to avoid the load to corrpute the file
storeintegrtStrctSensPert <- integrtStrctSensPert
#integrCISNFLoad <- load("Data/ctrpv2-CISNFIntegrated.RData")
contexLoad <- load("Data/rerkSimilarityMat-CTRP.Rdata")
contex <- get(contexLoad)

integrRSNFLoad <- load("Data/RSNF-ctrpv2.RData")
integrRSNF <-get(integrRSNFLoad)

integrRNCILoad <- load("Data/RNCI-ctrpv2.RData")
integrRNCI <-get(integrRNCILoad)

integrtStrctSensPert <- storeintegrtStrctSensPert

pairs <- generateDrugPairs(dataBench, strcAffMat, sensAffMat, pertAffMat, integrtStrctSensPert, average, finalIskarScore, NULL, drugERankSimil2,list(integrRSNF,contex,integrRNCI))


cutoffCML<-customRocPlot(pairs, d1Name="ctrpv2", d2Name="layertest", benchNam="drug-target",datestr=datestr,dnfnames=list("KRNN-SNF","KRNN-context","Both"))
customPRPlot(pairs, d1Name="ctrpv2", d2Name="layertest", benchNam="drug-target",datestr=datestr,dnfnames=list("KRNN-SNF","KRNN-context","Both"))

## 2- ATC
##load "superPred" post-processed results here ... (see compareTo_suprPred.R)
load("Data/superPredSimil-ctrp-kendall.Rdata")
load("Data/averageIorioPGX-all.RData") ## drug similarity matrix calculated based on Iorio snd Iskar et al. score (see iorioDIPS_PGX.R)
load("Data/averageIskarFinal.RData") ##load "iskar" results here ... (see iskar.R)
load("Output/ATCBench-ctrpv2.RData")

pairs2 <- generateDrugPairs(dataBench3, strcAffMat, sensAffMat, pertAffMat, integrtStrctSensPert, average, finalIskarScore, superPredSimil2, NULL,list(integrRSNF,contex,integrRNCI))

cutoffATC<-customRocPlot(pairs2, d1Name="ctrpv2", d2Name="layertest", benchNam="ATC(CHEMBL)",datestr=datestr,dnfnames=list("KRNN-SNF","KRNN-context","Both"))
customPRPlot(pairs2, d1Name="ctrpv2", d2Name="layertest", benchNam="ATC(CHEMBL)",datestr=datestr,dnfnames=list("KRNN-SNF","KRNN-context","Both"))

