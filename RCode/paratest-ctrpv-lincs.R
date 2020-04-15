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
source('RCode/generateRocPlot.R')
source('RCode/generatePRPlot.R')
source('RCode/predPerf.R')
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
library(PRROC)
library(apcluster)
library(OptimalCutpoints)
library(VennDiagram)
library(devEMF)
library(umap)


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

params_k1 <- c(8,16,32,64,128,192)
params_k2 <- c(4,8,16,32)

# Find common drugs between CTRPV2 and LINCS dataset
cDrugs <- preprocessInput(dname="ctrpv2", "lincs",venn_plot=TRUE)
dim(cDrugs)  ##239 X 28

# Process Sensitivity, Perturbation, and Structure layers for set of common drugs
sensData <- sensitivityData("ctrpv2", cDrugs)  ## 645 X 239
dim(sensData)
pertData <- perturbationData("lincs", cDrugs, "ctrpv2")  ## 978 X 239
dim(pertData)
strcData <- structureData("lincs", cDrugs)  ## a vector  --> 239 elemnts
length(strcData)

#
## Get the common drugs (239) among the 3 datasets/layers
commonDrugs <- Reduce(intersect,list(sort(names(strcData)),sort(colnames(sensData)),
                                      sort(colnames(pertData))))
length(commonDrugs) ## 239 ..

strcData<- strcData[commonDrugs] # 239 drugs
sensData <- sensData[,commonDrugs] # 645 x 239 drugs
pertData<- pertData[,commonDrugs] #978 genes x  239 

## network layer construction and integration by SNF
strcAffMat <- constStructureLayer(strcData)
sensAffMat <- constSensitivityLayer(sensData)
pertAffMat <- constPerturbationLayer(pertData)


strcAffMat <- (normalize(strcAffMat))
sensAffMat <- (normalize(sensAffMat))
pertAffMat <- (normalize(pertAffMat))


## Benchmarking and validation
## 1- DRUG-TARGET 
## loading and cleaning benchmark dataset
dataBench <- drugTargetBench("ctrpv",  commonDrugs) # 141 x 141 drug-drug adjacency matrix --> 141
## 2- ATC
dataBench3 <- ATCBench("chembl-new", cDrugs)


ROCS.dtg <- list()
PRCS.dtg <- list()
ROCS.atc <- list()
PRCS.atc <- list()
paras <- list()

for(k1 in params_k1){
  
  for(k2 in params_k2){
    
    if (k2>=k1) {
      next
    }else{
      
  
      integrtStrctSensPert <- integrateStrctSensPert_RNCE(sensAffMat,strcAffMat, pertAffMat,k1=k1,k2=k2,saverr=FALSE)
      

      save(integrtStrctSensPert, file=paste("Data/integration",k1,k2,"ctrpv2.RData",sep="-"))
      

      pairs <- generateIntDrugPairs(dataBench, integrtStrctSensPert)
      
      
      ## ROC and PR plots
      AUROC.dtg<-predPerf(pairs$integrPairs$obs.combiall, pairs$benchPairs$bench)
      AUPRC.dtg<-predPerf(pairs$integrPairs$obs.combiall, pairs$benchPairs$bench,plotType="PR")
      
      ROCS.dtg<-append(ROCS.dtg,AUROC.dtg$auc)
      PRCS.dtg<-append(PRCS.dtg,AUPRC.dtg$auc.integral)
      
      
      pairs2 <- generateIntDrugPairs(dataBench3, integrtStrctSensPert)
      
      ## ROC and PR plots
      AUROC.atc<-predPerf(pairs2$integrPairs$obs.combiall, pairs2$benchPairs$bench)
      AUPRC.atc<-predPerf(pairs2$integrPairs$obs.combiall, pairs2$benchPairs$bench,plotType="PR")
      
      ROCS.atc<-append(ROCS.atc,AUROC.atc$auc)
      PRCS.atc<-append(PRCS.atc,AUPRC.atc$auc.integral)
      paras<-append(paras,paste(k1,k2,sep=","))
      paste(k1,k2,sep=",")
    }
  }
}

df.result <- data.frame(K=unlist(paras)
  ,ROCchembl=unlist(ROCS.dtg)
                        ,PRCchembl=unlist(PRCS.dtg)
                        ,ROCATC=unlist(ROCS.atc)
                        ,PRCATC=unlist(PRCS.atc))
write.csv(df.result, file = paste(mainDir,subDir,"ctrp_paratest.csv",sep = "/"))

