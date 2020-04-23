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
library(PRROC)
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

# Find common drugs between CTRPV2 and LINCS dataset
cDrugs <- preprocessInput(dname="ctrpv2", "lincs",venn_plot=FALSE)
dim(cDrugs)  ##239 X 28

df.druginfo <- read.csv("Data/ctrpfulldruginfo.csv")

drug.FDA <- unlist(df.druginfo[df.druginfo$FDA_APPROVED==TRUE,]$MOLECULE_NAME)

cDrugs <- cDrugs[cDrugs$pert_iname %in% drug.FDA,]

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
pertAffMat <- constPerturbationLayerDNF(pertData)




strcAffMat <- (normalize(strcAffMat))
sensAffMat <- (normalize(sensAffMat))
pertAffMat <- (normalize(pertAffMat))

# save(pertAffMat, file="Data/pertSimilarityMat-ctrpv2.RData")
# save(sensAffMat, file="Data/sensSimilarityMat-ctrpv2.RData")
# save(strcAffMat, file="Data/strcSimilarityMat-ctrpv2.RData")

integrtStrctSensPert <- integrateStrctSensPert_RNCE(sensAffMat,strcAffMat, pertAffMat,k1=8,k2=2,saverr=FALSE)
intDNF <- integrateStrctSensPert(sensAffMat,strcAffMat, pertAffMat)
intCISNF <- integrateStrctSensPert_cisnf(sensAffMat,strcAffMat, pertAffMat,k1=16,k2=3)



#save(integrtStrctSensPert, file="Data/integrationSimilarityMat-ctrpv2.RData")



## Benchmarking and validation
## 1- DRUG-TARGET 
## loading and cleaning benchmark dataset
dataBench <- drugTargetBench("ctrpv",  commonDrugs) # 141 x 141 drug-drug adjacency matrix --> 141


pairs <- generateDrugPairs(dataBench, strcAffMat, sensAffMat, pertAffMat, integrtStrctSensPert, 
                           integrtStrctSensPert, integrtStrctSensPert, NULL, integrtStrctSensPert,list(intDNF,intCISNF))


cutoffCML<-customRocPlot(pairs, d1Name="ctrpv2", d2Name="FDA", benchNam="drug-target",datestr=datestr,list("DNF","CISNF"))
customPRPlot(pairs, d1Name="ctrpv2", d2Name="FDA", benchNam="drug-target",datestr=datestr,list("DNF","CISNF"))


# 
# ## 2- ATC
dataBench3 <- ATCBench("chembl-new", cDrugs)
dim(dataBench3) ##[1]  51 51

pairs2 <- generateDrugPairs(dataBench3, strcAffMat, sensAffMat, pertAffMat, integrtStrctSensPert, integrtStrctSensPert, integrtStrctSensPert, integrtStrctSensPert, NULL,list(intDNF,intCISNF))


## ROC and PR plots
cutoffATC<-customRocPlot(pairs2, d1Name="ctrpv2", d2Name="FDA", benchNam="ATC(CHEMBL)",datestr=datestr,list("DNF","CISNF"))
customPRPlot(pairs2, d1Name="ctrpv2", d2Name="FDA", benchNam="ATC(CHEMBL)",datestr=datestr,list("DNF","CISNF"))

# dtgbench.ctrp <- dataBench
# name.ctrp.dtg <- rownames(dtgbench.ctrp)
# graph.dtgbench.ctrp <- graph.adjacency(dtgbench.ctrp, mode="undirected",diag=FALSE)
# ceb.dtgbench.ctrp<-cluster_louvain(graph.dtgbench.ctrp)
# loulabel.dtgbench.ctrp<-affinClustering(dtgbench.ctrp,method="louvain")
# modularity(ceb.dtgbench.ctrp)
# 
# temp.result<-affinClustering(integrtStrctSensPert,method="spectral",K=max(ceb.dtgbench.ctrp$membership))

