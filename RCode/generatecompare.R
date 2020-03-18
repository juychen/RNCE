## The code assumes that the working directory is "RNCE" folder
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
pertAffMat <- constPerturbationLayerDNF(pertData)


integrtStrctSensPert <- integrateStrctSensPert_cisnf(sensAffMat,strcAffMat, pertAffMat,k1=192,k2=4)
save(integrtStrctSensPert, file="Data/CISNF-ctrpv2.RData")


integrtStrctSensPert <- integrateStrctSensPertR(sensAffMat,strcAffMat, pertAffMat)
save(integrtStrctSensPert, file="Data/RSNF-ctrpv2.RData")

integrtStrctSensPert <- integrateStrctSensPert(sensAffMat,strcAffMat, pertAffMat)
save(integrtStrctSensPert, file="Data/ctrpv2-DNFIntegrated.RData")


pertAffMat <- constPerturbationLayer(pertData)
strcAffMat <- (normalize(strcAffMat))
sensAffMat <- (normalize(sensAffMat))
pertAffMat <- (normalize(pertAffMat))

integrtStrctSensPert <- integrateStrctSensPert_RSNF_RNCI(sensAffMat,strcAffMat, pertAffMat)
save(integrtStrctSensPert, file="Data/RNCI-ctrpv2.RData")



##########################################################
#NCI60

cDrugs <- preprocessInput(dname="nci60", "lincs")
dim(cDrugs$lincsboth)  ##  238 X 28
dim(cDrugs$nciboth)  ##  238 X 66


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
pertAffMat <- constPerturbationLayerDNF(pertData)


integrtStrctSensPert <- integrateStrctSensPert_cisnf(sensAffMat,strcAffMat, pertAffMat,k1=192,k2=4)

save(integrtStrctSensPert, file="Data/CISNF-nci60.RData")

integrtStrctSensPert <- integrateStrctSensPertR(sensAffMat,strcAffMat, pertAffMat)
save(integrtStrctSensPert, file="Data/RSNF-nci60.RData")

integrtStrctSensPert <- integrateStrctSensPert(sensAffMat,strcAffMat, pertAffMat)
save(integrtStrctSensPert, file="Data/nci60-DNFIntegrated.RData")


pertAffMat <- constPerturbationLayer(pertData,K = 20)
strcAffMat <- (normalize(strcAffMat))
sensAffMat <- (normalize(sensAffMat))
pertAffMat <- (normalize(pertAffMat))


integrtStrctSensPert <- integrateStrctSensPert_RSNF_RNCI(sensAffMat,strcAffMat, pertAffMat)
save(integrtStrctSensPert, file="Data/RNCI-nci60.RData")

