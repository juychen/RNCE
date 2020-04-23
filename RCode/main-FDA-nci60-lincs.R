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

intersc <- cDrugs$lincsboth
df.drugs <-data.frame(MOLECULE_NAM=intersc$pert_iname,SMILES=intersc$canonical_smiles,STRUCTURE=intersc$structure_url)
write.csv(df.drugs,"Data/druginfo_nci60.csv", row.names = TRUE,col.names = FALSE)


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


save(pertAffMat, file="Data/pertSimilarityMat-nci60.RData")
save(sensAffMat, file="Data/sensSimilarityMat-nci60.RData")
save(strcAffMat, file="Data/strcSimilarityMat-nci60.RData")


integrtStrctSensPert <- integrateStrctSensPert_RNCE(sensAffMat,strcAffMat, pertAffMat,k1=64,k2=4,saverr="nci60")

save(integrtStrctSensPert, file="Data/integrationSimilarityMat-nci60.RData")



# Sanity Check - should all have the same dimensions: 237 X 237 --> 238...
dim(strcAffMat)
dim(sensAffMat)
dim(pertAffMat)
dim(integrtStrctSensPert)



## Benchmarking and validation
## 1- DRUG-TARGET 
dataBench1.5 <- drugTargetBench("uniprot", commonDrugs)
dim(dataBench1.5) ##[1] 86
load("Data/averageIorioPGX-all.RData") ## drug similarity matrix calculated based on Iorio snd Iskar et al. score (see iorioDIPS_PGX.R)
load("Data/averageIskarFinal.RData") ##load "iskar" results here ... (see iskar.R)
load("Data/drugERankSimil-nci60-kendall.Rdata")
save(dataBench1.5, file="Output/dataBench-nci60.RData")


## Store the vairable to avoid the load to corrpute the file
storeintegrtStrctSensPert <- integrtStrctSensPert
integrtDNFLoad <- load("Data/nci60-DNFIntegrated.RData")
integrtDNF <- get(integrtDNFLoad)
integrCISNFLoad <- load("Data/CISNF-nci60.RData")
integrtCISNF <-get(integrCISNFLoad)
integrtStrctSensPert <- storeintegrtStrctSensPert

pairs1 <- generateDrugPairs(dataBench1.5, strcAffMat, sensAffMat, pertAffMat, integrtStrctSensPert, average, finalIskarScore, NULL, drugERankSimil2,list(integrtDNF,integrtCISNF))

## compare cindices of integration layer vs. all single layers (e.g., structure, perturbation, sensitivity )
res1 <- compConcordIndx(pairs1)

cat("c.indexes values from each layer vs. the benchmark: \n integration: ", res1$cindxLst$integrCindex, "\n structure: ", res1$cindxLst$structureLayerCindex,
    "\n perturbation: ",  res1$cindxLst$perturbationLayerCindex, "\n sensitivity: ", res1$cindxLst$sensitivityLayerCindex, "\n Iorio: ", res1$cindxLst$iorioCindex, 
    "\n Iskar: ", res1$cindxLst$iskarCindex
    ,"\n DNF: ", res1$cindxLst$dnf1Cindex, "\n CISNF:", res1$cindxLst$dnf2Cindex)

cat("p-vals from the c.index comparison of integration layer vs. \n structure: ", res1$pVals$intgrStrcPVal,"\n perturbation: ", res1$pVals$intgrPertPVal,
    "\n sensitivity: ", res1$pVals$intgrSensPVal, "\n Iorio: ", res1$pVals$intgrIorioPVal, "\n Iskar: ", res1$pVals$intgrIskarPVal
    ,"\n DNF: ", res1$pVals$dnf1Pval, "\n CISNF: ", res1$pVals$dnf2Pval)

## ROC and PR plots
cutoffCHEMBL<-generateRocPlot(pairs1, d1Name="nci60", d2Name="lincs", benchNam="drug-target(UNIPROT)",datestr=datestr,list("DNF","CISNF"))
generatePRPlot(pairs1, d1Name="nci60", d2Name="lincs", benchNam="drug-target(UNIPROT)",datestr=datestr,list("DNF","CISNF"))


## 2- ATC
dataBench3 <- ATCBench("chembl-new", cDrugs$lincsboth,benchname="NCI60")
dim(dataBench3) ##[1] 72 72
##load "superPred" post-processed results here ... (see compareTo_suprPred.R)
load("Data/superPredSimil-nci60-kendall.Rdata")
load("Data/averageIorioPGX-all.RData") ## drug similarity matrix calculated based on Iorio snd Iskar et al. score (see iorioDIPS_PGX.R)
##load "iskar" results here ... (see iskar.R)
load("Data/averageIskarFinal.RData")
save(dataBench3, file="Output/ATCBench-nci60.RData")

pairs2 <- generateDrugPairs(dataBench3, strcAffMat, sensAffMat, pertAffMat, integrtStrctSensPert, average, finalIskarScore, superPredSimil2, NULL,list(integrtDNF,integrtCISNF))

## compare cindices of combination layer vs. a single layer (e.g., structure)
res2 <- compConcordIndx(pairs2)

cat("c.indexes values from each layer vs. the benchmark: \n integration: ", res2$cindxLst$integrCindex, "\n structure: ", res2$cindxLst$structureLayerCindex,
    "\n perturbation: ",  res2$cindxLst$perturbationLayerCindex, "\n sensitivity: ", res2$cindxLst$sensitivityLayerCindex, "\n Iorio: ", res2$cindxLst$iorioCindex, 
    "\n Iskar: ", res2$cindxLst$iskarCindex, "\n superPred: ", res2$cindxLst$superPredCindex
    ,"\n DNF: ", res2$cindxLst$dnf1Cindex, "\n CISNF:", res2$cindxLst$dnf2Cindex)

    
cat("p-vals from the c.index comparison of integration layer vs. \n structure: ", res2$pVals$intgrStrcPVal,"\n perturbation: ", res2$pVals$intgrPertPVal,
    "\n sensitivity: ", res2$pVals$intgrSensPVal, "\n Iorio: ", res2$pVals$intgrIorioPVal, "\n Iskar: ", res2$pVals$intgrIskarPVal, "\n superPred: ", res2$pVals$intgrSuperPVal
    ,"\n DNF: ", res2$pVals$dnf1Pval, "\n CISNF: ", res2$pVals$dnf2Pval)

## ROC and PR plots
cutoffATC<-generateRocPlot(pairs2, d1Name="nci60", d2Name="lincs", benchNam="ATC(CHEMBL)-Zscore",datestr=datestr,list("DNF","CISNF"))
generatePRPlot(pairs2, d1Name="nci60", d2Name="lincs", benchNam="ATC(CHEMBL)-Zscore",datestr=datestr,list("DNF","CISNF"))


## generate communities

gmt_targ_temp <- drugTargetBench("chembl", commonDrugs)
load("./Output/gmt_targ_chembl.RData")
#communityGen(integrtStrctSensPert, "nci60", GMT_TARG)
communityGen(integrtStrctSensPert, "nci60", GMT_TARG,method = 'spectral',K=33)




## End(Not run)
# a more elaborate two-set Venn diagram with title and subtitle
venn.plot <- venn.diagram(
  x = list(
    "Intergative" = rownames(integrtStrctSensPert),
    "Drug target test" = rownames(dataBench1.5),
    "ATC test" = rownames(dataBench3)
  ),
  filename = NULL, 
  scaled = FALSE,
  fill = c("white", "purple","green"),
  cat.fontfamily = rep("sans", 3),
  cat.dist = c(0.05,0.02,0.02),
  cat.default.pos ="text",  
  #main = "Complex Venn Diagram",
  #main.fontfamily ="sans",
  #sub = "Featuring: rotation and external lines",
  sub.fontfamily ="sans",
  
  main.cex = 2,
  sub.cex = 1,
  fontfamily =  "sans"
  
)
pdf(file="Output/Venn_nci60int_bench.pdf")
grid.draw(venn.plot)
dev.off()


emf(file="Output/Venn_nci60int_bench.emf")
grid.draw(venn.plot)
dev.off()


pdf(file=paste("Output/",datestr,"/interg_nci60_heatmap.pdf",sep = ""))
heatmap.2(1/(integrtStrctSensPert+1e-55),trace = "none",labRow = FALSE, labCol = FALSE)
dev.off()

pdf(file=paste("Output/",datestr,"/sens_nci60_heatmap.pdf",sep = ""))
heatmap.2(1/(sensAffMat+1e-55),trace = "none",labRow = FALSE, labCol = FALSE)
dev.off()

pdf(file=paste("Output/",datestr,"/pert_nci60_heatmap.pdf",sep = ""))
heatmap.2(1/(pertAffMat+1e-55),trace = "none",labRow = FALSE, labCol = FALSE)
dev.off()

pdf(file=paste("Output/",datestr,"/strc_nci60_heatmap.pdf",sep = ""))
heatmap.2(1/(strcAffMat+1e-55),trace = "none",labRow = FALSE, labCol = FALSE)
dev.off()

pdf(file=paste("Output/",datestr,"/DNF_nci60_heatmap.pdf",sep = ""))
heatmap.2(1/(integrtDNF+1e-55),trace = "none",labRow = FALSE, labCol = FALSE)
dev.off()


