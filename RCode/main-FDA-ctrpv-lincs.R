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
cDrugs <- preprocessInput(dname="ctrpv2", "lincs",venn_plot=TRUE)
dim(cDrugs)  ##239 X 28

# Process Sensitivity, Perturbation, and Structure layers for set of common drugs
sensData <- sensitivityData("ctrpv2", cDrugs)  ## 645 X 239
dim(sensData)
pertData <- perturbationData("lincs", cDrugs, "ctrpv2")  ## 978 X 239
dim(pertData)
strcData <- structureData("lincs", cDrugs)  ## a vector  --> 239 elemnts
length(strcData)


intersc <- cDrugs
df.drugs <-data.frame(MOLECULE_NAM=intersc$pert_iname,SMILES=intersc$canonical_smiles,STRUCTURE=intersc$structure_url)
write.csv(df.drugs,"Data/druginfo_ctrp.csv", row.names = TRUE,col.names = FALSE)


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

save(pertAffMat, file="Data/pertSimilarityMat-ctrpv2.RData")
save(sensAffMat, file="Data/sensSimilarityMat-ctrpv2.RData")
save(strcAffMat, file="Data/strcSimilarityMat-ctrpv2.RData")


integrtStrctSensPert <- integrateStrctSensPert_RNCE(sensAffMat,strcAffMat, pertAffMat,k1=128,k2=4,saverr="CTRP")

save(integrtStrctSensPert, file="Data/integrationSimilarityMat-ctrpv2.RData")



## Benchmarking and validation
## 1- DRUG-TARGET 
## loading and cleaning benchmark dataset
dataBench <- drugTargetBench("ctrpv",  commonDrugs) # 141 x 141 drug-drug adjacency matrix --> 141
load("Data/averageIorioPGX-all.RData") ## drug similarity matrix calculated based on Iorio snd Iskar et al. score (see iorioDIPS_PGX.R)
load("Data/averageIskarFinal.RData") ##load "iskar" results here ... (see iskar.R)
load("Data/drugERankSimil-CTRPV2-kendall.Rdata")
save(dataBench, file="Output/dataBench-ctrpv2.RData")

## Store the vairable to avoid the load to corrpute the file
storeintegrtStrctSensPert <- integrtStrctSensPert
integrtDNFLoad <- load("Data/ctrpv2-DNFIntegrated.RData")
integrtDNF <- get(integrtDNFLoad)
integrCISNFLoad <- load("Data/CISNF-ctrpv2.RData")


integrtCISNF <-get(integrCISNFLoad)
integrtStrctSensPert <- storeintegrtStrctSensPert

pairs <- generateDrugPairs(dataBench, strcAffMat, sensAffMat, pertAffMat, integrtStrctSensPert, average, finalIskarScore, NULL, drugERankSimil2,list(integrtDNF,integrtCISNF))




## compare cindices of combiantion layer vs. a single layer (e.g., structure)
res <- compConcordIndx(pairs)

cat("c.indexes values from each layer vs. the benchmark: \n integration: ", res$cindxLst$integrCindex, "\n structure: ", res$cindxLst$structureLayerCindex,
    "\n perturbation: ",  res$cindxLst$perturbationLayerCindex, "\n sensitivity: ", res$cindxLst$sensitivityLayerCindex, "\n Iorio: ", res$cindxLst$iorioCindex, 
    "\n Iskar: ", res$cindxLst$iskarCindex
    ,"\n DNF: ", res$cindxLst$dnf1Cindex, "\n CISNF:", res$cindxLst$dnf2Cindex)
cat("p-vals from the c.index comparison of integration layer vs. \n structure: ", res$pVals$intgrStrcPVal,"\n perturbation: ", res$pVals$intgrPertPVal,
    "\n sensitivity: ", res$pVals$intgrSensPVal, "\n Iorio: ", res$pVals$intgrIorioPVal, "\n Iskar: ", res$pVals$intgrIskarPVal
    ,"\n DNF: ", res$pVals$dnf1Pval, "\n CISNF: ", res$pVals$dnf2Pval)

## ROC and PR plots
cutoffCML<-generateRocPlot(pairs, d1Name="ctrpv2", d2Name="lincs", benchNam="drug-target",datestr=datestr,dnfnames=list("DNF","CISNF"))
generatePRPlot(pairs, d1Name="ctrpv2", d2Name="lincs", benchNam="drug-target",datestr=datestr,dnfnames=list("DNF","CISNF"))

## 2- ATC
dataBench3 <- ATCBench("chembl-new", cDrugs)
dim(dataBench3) ##[1]  51 51
##load "superPred" post-processed results here ... (see compareTo_suprPred.R)
load("Data/superPredSimil-ctrp-kendall.Rdata")
load("Data/averageIorioPGX-all.RData") ## drug similarity matrix calculated based on Iorio snd Iskar et al. score (see iorioDIPS_PGX.R)
load("Data/averageIskarFinal.RData") ##load "iskar" results here ... (see iskar.R)
save(dataBench3, file="Output/ATCBench-ctrpv2.RData")

pairs2 <- generateDrugPairs(dataBench3, strcAffMat, sensAffMat, pertAffMat, integrtStrctSensPert, average, finalIskarScore, superPredSimil2, NULL,list(integrtDNF,integrtCISNF))

## compare cindices of combination layer vs. a single layer (e.g., structure)
res2 <- compConcordIndx(pairs2)

cat("c.indexes values from each layer vs. the benchmark: \n integration: ", res2$cindxLst$integrCindex, "\n structure: ", res2$cindxLst$structureLayerCindex,
    "\n perturbation: ",  res2$cindxLst$perturbationLayerCindex, "\n sensitivity: ", res2$cindxLst$sensitivityLayerCindex, "\n Iorio: ", res2$cindxLst$iorioCindex, 
    "\n Iskar: ", res2$cindxLst$iskarCindex, "\n superPred: ", res2$cindxLst$superPredCindex
    ,"\n DNF: ", res$cindxLst$dnf1Cindex, "\n CISNF:", res2$cindxLst$dnf2Cindex)

cat("p-vals from the c.index comparison of integration layer vs. \n structure: ", res2$pVals$intgrStrcPVal,"\n perturbation: ", res2$pVals$intgrPertPVal,
    "\n sensitivity: ", res2$pVals$intgrSensPVal, "\n Iorio: ", res2$pVals$intgrIorioPVal, "\n Iskar: ", res2$pVals$intgrIskarPVal, "\n superPred: ", res2$pVals$intgrSuperPVal
    ,"\n DNF: ", res$pVals$dnf1Pval, "\n CISNF: ", res2$pVals$dnf2Pval)

## ROC and PR plots
cutoffATC<-generateRocPlot(pairs2, d1Name="ctrpv2", d2Name="lincs", benchNam="ATC(CHEMBL)",datestr=datestr,dnfnames=list("DNF","CISNF"))
generatePRPlot(pairs2, d1Name="ctrpv2", d2Name="lincs", benchNam="ATC(CHEMBL)",datestr=datestr,dnfnames=list("DNF","CISNF"))


## generate communities
load("Output/gmt_targ_ctrpv.RData")
communityGen(integrtStrctSensPert, "ctrpv2", GMT_TARG)
communityGen(integrtStrctSensPert, "ctrpv2", GMT_TARG,method = 'spectral',K=35)

## End(Not run)
# a more elaborate two-set Venn diagram with title and subtitle
venn.plot <- venn.diagram(
  x = list(
    "Intergative" = rownames(integrtStrctSensPert),
    "Drug target test" = rownames(dataBench),
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
pdf(file="Output/Venn_ctrpint_bench.pdf")
grid.draw(venn.plot)
dev.off()

emf(file="Output/Venn_ctrpint_bench.emf")
grid.draw(venn.plot)
dev.off()

pdf(file=paste("Output/",datestr,"/interg_ctrp_heatmap.pdf",sep = ""))
heatmap.2(1/(integrtStrctSensPert+1e-55),trace = "none",labRow = FALSE, labCol = FALSE)
dev.off()

pdf(file=paste("Output/",datestr,"/sens_ctrp_heatmap.pdf",sep = ""))
heatmap.2(1/(sensAffMat+1e-55),trace = "none",labRow = FALSE, labCol = FALSE)
dev.off()

pdf(file=paste("Output/",datestr,"/pert_ctrp_heatmap.pdf",sep = ""))
heatmap.2(1/(pertAffMat+1e-55),trace = "none",labRow = FALSE, labCol = FALSE)
dev.off()

pdf(file=paste("Output/",datestr,"/strc_ctrp_heatmap.pdf",sep = ""))
heatmap.2(1/(strcAffMat+1e-55),trace = "none",labRow = FALSE, labCol = FALSE)
dev.off()

pdf(file=paste("Output/",datestr,"/DNF_ctrp_heatmap.pdf",sep = ""))
heatmap.2(1/(integrtDNF+1e-55),trace = "none",labRow = FALSE, labCol = FALSE)
dev.off()
