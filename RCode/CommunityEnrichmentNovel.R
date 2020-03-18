# Junyi Chen modified 
# Originated from Deena M.A. Gendoo
# Determining Community Enrichment for DNF communities, against Drug Targets and ATC classes. 

library(pheatmap)
library(xlsx)
source("RCode/communityGen.R")
library(devEMF)
library(svglite)
badchars <- "[\xb5]|[\n]|[,]|[;]|[:]|[-]|[+]|[*]|[%]|[$]|[#]|[{]|[}]|[[]|[]]|[|]|[\\^]|[/]|[\\]|[.]|[_]|[-]|[ ]"

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
# Load network inetegratrion

name.integrate.ctrp<- "Data/integrationSimilarityMat-ctrpv2.RData"
name.integrate.nc60<- "Data/integrationSimilarityMat-nci60.RData"

# Gen community For ctrp
load("Output/gmt_targ_ctrpv.RData")
load(name.integrate.ctrp)
ctrp <- integrtStrctSensPert
communityGen(ctrp, "ctrpv2", GMT_TARG,method = 'spectral',K=35)

# Gen communituy for nci60
load("./Output/gmt_targ_chembl.RData")
load(name.integrate.nc60)
nc60 <- integrtStrctSensPert
communityGen(nc60, "nci60", GMT_TARG,method = 'spectral',K=33)


# ### Load DNF communities
# A=load("CTRPv2comm.RData")
# CTRP_DNF = get(A)
# 
# B=load("NCI60comm.RData")
# NCI60_DNF = get(B)

### Load DNF communities
A=load("Data/CTRPv2commSpec.Rdata")
CTRP_DNF = get(A)

B=load("Data/NCI60commSpec.RData")
NCI60_DNF = get(B)

#Sanity Check
length(CTRP_DNF)  #should be 239
length(NCI60_DNF) #should be 238

### Load Benchmarks for Drug Targets
# CTRP - Internal benchmark (unchanged from previous iteration)
x=load("Output/gmt_targ_ctrpv.RData")
CTRP_GMT_TARG = get(x)

# NCI60 - Uniprot instead of Chembl
# Previously, y=load("./Output/gmt_targ_chembl.RData")
y=load("Output/gmt_targ_uniprot.RData")
NCI60_GMT_TARG= get(y)

### Load Benchmarks for ATC
# ATC - New Updated chembl as opposed to older chembl!
# NB: The GMT files (RData) for the ATC were externally renamed to reflect CTRP or NCI DNF

# Previously,z=load("./Output/gmt_atc_chembl_CTRP.RData")
z=load("Output/gmt_atc_chembl-new_CTRP.RData")
CTRP_GMT_ATC = get(z)

# Previously,k=load("./Output/gmt_atc_chembl_NCI60.RData")
k=load("Output/gmt_atc_chembl-new_NCI60.RData")
NCI60_GMT_ATC = get(k)

### THE FUNCTION!
communityEnrichment <- function(DNFcomm, DrugDNFName,BenchmarkName, GMT,threshold = 0.05) {

  HypergeomCalculations<-list()
  
  type <- typeof(DNFcomm)
  
  #Generate a list of drug communities from AP Clustering
  DrugClusters <- list()
  
  # If input is a apclustering result method
  if(type == "S4"){
    
    ClusterList <- list()
    
    for(i in 1:length(DNFcomm)){
      xx <- names(DNFcomm[[i]])
      DrugClusters[[i]] <- xx
      ClusterList[[i]] <- DNFcomm[[i]]
    }
    
    #ClusterList <- DrugClusters
  }else{ # If input community is a named list of labels.
    
    ClusterList <- list()
    
    for(i in 1:max(DNFcomm)){
      xx <- names(DNFcomm[DNFcomm==i])
      DrugClusters[[i]] <- xx
      ClusterList[[i]] <- DNFcomm[DNFcomm==i]
    }
    
  }

  #Add community numbers to the AP clusters
  names(DrugClusters)<-paste("C",rep(1:length(DrugClusters)),sep="")
  #Number of drugs per community
  indx <- sapply(DrugClusters, length)

  #Check if the genesets from the Drug Clusters intersect with the drugs in the GMT
  # If not, filter out the non-intersecting drug communities
  #ClusterList <-DNFcomm@clusters
  #names(ClusterList)<-paste("C",rep(1:length(DNFcomm@clusters)),sep="")
  names(ClusterList)<-paste("C",rep(1:length(ClusterList)),sep="")
  # Comment this line since all drug are selcted in our test
  #Bkrd_Drugs<-length((unique(rownames(DNFcomm))))
  Bkrd_Drugs<- 0
  
  Bkrd_Benchmark<- length(unique(c(as.character(unlist(GMT))))) #All the drugs in the GMT benchmark, N=141
  Bkrd_Benchmark_Drugs<- unique(c(as.character(unlist(GMT)))) #All the drugs in the GMT benchmark, N=141  
  
  # Reduce the community to keep only the drugs, in each community, which intersect with the Bkrd Benchmark
    ClusterList <- lapply(ClusterList,function(x){if (length(intersect(Bkrd_Benchmark_Drugs,names(x[]))) >= 1) {intersect(Bkrd_Benchmark_Drugs,names(x[]))} else {NULL}})
    ClusterByCondition<- sapply(ClusterList, function(x) length(x) > 1) 

    EnrichmentMatrix <- matrix(NA, nrow=length(DrugClusters), ncol=length(GMT), 
                             dimnames=list(names(DrugClusters), names(GMT)))
  
  ###### hypergeometric
  for (i in 1:length(DrugClusters)) { #For every drug cluster/community
    for (j in 1:length(GMT)) { # For every geneset of drugs associated with a particular target
      
      DrugClus <- unique(unlist(DrugClusters[[i]])) #drugs per population
      DrugTarg <-  unique(as.character(unlist(GMT[[j]]))) #drugs per geneset
      
      LenClus <- length(DrugClus)
      LenTarg <- length(DrugTarg)
      
      ClusAndTarg <- length(intersect(DrugClus, DrugTarg)) #intersected drugs between the two sets (+ves)
      ClusAndNotTarg<-length(DrugClus)-ClusAndTarg #exist in drug cluster, but not found in target geneset
      TargAndNotClus<-length(DrugTarg)-ClusAndTarg #exist in target geneset, but not found in drug cluster
      NotTargAndNotClus<-Bkrd_Benchmark-(ClusAndTarg+ClusAndNotTarg+TargAndNotClus)

      #Conduct the test only if there is at least 1 drug overlap between a community and drug target set
      if(ClusAndTarg>=1){
        minxx3 <- fisher.test(matrix(c(ClusAndTarg,ClusAndNotTarg,TargAndNotClus,NotTargAndNotClus),
                                     nrow=2,ncol=2),alternative="greater")$p.value
        EnrichmentMatrix[i, j] <- minxx3 }
    }
  }
  
  # Replace the 'NA' hits (ie, non-overlaps) with p-value of 1
  EnrichmentMatrix[is.na(EnrichmentMatrix)]<-1
  # Do a correction for multiple testing
  EnrichmentMatrix_FDRcorrected<-apply(EnrichmentMatrix,1,p.adjust,method="fdr")
  write.xlsx(EnrichmentMatrix_FDRcorrected,file=paste(subDir,"/","EnrichmentMatrix_FDRcorrected_",DrugDNFName,"_",BenchmarkName,".xls",sep=""))
  
  # Find out how many communities have at least one calculated p-value
  # (ie, how many communities have calculations on thems)
  EnrichmentMatrix_FDRcorrected<-EnrichmentMatrix_FDRcorrected[!apply(EnrichmentMatrix_FDRcorrected,1,function(x){all(x>=threshold)}),]
  EnrichmentMatrix_FDRcorrected<-EnrichmentMatrix_FDRcorrected[,!apply(EnrichmentMatrix_FDRcorrected,2,function(x){all(x>=threshold)})]
  dim(EnrichmentMatrix_FDRcorrected) #rows=targets, columns = communities
  
  #Final Cleanup for publication
  gsub(rownames(EnrichmentMatrix_FDRcorrected),pattern = "Target_",replacement = "")
  gsub(rownames(EnrichmentMatrix_FDRcorrected),pattern = "ATC_",replacement = "")

  # Use only a few color codes: FDR < 0.05; FDR < 0.01; FDR < 0.001
  # FYI: -log10(0.01)=2 and -log10(0.1)=1 . Ergo, higher number is more significant here!
  svglite(paste(subDir,"/","Community_Enrichment_FDR_5e-2_",DrugDNFName,"_",BenchmarkName,".svg",sep = ""),width = 5,height = 5.5)
  pheatmap(-log10((EnrichmentMatrix_FDRcorrected)),cluster_rows = F,cluster_cols = F,
           color = c('#E6E6FA',rev(c('#9e0142','#d53e4f','#f46d43','#fdae61','#fee08b','#ffffbf','#e6f598','#abdda4','#66c2a5','#3288bd','#5e4fa2')))) #higher number = more significant
  dev.off()
  
  svglite(paste(subDir,"/","Community_Enrichment_",DrugDNFName,"_",BenchmarkName,".svg",sep = ""),width = 12,height = 10)
  pheatmap(-log10(t(EnrichmentMatrix)),cluster_rows = F,cluster_cols = F,
           color = c('#E6E6FA',rev(c('#9e0142','#d53e4f','#f46d43','#fdae61','#fee08b','#ffffbf','#e6f598','#abdda4','#66c2a5','#3288bd','#5e4fa2')))) #higher number = more significant
  dev.off()
} 

### GENERATE PLOTS

#CTRPv2 with Internal Drug Target Benchmark
communityEnrichment(DNFcomm = CTRP_DNF,GMT = CTRP_GMT_TARG,DrugDNFName = "CTRP",BenchmarkName = "Target")
#NCI60 with Drug Target Benchmark (Chembl)
communityEnrichment(DNFcomm = NCI60_DNF,GMT = NCI60_GMT_TARG,DrugDNFName = "NCI60",BenchmarkName = "Target")

# CTRP with ATC Benchmark (Chembl)
communityEnrichment(DNFcomm = CTRP_DNF,GMT = CTRP_GMT_ATC,DrugDNFName = "CTRP",BenchmarkName = "ATC")
# NCI60 with ATC Benchmark (Chembl)
communityEnrichment(DNFcomm = NCI60_DNF,GMT = NCI60_GMT_ATC,DrugDNFName = "NCI60",BenchmarkName = "ATC")
