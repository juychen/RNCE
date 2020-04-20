setwd("D:/rws/RNCE/")
library(netcom)  
library(apcluster)
library(xlsx)
library(igraph)
library(mclust)
library(netbiov)
library(devEMF)
library(ggplot2)
library(cluster)
library(dplyr)

library(rcdk)
library(ROCR)
library(survcomp)
library(reshape2)
library(proxy)
library(PRROC)


source('RCode/utilities.R')
source('RCode/affinityClustering.R')
source('Rcode/predPerf.R')
source('Rcode/generateDrugPairs.R')

# Deinfine unexpected charaters

badchars <- "[\xb5]|[\n]|[,]|[;]|[:]|[-]|[+]|[*]|[%]|[$]|[#]|[{]|[}]|[[]|[]]|[|]|[\\^]|[/]|[\\]|[.]|[_]|[ ]"
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

# Load data files
load("Output/dataBench-ctrpv2.RData") ##load "drug-target bench for ctrpv2" here ...
dtgbench.ctrp <- dataBench
load("Output/ATCBench-ctrpv2.RData") ##load "nci60" results here ... 
atcbench.ctrp <- dataBench3
load("Output/dataBench-nci60.RData") ##load "iskar" results here ... (see iskar.R)
dtgbench.nci60 <- dataBench1.5
load("Output/ATCBench-nci60.RData") ##load "iskar" results here ... (see iskar.R)
atcbench.nci60 <- dataBench3

# Load the result of the integration
name.integrate.ctrp<- "Data/integrationSimilarityMat-ctrpv2.RData"
name.integrate.nci60<- "Data/integrationSimilarityMat-nci60.RData"

# Load the result of the integration
name.integrate.ctrp.dnf<- "Data/ctrpv2-DNFIntegrated.RData"
name.integrate.nci60.dnf<- "Data/nci60-DNFIntegrated.RData"
name.integrate.ctrp.cisnf<- "Data/CISNF-ctrpv2.Rdata"
name.integrate.nci60.cisnf<- "Data/CISNF-nci60.RData"

# Network in each step
name.initial.nci60.rnce <- "Data/initSimilarityMat-nci60.Rdata"
name.initial.ctrp.rnce <- "Data/initSimilarityMat-CTRP.Rdata"
name.ci.ctrp.rnce<-"Data/rerkSimilarityMat-CTRP.Rdata"
name.ci.nci60.rnce<-"Data/rerkSimilarityMat-nci60.Rdata" 

name.save.ctrp<- gsub(badchars,"",name.integrate.ctrp)
name.save.nci60<- gsub(badchars,"",name.integrate.nci60)
name.save.ctrp.dnf<- gsub(badchars,"",name.integrate.ctrp.dnf)
name.save.nci60.dnf<- gsub(badchars,"",name.integrate.nci60.dnf)

# Load inputs of ctrpv2

load("Data/pertSimilarityMat-ctrpv2.RData")
pert.ctrp<-pertAffMat
load("Data/sensSimilarityMat-ctrpv2.RData")
sens.ctrp<-sensAffMat
load("Data/strcSimilarityMat-ctrpv2.RData")
strc.ctrp<-strcAffMat

# Load inputs of nci60
load("Data/pertSimilarityMat-nci60.RData")
pert.nci60<-pertAffMat
load("Data/sensSimilarityMat-nci60.RData")
sens.nci60<-sensAffMat
load("Data/strcSimilarityMat-nci60.RData")
strc.nci60<-strcAffMat

benchs<-c("dtg","atc")
datasources<-c("ctrp","nci60")
inputs<-c("strc","pert","sens")
methods<-c("","dnf","cisnf","in","ci")


# Funtion to assign names to the matrix
submatbyname<-function(mat,name){
  rmat <- mat[name,name]
  colnames(rmat)<-rownames(rmat) <- name
  return(rmat)
}  

# Get our data and comparing data
load(name.integrate.ctrp)
ctrp <- integrtStrctSensPert
diag(ctrp) <- 1

load(name.integrate.nci60)
nci60 <- integrtStrctSensPert
diag(nci60) <- 1

load(name.integrate.ctrp.dnf)
ctrp.dnf <- integrtStrctSensPert
diag(ctrp.dnf) <- 1

load(name.integrate.nci60.dnf)
nci60.dnf <- integrtStrctSensPert
diag(nci60.dnf) <- 1

load(name.integrate.ctrp.cisnf)
ctrp.cisnf <- integrtStrctSensPert
diag(ctrp.cisnf) <- 1

load(name.integrate.nci60.cisnf)
nci60.cisnf <- integrtStrctSensPert
diag(nci60.cisnf) <- 1

load(name.integrate.nci60.cisnf)
nci60.cisnf <- integrtStrctSensPert
diag(nci60.cisnf) <- 1

load(name.initial.ctrp.rnce)
ctrp.in<- integration
diag(ctrp.in) <- 1

load(name.initial.nci60.rnce)
nci60.in<- integration
diag(nci60.in) <- 1

load(name.ci.ctrp.rnce)
ctrp.ci<- re_ranked
diag(ctrp.ci) <- 1

load(name.ci.nci60.rnce)
nci60.ci<- re_ranked
diag(nci60.ci) <- 1

# process inputs diags
for(i in inputs){for(d in datasources){
  temp.sim<-get(paste(i,d,sep = "."))
  diag(temp.sim)<-1
  assign(paste(i,d,sep = "."),temp.sim)
}}


# Asssign names to the similarity matrix
name.ctrp.atc <- rownames(atcbench.ctrp)
name.ctrp.dtg <- rownames(dtgbench.ctrp)
name.nci60.atc <- rownames(atcbench.nci60)
name.nci60.dtg <- rownames(dtgbench.nci60)


# Select the benchmarking subsets ON CTRP AND NCI60

for(m in methods){for(d in datasources){for(b in benchs){
  
  if(m==""){temp.mat<-get(d)}
  else{temp.mat<-get(paste(d,m,sep = "."))}
  
  temp.name<-get(paste("name",d,b,sep="."))
  temp.result<-submatbyname(temp.mat,temp.name)
  assign(paste(d,".sub",m,".",b,sep=''),temp.result)
}}}

## Subset three Inputs 
# process inputs diags
for(i in inputs){for(d in datasources){for(b in benchs){
  temp.mat<-get(paste(i,d,sep = "."))
  temp.name<-get(paste("name",d,b,sep="."))
  temp.result<-submatbyname(temp.mat,temp.name)
  assign(paste(d,i,b,sep='.'),temp.result)
}}}


# Create the bench graph
graph.dtgbench.ctrp <- graph.adjacency(dtgbench.ctrp, mode="undirected",diag=FALSE)
graph.atcbench.ctrp <- graph.adjacency(atcbench.ctrp, mode="undirected",diag=FALSE)

ceb.dtgbench.ctrp<-cluster_louvain(graph.dtgbench.ctrp)
ceb.atcbench.ctrp<-cluster_louvain(graph.atcbench.ctrp)

graph.dtgbench.nci60 <- graph.adjacency(dtgbench.nci60, mode="undirected",diag=FALSE)
graph.atcbench.nci60 <- graph.adjacency(atcbench.nci60, mode="undirected",diag=FALSE)

ceb.dtgbench.nci60<-cluster_louvain(graph.dtgbench.nci60)
ceb.atcbench.nci60<-cluster_louvain(graph.atcbench.nci60)

# Apply the louvain on the benchmark ctrp dataset
loulabel.dtgbench.ctrp<-affinClustering(dtgbench.ctrp,method="louvain")
loulabel.atcbench.ctrp<-affinClustering(atcbench.ctrp,method="louvain")

# Apply the louvain on the benchmark nci60 dataset
loulabel.dtgbench.nci60<-affinClustering(dtgbench.nci60,method="louvain")
loulabel.atcbench.nci60<-affinClustering(atcbench.nci60,method="louvain")

# Apply the apcluster on the benchmark ATC dataset
apcomb.dtgbench.ctrp <- apcluster(dtgbench.ctrp, q=0.9)
apcomb.atcbench.ctrp <- apcluster(atcbench.ctrp, q=0.5)
apclabel.dtgbench.ctrp<-labels(apcomb.dtgbench.ctrp,"enum")
apclabel.atcbench.ctrp<-labels(apcomb.atcbench.ctrp,"enum")

# Apply the clustering each steps

## Spectral clustering for all tested networks 
for(m in methods){for(d in datasources){for(b in benchs){
  temp.mat<-get(paste(d,".sub",m,".",b,sep = ""))
  temp.label<-get(paste("loulabel.",b,"bench.",d,sep=""))
  temp.result<-affinClustering(temp.mat,method="spectral",K=max(temp.label))
  assign(paste("speclabel.",d,m,".",b,sep=''),temp.result)
}}}

# Three inputs spectral clustering labels
for(i in inputs){for(d in datasources){for(b in benchs){
  temp.mat<-get(paste(d,i,b,sep = "."))
  temp.label<-get(paste("loulabel.",b,"bench.",d,sep=""))
  temp.result<-affinClustering(temp.mat,method="spectral",K=max(temp.label))
  assign(paste("speclabel.",d,i,".",b,sep=''),temp.result)
}}}

# Socres

scores <- c()
measurement <- c()
benchdata <- c()
sensdata <- c()
conditions <- c()


# Modularities

# Modularity for the groond truth network
modularity(ceb.dtgbench.ctrp)
modularity(ceb.atcbench.ctrp)
modularity(ceb.dtgbench.nci60)
modularity(ceb.atcbench.nci60)

# Testing network modularities
for(m in methods){for(d in datasources){for(b in benchs){
  temp.graph <- get(paste("graph.",b,"bench.",d,sep=""))
  temp.label<-get(paste("speclabel.",d,m,".",b,sep=""))
  temp.result<-modularity(temp.graph,temp.label)
  assign(paste("modularity.",d,m,".",b,sep=''),temp.result)
  
  scores <- append(scores,temp.result)
  benchdata <- append(benchdata,b)
  sensdata <- append(sensdata,d)
  measurement <-append(measurement,'modularities')
  conditions<-append(conditions,m)
  
}}}


# Three inputs modluarities
for(i in inputs){for(d in datasources){for(b in benchs){
  temp.graph <- get(paste("graph.",b,"bench.",d,sep=""))
  temp.label<-get(paste("speclabel.",d,i,".",b,sep=""))
  temp.result<-modularity(temp.graph,method="spectral",temp.label)
  assign(paste("modularity.",d,i,".",b,sep=''),temp.result)
  
  
  scores <- append(scores,temp.result)
  benchdata <- append(benchdata,b)
  sensdata <- append(sensdata,d)
  measurement <-append(measurement,'modularities')
  conditions<-append(conditions,i)
}}}

# ARI

# Testing network modularities
for(m in methods){for(d in datasources){for(b in benchs){
  temp.spectral <- get(paste("speclabel.",d,m,".",b,sep=""))
  temp.label<-get(paste("loulabel.",b,"bench.",d,sep=""))
  temp.result<-adjustedRandIndex(temp.label,temp.spectral)
  assign(paste("ari.",d,m,".",b,sep=''),temp.result)
  
  
  scores <- append(scores,temp.result)
  benchdata <- append(benchdata,b)
  sensdata <- append(sensdata,d)
  measurement <-append(measurement,'ari')
  conditions<-append(conditions,m)
}}}


# Three inputs modluarities
for(i in inputs){for(d in datasources){for(b in benchs){
  temp.spectral <- get(paste("speclabel.",d,i,".",b,sep=""))
  temp.label<-get(paste("loulabel.",b,"bench.",d,sep=""))
  temp.result<-adjustedRandIndex(temp.label,temp.spectral)
  assign(paste("ari.",d,i,".",b,sep=''),temp.result)
  
  
  scores <- append(scores,temp.result)
  benchdata <- append(benchdata,b)
  sensdata <- append(sensdata,d)
  measurement <-append(measurement,'ari')
  conditions<-append(conditions,i)
}}}

# AUROC

for(m in methods){for(d in datasources){for(b in benchs){
  
  temp.bench <- get(paste(b,"bench.",d,sep=""))
  
  if(m==""){temp.mat<-get(d)}
  else{temp.mat <- get(paste(d,m,sep="."))}
  temp.pairs <- generateIntDrugPairs(temp.bench, temp.mat)
  temp.result<-predPerf(temp.pairs$integrPairs$obs.combiall, temp.pairs$benchPairs$bench)$auc
  assign(paste("auroc.",d,m,".",b,sep=''),temp.result)
  
  
  scores <- append(scores,temp.result)
  benchdata <- append(benchdata,b)
  sensdata <- append(sensdata,d)
  measurement <-append(measurement,'auroc')
  conditions<-append(conditions,m)
}}}


# Three inputs auroc
for(i in inputs){for(d in datasources){for(b in benchs){
  
  temp.bench <- get(paste(b,"bench.",d,sep=""))
  temp.mat <- get(paste(i,d,sep="."))
  temp.pairs <- generateIntDrugPairs(temp.bench, temp.mat)
  temp.result<-predPerf(temp.pairs$integrPairs$obs.combiall, temp.pairs$benchPairs$bench)$auc
  assign(paste("auroc.",d,m,".",b,sep=''),temp.result)
  
  
  scores <- append(scores,temp.result)
  benchdata <- append(benchdata,b)
  sensdata <- append(sensdata,d)
  measurement <-append(measurement,'auroc')
  conditions<-append(conditions,i)
}}}



# AUPRC

for(m in methods){for(d in datasources){for(b in benchs){
  
  temp.bench <- get(paste(b,"bench.",d,sep=""))
  if(m==""){temp.mat<-get(d)}
  else{temp.mat <- get(paste(d,m,sep="."))}
  temp.pairs <- generateIntDrugPairs(temp.bench, temp.mat)
  temp.result<-predPerf(temp.pairs$integrPairs$obs.combiall, temp.pairs$benchPairs$bench,plotType="PR")$auc.integral
  assign(paste("auroc.",d,m,".",b,sep=''),temp.result)
  
  
  scores <- append(scores,temp.result)
  print(temp.result)
  print(length(scores))
  benchdata <- append(benchdata,b)
  sensdata <- append(sensdata,d)
  measurement <-append(measurement,'auprc')
  conditions<-append(conditions,m)
}}}


# Three inputs auroc
for(i in inputs){for(d in datasources){for(b in benchs){
  
  temp.bench <- get(paste(b,"bench.",d,sep=""))
  temp.mat <- get(paste(i,d,sep="."))
  temp.pairs <- generateIntDrugPairs(temp.bench, temp.mat)
  temp.result<-predPerf(temp.pairs$integrPairs$obs.combiall, temp.pairs$benchPairs$bench,plotType="PR")$auc.integral
  assign(paste("auroc.",d,i,".",b,sep=''),temp.result)
  
  
  scores <- append(scores,temp.result)
  benchdata <- append(benchdata,b)
  sensdata <- append(sensdata,d)
  measurement <-append(measurement,'auprc')
  conditions<-append(conditions,i)
}}}

# scores <- c()
# measurement <- c()
# benchdata <- c()
# sensdata <- c()
# conditions <- c()

df.score <-data.frame(score=scores,metric=measurement,benchdata=benchdata,sensdata=sensdata,condition=conditions)