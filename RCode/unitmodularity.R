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
source('RCode/utilities.R')
source('RCode/affinityClustering.R')

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
ATCbench.ctrp <- dataBench3
load("Output/dataBench-nci60.RData") ##load "iskar" results here ... (see iskar.R)
dtgbench.nci60 <- dataBench1.5
load("Output/ATCBench-nci60.RData") ##load "iskar" results here ... (see iskar.R)
ATCbench.nci60 <- dataBench3

name.integrate.ctrp<- "Data/integrationSimilarityMat-ctrpv2.RData"
name.integrate.nc60<- "Data/integrationSimilarityMat-nci60.RData"

name.integrate.ctrp.dnf<- "Data/ctrpv2-DNFIntegrated.RData"
name.integrate.nc60.dnf<- "Data/nci60-DNFIntegrated.RData"


name.integrate.ctrp.cisnf<- "Data/CISNF-ctrpv2.Rdata"
name.integrate.nc60.cisnf<- "Data/CISNF-nci60.RData"

# Network in each step
name.initial.nci60.rnce <- "Data/initSimilarityMat-nci60.Rdata"
name.initial.ctrp.rnce <- "Data/initSimilarityMat-CTRP.Rdata"
name.ci.ctrp.rnce<-"Data/rerkSimilarityMat-CTRP.Rdata"
name.ci.nci60.rnce<-"Data/rerkSimilarityMat-nci60.Rdata" 

name.save.ctrp<- gsub(badchars,"",name.integrate.ctrp)
name.save.nc60<- gsub(badchars,"",name.integrate.nc60)

name.save.ctrp.dnf<- gsub(badchars,"",name.integrate.ctrp.dnf)
name.save.nc60.dnf<- gsub(badchars,"",name.integrate.nc60.dnf)

#load("Data/integrationSimilarityMat-nci60.RData") ##load "iskar" results here ... (see iskar.R)\
#load("Data/nci60-DNFIntegrated.RData") ##load "iskar" results here ... (see iskar.R)

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

load(name.integrate.nc60)
nc60 <- integrtStrctSensPert
diag(nc60) <- 1

load(name.integrate.ctrp.dnf)
ctrp.dnf <- integrtStrctSensPert
diag(ctrp.dnf) <- 1

load(name.integrate.nc60.dnf)
nc60.dnf <- integrtStrctSensPert
diag(nc60.dnf) <- 1

load(name.integrate.ctrp.cisnf)
ctrp.cisnf <- integrtStrctSensPert
diag(ctrp.cisnf) <- 1

load(name.integrate.nc60.cisnf)
nc60.cisnf <- integrtStrctSensPert
diag(nc60.cisnf) <- 1

load(name.integrate.nc60.cisnf)
nc60.cisnf <- integrtStrctSensPert
diag(nc60.cisnf) <- 1

load(name.initial.ctrp.rnce)
ctrp.inrnce<- integration
diag(ctrp.inrnce) <- 1

load(name.initial.nci60.rnce)
nci60.inrnce<- integration
diag(nci60.inrnce) <- 1

load(name.ci.ctrp.rnce)
ctrp.cirnce<- re_ranked
diag(ctrp.cirnce) <- 1

load(name.ci.nci60.rnce)
nci60.cirnce<- re_ranked
diag(nci60.cirnce) <- 1

# Asssign names to the similarity matrix
name.ctrp.atc <- rownames(ATCbench.ctrp)
name.ctrp.dtg <- rownames(dtgbench.ctrp)
name.nci60.atc <- rownames(ATCbench.nci60)
name.nci60.dtg <- rownames(dtgbench.nci60)



# Select the benchmarking subsets ON CTRP AND NCI60
ctrp.sub.dtg <- submatbyname(ctrp,name.ctrp.dtg)
ctrp.sub.atc <- submatbyname(ctrp,name.ctrp.atc)

nc60.sub.dtg <- submatbyname(nc60,name.nci60.dtg)
nc60.sub.atc <- submatbyname(nc60,name.nci60.atc)

ctrp.subdnf.dtg <- submatbyname(ctrp.dnf,name.ctrp.dtg)
ctrp.subdnf.atc <- submatbyname(ctrp.dnf,name.ctrp.atc)

nc60.subdnf.dtg <- submatbyname(nc60.dnf,name.nci60.dtg)
nc60.subdnf.atc <- submatbyname(nc60.dnf,name.nci60.atc)

ctrp.subcisnf.dtg <- submatbyname(ctrp.cisnf,name.ctrp.dtg)
ctrp.subcisnf.atc <- submatbyname(ctrp.cisnf,name.ctrp.atc)

nc60.subcisnf.dtg <- submatbyname(nc60.cisnf,name.nci60.dtg)
nc60.subcisnf.atc <- submatbyname(nc60.cisnf,name.nci60.atc)


ctrp.subci.dtg <- submatbyname(ctrp.cirnce,name.ctrp.dtg)
ctrp.subci.atc <- submatbyname(ctrp.cirnce,name.ctrp.atc)

nc60.subci.dtg <- submatbyname(nci60.cirnce,name.nci60.dtg)
nc60.subci.atc <- submatbyname(nci60.cirnce,name.nci60.atc)


ctrp.subin.dtg <- submatbyname(ctrp.inrnce,name.ctrp.dtg)
ctrp.subin.atc <- submatbyname(ctrp.inrnce,name.ctrp.atc)

nc60.subin.dtg <- submatbyname(nci60.inrnce,name.nci60.dtg)
nc60.subin.atc <- submatbyname(nci60.inrnce,name.nci60.atc)


# Create the bench graph
graph.dtgbench.ctrp <- graph.adjacency(dtgbench.ctrp, mode="undirected",diag=FALSE)
graph.actbench.ctrp <- graph.adjacency(ATCbench.ctrp, mode="undirected",diag=FALSE)


ceb.dtgbench.ctrp<-cluster_louvain(graph.dtgbench.ctrp)
ceb.actbench.ctrp<-cluster_louvain(graph.actbench.ctrp)


graph.dtgbench.nci60 <- graph.adjacency(dtgbench.nci60, mode="undirected",diag=FALSE)
graph.actbench.nci60 <- graph.adjacency(ATCbench.nci60, mode="undirected",diag=FALSE)

ceb.dtgbench.nci60<-cluster_louvain(graph.dtgbench.nci60)
ceb.actbench.nci60<-cluster_louvain(graph.actbench.nci60)

# Apply the louvain on the benchmark ctrp dataset
loulabel.dtgbench.ctrp<-affinClustering(dtgbench.ctrp,method="louvain")
loulabel.atcbench.ctrp<-affinClustering(ATCbench.ctrp,method="louvain")


# Apply the louvain on the benchmark nci60 dataset
loulabel.dtgbench.nci60<-affinClustering(dtgbench.nci60,method="louvain")
loulabel.atcbench.nci60<-affinClustering(ATCbench.nci60,method="louvain")

# Apply the apcluster on the benchmark ATC dataset
apcomb.dtgbench.ctrp <- apcluster(dtgbench.ctrp, q=0.9)
apcomb.ATCbench.ctrp <- apcluster(ATCbench.ctrp, q=0.5)
apclabel.dtgbench.ctrp<-labels(apcomb.dtgbench.ctrp,"enum")
apclabel.ATCbench.ctrp<-labels(apcomb.ATCbench.ctrp,"enum")

# Apply the clustering each steps
## One the initial step:

speclabel.ctrpin.dtg<-affinClustering(ctrp.subin.dtg,method="spectral",K=max(loulabel.dtgbench.ctrp))
speclabel.nci60in.dtg<-affinClustering(nc60.subin.dtg,method="spectral",K=max(loulabel.dtgbench.nci60))

speclabel.ctrpin.atc<-affinClustering(ctrp.subin.atc,method="spectral",K=max(loulabel.atcbench.ctrp))
speclabel.nci60in.atc<-affinClustering(nc60.subin.atc,method="spectral",K=max(loulabel.atcbench.nci60))



## On the ci step:
speclabel.ctrpci.dtg<-affinClustering(ctrp.subci.dtg,method="spectral",K=max(loulabel.dtgbench.ctrp))
speclabel.nci60ci.dtg<-affinClustering(nc60.subci.dtg,method="spectral",K=max(loulabel.dtgbench.nci60))

speclabel.ctrpci.atc<-affinClustering(ctrp.subci.atc,method="spectral",K=max(loulabel.atcbench.ctrp))
speclabel.nci60ci.atc<-affinClustering(nc60.subci.dtg,method="spectral",K=max(loulabel.atcbench.nci60))

modularity(graph.dtgbench.ctrp, speclabel.ctrpci.dtg)
modularity(graph.dtgbench.nci60, speclabel.nci60ci.dtg)
modularity(graph.actbench.ctrp, speclabel.ctrpci.atc)
modularity(graph.actbench.nci60, speclabel.nci60ci.atc)


