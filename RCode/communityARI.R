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

# Create the bench graph
graph.dtgbench.ctrp <- graph.adjacency(dtgbench.ctrp, mode="undirected",diag=FALSE)
graph.actbench.ctrp <- graph.adjacency(ATCbench.ctrp, mode="undirected",diag=FALSE)


ceb.dtgbench.ctrp<-cluster_louvain(graph.dtgbench.ctrp)
ceb.actbench.ctrp<-cluster_louvain(graph.actbench.ctrp)

emf(file="Output/Graph_cheMbl_bench_ctrp.emf")
plot(ceb.dtgbench.ctrp,graph.dtgbench.ctrp,vertex.size=3,edge.arrow.size=0, edge.curved=0,edge.size=0,
     vertex.label=NA
     
     #vertex.color="orange", vertex.frame.color="#555555"#, #vertex.label.color="black",
     
     #,vertex.label = c()
     #,vertex.label.cex=0
     )
dev.off()


emf(file="Output/Graph_atc_bench_ctrp.emf")

plot(ceb.actbench.ctrp,graph.actbench.ctrp,vertex.size=3,edge.arrow.size=0, edge.curved=0,edge.size=0,
     vertex.label=NA
     
     #vertex.color="orange", vertex.frame.color="#555555"#, vertex.label.color="black",
     
     #vertex.label.cex=.7)
     #vertex.label.cex=0
     )
dev.off()

graph.dtgbench.nci60 <- graph.adjacency(dtgbench.nci60, mode="undirected",diag=FALSE)
graph.actbench.nci60 <- graph.adjacency(ATCbench.nci60, mode="undirected",diag=FALSE)

ceb.dtgbench.nci60<-cluster_louvain(graph.dtgbench.nci60)
ceb.actbench.nci60<-cluster_louvain(graph.actbench.nci60)


emf(file="Output/Graph_cheMbl_bench_nci60.emf")
plot(ceb.dtgbench.nci60,graph.dtgbench.nci60,vertex.size=3,edge.arrow.size=0, edge.curved=0,edge.size=0,
     vertex.label=NA)
   
dev.off()

emf(file="Output/Graph_atc_bench_nci60.emf")

plot(ceb.actbench.nci60,graph.actbench.nci60,vertex.size=3,edge.arrow.size=0, edge.curved=0,edge.size=0,
     vertex.label=NA)
     
dev.off()


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

# 10 RUNS my and dnf

folds = 10

# label <-ceb.dtgbench.ctrp$membership
# names(label) <- names(loulabel.dtgbench.ctrp)
# spec10.ctrp1.dtg <- nRunClusterPerformance(ctrp.sub.dtg,label,folds=folds ,method = "spectral",K=max(label))
# spec10.ctrpdnf1.dtg <- nRunClusterPerformance(ctrp.subdnf.dtg,label,folds=folds ,method = "spectral",K=max(label))


spec10.ctrp.dtg <- nRunClusterPerformance(ctrp.sub.dtg,loulabel.dtgbench.ctrp,folds=folds ,method = "spectral",K=max(loulabel.dtgbench.ctrp))
spec10.ctrp.atc <- nRunClusterPerformance(ctrp.sub.atc,loulabel.atcbench.ctrp,folds=folds ,method = "spectral",K=max(loulabel.atcbench.ctrp))
spec10.nc60.dtg <- nRunClusterPerformance(nc60.sub.dtg,loulabel.dtgbench.nci60,folds=folds ,method = "spectral",K=max(loulabel.dtgbench.nci60))
spec10.nc60.atc <- nRunClusterPerformance(nc60.sub.atc,loulabel.atcbench.nci60,folds=folds ,method = "spectral",K=max(loulabel.atcbench.nci60))

spec10.ctrpdnf.dtg <- nRunClusterPerformance(ctrp.subdnf.dtg,loulabel.dtgbench.ctrp,folds=folds ,method = "spectral",K=max(loulabel.dtgbench.ctrp))
spec10.ctrpdnf.atc <- nRunClusterPerformance(ctrp.subdnf.atc,loulabel.atcbench.ctrp,folds=folds ,method = "spectral",K=max(loulabel.atcbench.ctrp))
spec10.nc60dnf.dtg <- nRunClusterPerformance(nc60.subdnf.dtg,loulabel.dtgbench.nci60,folds=folds ,method = "spectral",K=max(loulabel.dtgbench.nci60))
spec10.nc60dnf.atc <- nRunClusterPerformance(nc60.subdnf.atc,loulabel.atcbench.nci60,folds=folds ,method = "spectral",K=max(loulabel.atcbench.nci60))

spec10.ctrpcisnf.dtg <- nRunClusterPerformance(ctrp.subcisnf.dtg,loulabel.dtgbench.ctrp,folds=folds ,method = "spectral",K=max(loulabel.dtgbench.ctrp))
spec10.ctrpcisnf.atc <- nRunClusterPerformance(ctrp.subcisnf.atc,loulabel.atcbench.ctrp,folds=folds ,method = "spectral",K=max(loulabel.atcbench.ctrp))
spec10.nc60cisnf.dtg <- nRunClusterPerformance(nc60.subcisnf.dtg,loulabel.dtgbench.nci60,folds=folds ,method = "spectral",K=max(loulabel.dtgbench.nci60))
spec10.nc60cisnf.atc <- nRunClusterPerformance(nc60.subcisnf.atc,loulabel.atcbench.nci60,folds=folds ,method = "spectral",K=max(loulabel.atcbench.nci60))

print(wilcox.test(spec10.ctrp.dtg,spec10.ctrpdnf.dtg))
print(wilcox.test(spec10.ctrp.atc,spec10.ctrpdnf.atc))
print(wilcox.test(spec10.nc60.dtg,spec10.nc60dnf.dtg))
print(wilcox.test(spec10.nc60.atc,spec10.nc60dnf.atc))
print(wilcox.test(spec10.ctrp.dtg,spec10.ctrpcisnf.dtg))
print(wilcox.test(spec10.ctrp.atc,spec10.ctrpcisnf.atc))
print(wilcox.test(spec10.nc60.dtg,spec10.nc60cisnf.dtg))
print(wilcox.test(spec10.nc60.atc,spec10.nc60cisnf.atc))



#print(wilcox.test(ari.clusters[1:40,1],ari.clusters[40:80,1]))

# 10 APC

apc10.ctrpdnf.dtg <- nRunClusterPerformance(ctrp.subdnf.dtg,loulabel.dtgbench.ctrp,method = "apcluster")
apc10.ctrpdnf.atc <- nRunClusterPerformance(ctrp.subdnf.atc,loulabel.atcbench.ctrp,method = "apcluster")
apc10.nc60dnf.dtg <- nRunClusterPerformance(nc60.subdnf.dtg,loulabel.dtgbench.nci60,method = "apcluster")
apc10.nc60dnf.atc <- nRunClusterPerformance(nc60.subdnf.atc,loulabel.atcbench.nci60,method = "apcluster")

print(wilcox.test(spec10.ctrp.dtg,apc10.ctrpdnf.dtg))
print(wilcox.test(spec10.ctrp.atc,apc10.ctrpdnf.atc))
print(wilcox.test(spec10.nc60.dtg,apc10.nc60dnf.dtg))
print(wilcox.test(spec10.nc60.atc,apc10.nc60dnf.atc))

ari.clusters <- data.frame( 
  ARI =c(spec10.ctrp.dtg,spec10.ctrp.atc,spec10.nc60.dtg,spec10.nc60.atc,
         spec10.ctrpdnf.dtg,spec10.ctrpdnf.atc,spec10.nc60dnf.dtg,spec10.nc60dnf.atc,
         apc10.ctrpdnf.dtg,apc10.ctrpdnf.atc,apc10.nc60dnf.dtg,apc10.nc60dnf.atc,
         spec10.ctrpcisnf.dtg,spec10.ctrpcisnf.atc,spec10.nc60cisnf.dtg,spec10.nc60cisnf.atc), 
  
  Benchmark = rep(rep(c("CTRP+Target","CTRP+ATC","NCI60+Target","NCI60+ATC"), each=folds),4), 
  Method=rep(c("Our method","DNF+spectral","DNF+Apcluster","CISNF+spectral"), each=folds*4) 
)

for(bench in c("CTRP+Target","CTRP+ATC","NCI60+Target","NCI60+ATC")){
  ari.sub <- ari.clusters[ari.clusters$Benchmark==bench,]
  ari.grouped<-group_by(ari.sub,Method)
  ari.summary <- summarise(ari.grouped,IQR=IQR(ARI),median=median(ARI))
  ari.summary$median<- round(ari.summary$median,digits = 4)
  ari.summary$IQR<- round(ari.summary$IQR,digits = 4)
  
  emf(file=paste("Output/",datestr,"/",bench,"ari_violin.emf",sep = ""), width=3, height = 2)
  
  
  p <- ggplot(ari.sub, aes(x=factor(Method), y=ARI))+geom_violin(aes(fill = factor(Method)),show.legend = FALSE)+
    geom_text(data = ari.summary, aes(x =factor(Method), y = median, label = median))+ 
    theme( axis.title.y=element_text(size=12),legend.title = element_text(size=12),
           legend.text = element_text(size=10),axis.text.y  = element_text(size=10))+
    scale_x_discrete("IQR",labels= ari.summary$IQR)
  print(p)
  
  dev.off()
}

emf(file=paste("Output/",datestr,"/ari_Compare_boxpolt.emf",sep = ""), width=2.5, height = 2.5)
p <- ggplot(ari.clusters, aes(x=factor(Benchmark), y=ARI))
p<-p + geom_boxplot(aes(fill = factor(Method)),show.legend = FALSE)

p + theme( axis.title.y=element_text(size=12),legend.title = element_text(size=12),
          legend.text = element_text(size=10),axis.text.y  = element_text(size=10),
          axis.title.x=element_blank(),
          axis.text.x=element_blank(),
          axis.ticks.x=element_blank())

dev.off()



PK <- seq(from = 20, to = 60, by = 1)


kper.ctrp.dtg<-clusterParamTune(ctrp,loulabel.dtgbench.ctrp,paramK = PK)
kper.ctrp.atc<-clusterParamTune(ctrp,loulabel.atcbench.ctrp,paramK = PK)

kper.nc60.dtg<-clusterParamTune(nc60,loulabel.dtgbench.nci60,paramK = PK)
kper.nc60.atc<-clusterParamTune(nc60,loulabel.atcbench.nci60,paramK = PK)


dat.paramk.ctrp <- data.frame( 
  performance =c(kper.ctrp.dtg$ari,kper.ctrp.dtg$silh,kper.ctrp.atc$ari), 
  K = c(kper.ctrp.dtg$params,kper.ctrp.dtg$params,kper.ctrp.atc$params), 
  benchmark=c(rep(c("ari-drug target","silhoutte","ari-ATC"), c(length(kper.ctrp.dtg$params),length(kper.ctrp.dtg$params),length(kper.ctrp.atc$params))) )
)

emf(file=paste("Output/",datestr,"/ctrp_sil_ari_curve_K.emf",sep = ""), width=2.5, height = 2.5)
# ggplot() + 
#   geom_line(aes(x=kper.ctrp.dtg$params,y=kper.ctrp.dtg$silh),color='orange',size=1.2) + 
#   geom_line(aes(x=kper.ctrp.dtg$params,y=kper.ctrp.dtg$ari),color='blue',size=1.2) + 
#   geom_line(aes(x=kper.ctrp.atc$params,y=kper.ctrp.atc$ari),color='green',size=1.2) + 
#   ylab('performance')+xlab('K')
p<-ggplot(data=dat.paramk.ctrp,
       aes(x=K, y=performance, colour=benchmark),size=1.1) +
        geom_line()

p + theme(axis.title.x =element_text(size=12), axis.title.y=element_text(size=12),legend.title = element_text(size=12),
          legend.text = element_text(size=10),axis.text.x  = element_text(size=10),axis.text.y  = element_text(size=10),
          legend.position ="none")
dev.off()


dat.paramk.nc60 <- data.frame( 
  performance =c(kper.nc60.dtg$ari,kper.nc60.dtg$silh,kper.nc60.atc$ari), 
  K = c(kper.nc60.dtg$params,kper.nc60.dtg$params,kper.nc60.atc$params), 
  benchmark=c(rep(c("ari-drug target","silhoutte","ari-ATC"), c(length(kper.nc60.dtg$params),length(kper.nc60.dtg$params),length(kper.nc60.atc$params))) )
)

emf(file=paste("Output/",datestr,"/nci60_sil_ari_curve_K.emf",sep = ""), width=2.5, height = 2.5)

# ggplot() + 
#   geom_line(aes(x=kper.nc60.dtg$params,y=kper.nc60.dtg$silh),color='orange',size=1.2) + 
#   geom_line(aes(x=kper.nc60.dtg$params,y=kper.nc60.dtg$ari),color='blue',size=1.2) + 
#   geom_line(aes(x=kper.nc60.atc$params,y=kper.nc60.atc$ari),color='green',size=1.2) + 
#   ylab('performance')+xlab('K')
p<-ggplot(data=dat.paramk.nc60,
       aes(x=K, y=performance, colour=benchmark),size=1.1) +
  geom_line()
p + theme(axis.title.x =element_text(size=12), axis.title.y=element_text(size=12),legend.title = element_text(size=12),
          legend.text = element_text(size=10),axis.text.x  = element_text(size=10),axis.text.y  = element_text(size=10),
          legend.position ="none")
dev.off()

# Apply the louvain on the benchmark nci60 dataset
loulabel.dtgbench.nci60<-affinClustering(dtgbench.nci60,method="louvain")
loulabel.atcbench.nci60<-affinClustering(ATCbench.nci60,method="louvain")
