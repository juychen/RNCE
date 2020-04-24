#################
### Code to construct a network from community clustering and keeping exemplar drugs as nodes
#################

library(igraph)
library(netbiov)
library(apcluster)
library(xlsx)
library(RColorBrewer)
library(devEMF)

source("Rcode/affinityClustering.R")



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

## note: the main scripts: main-nci10-lincs.R and main-ctrpv-lincs.R must be run prior to running this script
#load("Data/ctrpv2-Integrated.RData")
load("Data/integrationSimilarityMat-ctrpv2.RData")
load("Data/sensSimilarityMat-ctrpv2.Rdata")
load("Data/strcSimilarityMat-ctrpv2.RData")
load("Data/pertSimilarityMat-ctrpv2.RData")


combi <- integrtStrctSensPert

############### transform the similarity network into a distance matrix where the diagonals are 0s and similarity values are transformed into shortest distances
combiall <- combi
diag(combiall) <- 0
cormd <- 0.01/combiall
#cormd[cormd>50] <-50

#cormd <- exp(0.2/combiall)

diag(cormd) <- 0
colnames(cormd) <- colnames(combiall)
rownames(cormd) <- rownames(combiall)

##### if you have forget to run the apclust object
#set.seed(12345)
apcomb.back <- apcluster(combi, q=0.9)



# apcomb <- apcluster(combi, q=0.9)
# apcombstrc <- apcluster(strcAffMat, q=0.9)
# apcombpert <- apcluster(pertAffMat, q=0.9)
# apcombsens <- apcluster(sensAffMat, q=0.9)


apcomb<-affinClustering(combi,K=35)
cluster.labels <- labeltoCluster(apcomb)

exemplars<-exeamplarCluster(cluster.labels$nameList,combi) 
###### names of the exemplar drugs
#exemp <- names(apcomb@exemplars)
exemp <- unlist(exemplars)
########## build the network

graph <- graph.adjacency(cormd, weighted=TRUE, mode="undirected")
dfd <- get.data.frame(graph, what="edges")
# vts <- get.data.frame(graph, what="vertices")
# 
# write.csv(dfd, file = "Output//edges.csv")
# write.csv(vts, file = "Output//vertices.csv")

nexemp <- dfd$from %in% exemp   &  dfd$to %in%  exemp 
dfd <- dfd[nexemp,]

dataframe2adjacency<-function(dfd, names.var=NULL){
  if(is.null(names.var)){
    names.var <- sort(union(dfd[,1],dfd[,2]))
  }
  adj <- matrix(0,nrow=length(names.var),ncol=length(names.var),dimnames=list(names.var,names.var))
  for(i in 1:nrow(dfd)){
    # print(dfd[i,3])
    adj[as.character(dfd[i,1]),as.character(dfd[i,2])]<-dfd[i,3]
  }
  return(adj)
}
mynet <- dataframe2adjacency(dfd)
g <- graph.adjacency(mynet,weighted=T,mode="undirected")

# Edited by jy to generate colorful graphs
#mycols <- rep("darkgrey",nrow(mynet))
ncolors <- nrow(mynet)
qual_col_pals <- brewer.pal.info[brewer.pal.info$category == 'qual',]
col_vector <- unlist(mapply(brewer.pal, qual_col_pals$maxcolors, rownames(qual_col_pals)))
#pie(rep(1,ncolors), col=sample(col_vector, ncolors))
mycols<-col_vector[1:ncolors]

names(mycols) <- rownames(mynet) 
v.colors <- mycols

###### order the names of exemplars according to network g and assign size of nodes based on number of drugs in each community
#names(apcomb@clusters) <- exemp
#apres2 <- apcomb@clusters[order(names(apcomb@clusters))]

names(cluster.labels$nameList) <- exemp
apres2 <- cluster.labels$nameList[order(names(cluster.labels$nameList))]

ll <- unlist(lapply(apres2, function(x) length(apres2[x])))

set.seed(3878)
#windowsFonts(A = windowsFont("Sans"))
pdf.options(family = "Helvetica")

emf(paste(subDir,"plot_network_ctrpv2_targchembl_exemp.emf",sep="/"),height=5.5,width=5.5,family  = "Times")
#emf("plot_network_ctrpv2_targchembl_exemp.emf",height=10,width=10,emfPlus = FALSE)
par(mar=c(0,0,0,0.5))

plot<-mst.plot(g, mst.edge.col="darkgrey",
         mst.e.size=1.2,
         vertex.color=mycols, colors="gray97",tkplot=FALSE, bg=NA, v.size=ll,
         layout.function=layout_with_lgl,
         #layout.function=layout.auto,
         e.arrow=0, e.size=(dfd$weight <= 15),#v.lab=TRUE,
         v.lab =names(apres2),
         v.lab.col="black", lab.dist=0.5, v.lab.cex=0.8)

#par(family="Helvetica", font=1)   


#mst.plot(g)
dev.off()

emf(paste(subDir,"plot_network_ctrpv2_targchembl_exemp2.emf",sep="/"),height=5.5,width=5.5,family  = "Times")
#emf("plot_network_ctrpv2_targchembl_exemp.emf",height=10,width=10,emfPlus = FALSE)
par(mar=c(0,0,0,0.5))

plot<-mst.plot(g, mst.edge.col="darkgrey",
               mst.e.size=1.2,
               vertex.color=mycols, colors="gray97",tkplot=FALSE, bg=NA, v.size=ll,
               layout.function=layout_with_lgl,
               #layout.function=layout.auto,
               e.arrow=0, e.size=(dfd$weight <= 15),#v.lab=TRUE,
               v.lab = FALSE,
               v.lab.col="black", lab.dist=0.5, v.lab.cex=0.8)

#par(family="Helvetica", font=1)   


#mst.plot(g)
dev.off()
