library(apcluster)
source('RCode/affinityClustering.R')

communityGen <- function(combi, dname, GMT_TARG, method = "apcluster",...) {
  
  kwargs <- list(...)
  
## all communities from integrative analysis results for nci60/l1000 or ctrpv2/l1000 
## Modified by JY to accept spectral result
if(method == "apcluster"){
  apcomb <- apcluster(combi, q=0.9)
  
  if(dname =="ctrpv2"){
    save(apcomb, file="Data/CTRPv2comm.Rdata")
  }else if(dname =="nci60"){
    save(apcomb, file="Data/NCI60comm.RData")
  }
  
  ll <- list()
  for(i in 1:length(apcomb)){
    xx <- names(apcomb[[i]])
    ll[[i]] <- xx
  }
  
  Clust <- apcomb@clusters
}else{
  
  comm <- affinClustering(combi,method = method,...)
  
  if(dname =="ctrpv2"){
    save(comm, file="Data/CTRPv2commSpec.Rdata")
  }else if(dname =="nci60"){
    save(comm, file="Data/NCI60commSpec.RData")
  }
  
  
  ll <- list()
  names <- list()
  Clust <- list()
  
  for(i in 1:max(comm)){
    xx <- names(comm[comm==i])
    ll[[i]] <- xx
    #ll[[i]] <- comm[comm==i]
    #names[[i]] <- xx[[1]]
    Clust[[i]] <- comm[comm==i]
    
  }
}  


indx <- sapply(ll, length)
#indx <- lengths(lst) 
res <- as.data.frame(do.call(rbind,lapply(ll, `length<-`,max(indx))))
llx <- data.frame("population"=1:length(rownames(res)), "number of drugs"=indx,res)
llx <- llx[llx$number.of.drugs>0,]
if(method=="apcluster"){
  row.names(llx) <- names(apcomb@exemplars)
  
}else{
  # exemplars
  communities<-labeltoCluster(comm)
  exemplars<-exeamplarCluster(communities$nameList,combi)
  row.names(llx) <- unlist(exemplars)
  
}
filename = paste(getwd(), "/Output/", "communities_combi_", dname, ".csv", sep="")
write.csv(llx, filename, row.names=TRUE)
dim(res)

ll <- list()
for(i in 1:length(GMT_TARG)){
   xx <- GMT_TARG[[i]]
   ll[[i]] <- xx
}

indx <- sapply(ll, length)
#indx <- lengths(lst) 
res <- as.data.frame(do.call(rbind,lapply(ll, `length<-`,max(indx))))
llx <- data.frame("population"=1:length(rownames(res)), "number of drugs"=indx,res)
row.names(llx) <- unlist(names(GMT_TARG))
#write.csv(llx," GMT_targ_chembl.csv", row.names=TRUE)
dim(res)


## Keep communities with at least 2 drugs showing a known mechanism of action from GMT
GMT_TARG2 <- c(as.character(unlist(GMT_TARG)))

#Refine the Clusters: exclude drugs from each cluster which are not found in the GMT benchmark drugs
#Clust <- apcomb@clusters
ClusterList <- lapply(Clust,function(x){if (length(intersect(GMT_TARG2 ,names(x[]))) >= 2) {names(x[])} else {NULL}})
ClusterByCondition<- sapply(ClusterList, function(x) length(x) > 1) 
ClusterListRefined <- ClusterList[ClusterByCondition]

indx <- sapply(ClusterListRefined, length)
res <- as.data.frame(do.call(rbind,lapply(ClusterListRefined, `length<-`,max(indx))))
llx <- data.frame("population"=1:length(rownames(res)), "number of drugs"=indx,res)
filename = paste(getwd(), "/Output/", "clusterListRefined_combi_", dname, ".csv", sep="")
write.csv(llx, filename, row.names=TRUE)
}

