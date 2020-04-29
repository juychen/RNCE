library(netcom)
library(apcluster)
library(xlsx)
library(igraph)
library(mclust)
library(netbiov)
library(devEMF)
library(SNFtool)
library(caret)

default<-function(value,default){
  # Funtion to get default params
  return(if(is.null(value)) default else value)
}

funTuneaspec<-function(x,k){
  result <- list(cluster=spectralClustering(x,K=k))
  
  return(result)
}

affinClustering<-function(affinity,method="spectral",...){
  # Clustering methods for affinity matrices
  kwargs <- list(...)
  
  # Generate the clustering using apclustering
  if(method=="apcluster"){
    
    # if(kwargs$K !=NULL){
    #   apresult <-apclusterK(affinity,K=kwargs$K)
    #   result<-labels(apresult,"enum")
    # }
    arg.q = default(kwargs$q,0.9)
    apresult <- apcluster(affinity, q=arg.q)
    result<-labels(apresult,"enum")
    
  # Generate the clustering using louvain  
  }else if(method=="louvain"){
    
    arg.mode = default(kwargs$mode,"undirected")
    arg.diag = default(kwargs$diag,FALSE)
    
    graph <- graph.adjacency(affinity, mode=arg.mode,diag=arg.diag)
    louvain<-cluster_louvain(graph)
    result<-louvain$membership
    
    
  # Generated the clustering using the spectural clustering  
  }else if(method=="spectral"){
    
    arg.K = default(kwargs$K,10)
    arg.type = default(kwargs$spectype,2)
    result <- spectralClustering(affinity,K=arg.K,type=arg.type)
    
  }
  names(result) <- rownames(affinity)
  return(result)
}

nRunClusterPerformance<-function(affinity,ground_truth,folds=10,method="spectral",truth_method = "louvain",verbose=FALSE,...){
  ### Compare clustering methods in multiple runs
  results <- list()
  kwargs <- list(...)
  names <- rownames(affinity)
  arg.replace = default(kwargs$replace,FALSE)
  
  # If the groud truth is labels
  if(is.null(dim(ground_truth))){
    ground_label <- ground_truth
  # If the groud truth is a adjecency matrix, the generate the result
  }else if(length(dim(ground_truth))>1){
    if(truth_method == "louvain"){
      graph.truth <- graph.adjacency(ground_truth, mode="undirected",diag=FALSE)
      louvain.truth<-cluster_louvain(graph.truth)
      ground_label<-louvain.truth$membership
    }else if(truth_method=="apcluster"){
      arg.truthq = default(kwargs$truthq,0.9)
      aptruth <- apcluster(ground_truth, q=arg.truthq)
      ground_label<-labels(aptruth,"enum")
    }
  }
  
  set.seed(default(kwargs$seed,0))
  if(folds!=FALSE){
    flds <- createFolds(names, k = folds, list = TRUE, returnTrain = TRUE)
  }
  
  # N runs of the clustering
  for(i in c(1:folds)){
    
    names.fold <- names[flds[[i]]]
    affin.fold <- affinity[names.fold,names.fold]
    label.fold <- ground_truth[names.fold]

    label <- affinClustering(affinity = affin.fold,method=method,...)
    ari<- adjustedRandIndex(label,label.fold)
    results[i] <- ari
  }
  
  if (verbose==TRUE) {
    print(results)
  }
  
  return(unlist(results))
}

labeltoCluster<-function(labels){
  # Convert label to cluster lists
  ClusterIdList <- list()
  ClusterNameList<- list()
  for(i in 1:max(labels)){
    xx <- names(labels[labels==i])
    ClusterNameList[[i]] <- xx
    ClusterIdList[[i]] <- labels[labels==i]
  }
  
  ClusterIdList<-Filter(length, ClusterIdList)
  ClusterNameList<-Filter(length, ClusterNameList)
  
  return(list(nameList=ClusterNameList,idList=ClusterIdList))
}


exeamplarCluster<-function(clusters,simMat){
  # Convert exemplars to clusters
  exemplars<-list()
  for(c in clusters){
    
    if (length(c)>1){
      
      clusterSim <- simMat[c,c]
      sumSim<-colSums(clusterSim)
      exemplar<-names(sumSim[sumSim == max(sumSim)])
      exemplars<-append(exemplars,exemplar[1])
    
    }else{
      exemplars<-append(exemplars,c[[1]])
    }
  }
  
  return(exemplars)
}



clusterParamTune<-function(affinity,ground_truth,paramK=c(30,40,50,60),method="spectral",truth_method = "louvain",verbose=FALSE,...){
  # Tune parameters
  silhs <- list()
  aris <- list()
  
  kwargs <- list(...)
  names <- rownames(affinity)
  set.seed(default(kwargs$seed,0))
  # N runs of the clustering
  for(k in paramK){
    
    
    label <- affinClustering(affinity = affinity,method=method,K=k,...)
    label.sub <- label[names(ground_truth)]
    ari<- adjustedRandIndex(label.sub,ground_truth)
    silh <- silhouette(label,1/affinity)
    aris <- append( aris,ari)
    silhs <- append( silhs,mean(silh[,3]))
    
    #print()
    
  }
  
  if (verbose==TRUE) {
    print((unlist(aris)))
    print((unlist(silhs)))
    print((unlist(silhs))*(unlist(aris)))
    
  }
  results <- data.frame(ari=(unlist(aris)),silh=(unlist(silhs)),params=paramK)
  return((results))
}
