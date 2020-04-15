###############################################################################################################
## Function reads in the affinity matrices or network layers of sensitivity, structure and perturbation data, 
## and integrates them via SNF (similarity network fusion) method 
## input: 
##     sensAff: affinity matrix generated in the "constSensitivityLayer"
##     strcAff: affinity matrix generated in the "constStructureLayer"
##     pertAff: affinity matrix generated in the "constPerturbationLayer"
##     
## output: 
##     integartion result (matrix)  
##
## 
###############################################################################################################


integrateStrctSensPert <- function(sensAff, strcAff, pertAff) {
  ## Original DNF
  integration <- SNFtool::SNF(list(sensAff, strcAff, pertAff), K = 20, t =20)
  colnames(integration) <- rownames(integration) <- colnames(strcAff)
  
  return(integration)
  
}

integrateStrctSensPertR <- function(sensAff, strcAff, pertAff,K=20,t=20) {
  ## Reciprocal neighbor DNF
  integration <- RSNF(list(sensAff, strcAff, pertAff), K = K, t =t)
  colnames(integration) <- rownames(integration) <- colnames(strcAff)
  
  return(integration)
  
}

integrateStrctSensPert_RNCE <- function(sensAff, strcAff, pertAff,k1=100,k2=5,saverr=TRUE) {
  # RNCE
  start.time <- Sys.time()
  
  # Get the integration of the three types of distances
  integration <- integrateStrctSensPert(sensAff, strcAff, pertAff)
  
  if(saverr!=FALSE){
    # Save contextual matrix
    badcs <- "[\xb5]|[\n]|[,]|[;]|[:]|[-]|[+]|[*]|[%]|[$]|[#]|[{]|[}]|[[]|[]]|[|]|[\\^]|[/]|[\\]|[.]|[_]|[ ]"
    end.time <- gsub(badcs, "",  toString(end.time))
    save(re_ranked, file=paste("Data/initSimilarityMat-",saverr#,end.time
                               ,".Rdata",sep = ""))
    
  }

  # Rranking accoding to the integrated result
  re_ranked <- 1-k_renei_rerank(integration,k1=k1,k2=k2)
  colnames(re_ranked) <- rownames(re_ranked) <- colnames(integration)
  re_ranked <- normalize(re_ranked)
  
  
  # Fusion the reranking result and the integration result
  result <- SNFtool::SNF(list(integration,re_ranked), K = 20, t =20)

  colnames(result) <- rownames(result) <- colnames(integration)
  
  end.time <- Sys.time()
  
  print(end.time - start.time)
  
  if(saverr!=FALSE){
    # Save contextual matrix
    badcs <- "[\xb5]|[\n]|[,]|[;]|[:]|[-]|[+]|[*]|[%]|[$]|[#]|[{]|[}]|[[]|[]]|[|]|[\\^]|[/]|[\\]|[.]|[_]|[ ]"
    end.time <- gsub(badcs, "",  toString(end.time))
    save(re_ranked, file=paste("Data/rerkSimilarityMat-",saverr#,end.time
                               ,".Rdata",sep = ""))
    
  }
  
  
  return(result)
  
}

integrateStrctSensPert_RSNF_RNCI <- function(sensAff, strcAff, pertAff,k1=100,k2=5,saverr=TRUE) {
  # Recrprocal DNF + Contextual information
  start.time <- Sys.time()
  
  
  # Get the integration of the three types of distances
  integration <- integrateStrctSensPertR(sensAff, strcAff, pertAff)

  # Rranking accoding to the integrated result
  re_ranked <- 1-k_renei_rerank(integration,k1=k1,k2=k2)
  colnames(re_ranked) <- rownames(re_ranked) <- colnames(integration)
  re_ranked <- normalize(re_ranked)
  
  
  # Fusion the reranking result and the integration result
  result <- SNFtool::SNF(list(integration,re_ranked), K = 20, t =20)

  colnames(result) <- rownames(result) <- colnames(integration)
  
  end.time <- Sys.time()
  
  print(end.time - start.time)
  
  if(saverr!=FALSE){
    badcs <- "[\xb5]|[\n]|[,]|[;]|[:]|[-]|[+]|[*]|[%]|[$]|[#]|[{]|[}]|[[]|[]]|[|]|[\\^]|[/]|[\\]|[.]|[_]|[ ]"
    end.time <- gsub(badcs, "",  toString(end.time))
    save(re_ranked, file=paste("Data/rerkSimilarityMat-",saverr#,end.time
                               ,".Rdata",sep = ""))
    
  }
  
  
  return(result)
  
}


integrateStrctSensPert_cisnf <- function(sensAff, strcAff, pertAff,k1=100,k2=5) {
  # CISNF
  integration <- integrateStrctSensPert(sensAff, strcAff, pertAff)
  cisnf_i <- CISNF(integration,k1=k1,k2=k2)
  colnames(cisnf_i) <- rownames(cisnf_i) <- colnames(integration)
  return(cisnf_i)
  
}
