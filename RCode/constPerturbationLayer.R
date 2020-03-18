###############################################################################################################
## Function reads in the perturbation data, and generates "affinity matrix" for the perturbation layer to be used either 
## solely or in combination of other layers (see function "integrateStrctSensPert") 
## input: 
##     sensDat: perturbation data processed in "perturbationData" function 
##     
## output: 
##     affinity matrix 
##
## 
###############################################################################################################


constPerturbationLayerDNF <- function(pertDat,method="pearson") {

   # Correlation for Perturbation
   pertCor <- cor(pertDat, method = method, use = "pairwise.complete.obs")
   ## Calculate affinity matrix (from generic distance) as described in SNFtool package with default values
   pertAff <- SNFtool::affinityMatrix(1-pertCor, 20, 0.5)

   return(pertAff)

}


constPerturbationLayer <- function(pertDat,metric="eJaccard",K = 20) {
  
  pertCor <- as.matrix(dist(t(pertDat),method = metric, diag = TRUE, upper = FALSE, p = 1))
  #expCor <- exp(-pertCor)
  expCor <- SNFtool::affinityMatrix(pertCor, K, 0.5)
  
  #pertAff <- k_reciprocal_nn(expCor,KK = K) * expCor
  
  pertAff <- expCor
  #pertAff <- SNFtool::affinityMatrix(1-pertCor, K, 0.5)
  
  
  return(pertAff)
  
}

