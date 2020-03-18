# # #### Testing snippet
# a <- t(matrix(c(1.0,0.4,0.2,0.2,0.7,0.8,
#                 0.4,1.0,0.2,0.5,0.4,0.2,
#                 0.2,0.2,1.0,0.6,0.5,0.3,
#                 0.2,0.5,0.6,1.0,0.6,0.1,
#                 0.7,0.4,0.5,0.6,1.0,0.7,
#                 0.8,0.2,0.3,0.1,0.7,1.0),nrow=6,ncol=6))
# #### Testing snippet
# a <- t(matrix(c(1.0,0.4,0.2,0.2,0.7,0.8,
#                 1.0,0.4,0.2,0.2,0.7,0.9,
#                 0.2,0.2,1.0,0.6,0.5,0.3,
#                 0.2,0.5,0.6,1.0,0.6,0.1,
#                 0.7,0.4,0.5,0.6,1.0,0.7,
#                 0.8,0.2,0.3,0.1,0.7,1.0),nrow=6,ncol=6))
# a <- t(matrix(c(1,1,1,0,0,
#                 1,1,0,0,1,
#                 1,0,1,1,0,
#                 0,0,1,1,0,
#                 0,1,0,0,1
#                 ),nrow=5,ncol=5))
# b <- t(matrix(c(1.0,0.8,0.7,0.0,0.1,
#                 0.8,1.0,0,0,1,
#                 0.7,0,1,1,0,
#                 0.0,0,1,1,0,
#                 0.1,1,0,0,1
# ),nrow=5,ncol=5))
#options(digits=5,scipen = 15)

# scaa <- sca(a, k=3)
# wlcea <-  wlce(scaa, k=3)
# D <- modjaccard(lcea)
# round(D,5) 


.kneighbysim <- function(xx,KK=20) {
  ###This function outputs the KK nearest neighbors without normalization.	
  
  zero <- function(x) {
    s = sort(x, index.return=TRUE)
    x[s$ix[1:(length(x)-KK)]] = 0
    return(x)
  }
  #normalize <- function(X) X / rowSums(X)
  A = matrix(0,nrow(xx),ncol(xx));
  for(i in 1:nrow(xx)){
    A[i,] = zero(xx[i,]);
    
  }
  
  return(A)
}

.argkneighbysim <- function(xx,KK=20) {
  ###This function outputs the indicator matrix of top KK neighbors.	
  A = matrix(FALSE,nrow(xx),ncol(xx));
  
  # zero <- function(x) {
  #   # Sort is in decending order
  #   s = sort(x, index.return=TRUE)
  #   x[s$ix[1:(length(x)-KK)]] = 0
  #   x[s$ix[(length(x)-KK)+1:(length(x)-1)]] = 1
  # 
  #   return(x)
  # }
  
  for(i in 1:nrow(xx)){
    s <- order(xx[i,],decreasing = TRUE)
    #s <- sort(xx[i,], index.return=TRUE)
    A[i,s[1:KK]] <- TRUE
  }
  
  return(A)
}

.kre_neighbysim <- function(x,KK=20) {
  ###This function outputs the indicator matrix of top KK R-neighbors.	
  
  xx <- matrix(x,nrow = nrow(x),ncol = ncol(x))
  #diag(xx) <- 0
  
  argknxx <- .argkneighbysim(xx,KK=KK)
  krn_indicator <- argknxx * t(argknxx)
  return(krn_indicator)
}

# normalize <- function(x,type='mean'){
#   if(type == 'mean'){
#     return(x/ sum(x))
#   }else if(type == 'softmax'){
#     return(exp(x)/sum(exp(x)))
#   }
# }

normalize <- function(X,type='mean'){
  #Normalization method for affinity matrices
  
  
  row.sum.mdiag <- rowSums(X) - diag(X) 
  #If rowSumx(X) == diag(X), set row.sum.mdiag to 1 to avoid div by zero
  row.sum.mdiag[row.sum.mdiag == 0] <- 1   
  
  if(type =='mean'){

    X <- X/(2*(row.sum.mdiag))
    diag(X) <- 0.5
  }else if(type == 'softmax'){
    X <- exp(X)/(2*exp((row.sum.mdiag)))
    diag(X) <- 0.5
  }

  return(X)
}

sca <- function(amat, k , neighbor = 'knn'){
    ###This function outputs the sca of the affinenaty matrix

    exp_amat <- exp(-1/amat)
    if(neighbor =='knn'){
      neig_mask_amat <- .argkneighbysim(amat, KK=k)
      
    }
    else if(neighbor =='krnn'){
      neig_mask_amat <- k_reciprocal_nn(amat,k)
    }

    return((neig_mask_amat*exp_amat))
}

lce <- function(fmat,k){
    ## This is a local consistency enhancement snippet
    
    #normalize <- function(X) X / rowSums(X)
    #diag(fmat) <-0

    G = matrix(0,nrow(fmat),ncol(fmat));

    for(i in 1:nrow(fmat)){
        sorted_row = sort(fmat[i,], index.return=TRUE)
        idx_neigh <- tail(sorted_row$ix, n=k)
        #print(idx_neigh)

        if(k==1){
            #G[i,] <- normalize(matrix(fmat[idx_neigh,],nrow=1),type = 'softmax')
            G[i,] <- (matrix(fmat[idx_neigh,],nrow=1))
        } else{
            gi <- colMeans(fmat[idx_neigh,])
            #G[i,] <- normalize(matrix(gi,nrow = 1),type = 'softmax')
            G[i,] <- (matrix(gi,nrow = 1))
        }
    }
    
    G <- (G + t(G))/2
    
    G<-normalize(G,'mean')
    return(G)
}

modjaccard <- function(mat,parallel=TRUE){
  ### This the modified jaccard distance of a similarity matrix
  D = matrix(0,nrow(mat),ncol(mat));
  
  if(parallel==FALSE){
    for(p in 1:nrow(mat)){
      for(q in 1:ncol(mat)){
        if(p==q){
          D[p,q] <- 0
        }else{
          D[p,q] <- 1-(sum(apply(mat[c(p,q),],2,min))/sum(apply(mat[c(p,q),],2,max)))
        }
      }
    }
  }else{
    
    print("Calculating distances in parallel")
    
    library(doParallel)
    cores=detectCores()
    cl <- makeCluster(cores[1]) #not to overload your computer
    registerDoParallel(cl)
    
    print("Clusters registered")
    
    i <- 1:(nrow(mat)*ncol(mat))
    mata <<- matrix(i,nrow=nrow(mat),ncol=ncol(mat))
    mat <<- mat
    
    fun <- function(i){
      idx2<-which(mata==i,arr.ind=TRUE)
      p <- idx2[1]
      q <- idx2[2]
      return(1-(sum(apply(mat[c(p,q),],2,min))/sum(apply(mat[c(p,q),],2,max))))
    }
    
    clusterExport(cl, varlist = c("mat","mata"))
    
    print("Start calculating")
    
    res <- parLapply(cl,i,fun)
    D <- matrix(unlist(res),nrow=nrow(mat),ncol=ncol(mat))
    print("Finished calculating distanc in parallel")
    
    stopCluster(cl)
    print("Clusters stopped")
    
    rm(mata,envir = globalenv())
    rm(mat,envir = globalenv())
    
  }
  
  return(D)
    
}


k_reciprocal_nn <-function(amat,KK){
  # Calculated the expaned k reciprocal neighbor
  xx <- amat

  Result = matrix(FALSE,nrow(amat),ncol(amat));
  
  
  # Get k and k/2 recoiprocal neighbors of xx
  mat_ind_krecipnn_xx <-.kre_neighbysim(xx,KK=KK) 
  mat_ind_halfkrecipnn_xx <-.kre_neighbysim(xx,KK=ceiling(KK/2))
  
  # Get the intersection for each nondes of the k and k/2 recoiprocal nearest neighbors
  mat_ind_intersect_krnn_xx <- mat_ind_krecipnn_xx %*% t(mat_ind_halfkrecipnn_xx)
  

  # Get the size of each half-k recoiprocal nearest neighbors set
  size_halfkrecipnn_xx <- (colSums(mat_ind_halfkrecipnn_xx))
  mat_size_halfkrecipnn_xx <- t(replicate(nrow(xx),size_halfkrecipnn_xx))
  
  # Record whihch neigbhour match the critetira of the k expanded recoiprocal nearest neighbors
  mat_ind_neigbor_match<-(mat_ind_intersect_krnn_xx > 2/3*mat_size_halfkrecipnn_xx)
  
  for(i in 1:nrow(mat_ind_halfkrecipnn_xx)){
      r <- mat_ind_halfkrecipnn_xx[i,]
      idx_krn <- mat_ind_neigbor_match[i,]
      
      if(sum(idx_krn)<1) {
        vec_i_rneigh_intersections <- mat_ind_krecipnn_xx[i,]
      }
      
      else if(sum(idx_krn)>1){
        vec_i_rneigh_intersections  <- colSums(mat_ind_halfkrecipnn_xx[idx_krn,])
        
      }
      else{
        vec_i_rneigh_intersections  <- (mat_ind_halfkrecipnn_xx[idx_krn,])
        
      }
      Result[i,] <- vec_i_rneigh_intersections | mat_ind_krecipnn_xx[i,]
  } 
  
  #diag(Result) <- FALSE
  return(Result)
  
}


range01 <- function(XX, ...){
  # Scale the matrix from 0 to 1
  X <- as.vector(XX)
  rangex <- ( (X - min(X, ...)) / (max(X, ...) - min(X, ...)))
  result <- matrix(rangex,nrow=nrow(XX),ncol=ncol(XX),byrow=TRUE)
  diag(result) <- 1 
  
  colnames(result) <- rownames(result) <- colnames(XX)
  return(result)
}

zscore <-function(XX, ...){
  # Nomralize by z score
  X <- as.vector(XX)
  rangex <- ( (X - mean(X, ...)) / (sd(X)))
  result <- matrix(rangex,nrow=nrow(XX),ncol=ncol(XX),byrow=TRUE)
  colnames(result) <- rownames(result) <- colnames(XX)
  diag(result) <- 1 
  colnames(result) <- rownames(result) <- colnames(XX)
  
  return(result)
}

affinityKRNNMatrix <- function(diff,K=20,sigma=0.5) {
  # Computes the affinity matrix for a given distance matrix
  #
  # Args:
  #   diff: Distance matrix
  #   K: Number of Rreciprocal nearest neighbours to sparsify similarity
  #   sigma: Variance for local model
  #
  # Returns:
  #   Affinity matrix using exponential similarity kernel scaled by k nearest
  #       neighbour average similarity
  #

  N <- nrow(diff)

  diff <- (diff + t(diff)) / 2
  diag(diff) <- 0
  #sortedColumns <- as.matrix(t(apply(diff,2,sort)))
  sortedColumns <- as.matrix(diff)
  
  mask_KRNN <-k_reciprocal_nn(diff,KK = K)
  mask_KRNN[mask_KRNN==0] <- NA
  
  finiteMean <- function(x){
    return(mean(x[is.finite(x)]))
  }
  
  sortedColumns <- sortedColumns * mask_KRNN
  means <- apply(sortedColumns,1,finiteMean)+.Machine$double.eps;

  avg <- function(x,y){
    return((x+y)/2)
  }
  Sig <- outer(means,means,avg)/3*2 + diff/3 + .Machine$double.eps;
  Sig[Sig <= .Machine$double.eps] <- .Machine$double.eps
  densities <- dnorm(diff, 0, sigma*Sig, log = FALSE)

  W <- (densities + t(densities)) / 2
  return(W)
}

k_renei_rerank <- function(mat,k1=20,k2=20){
  
  # RNCE method of calculate the contextual information using reciprocal neighbor and reranking 
  mask <- .kre_neighbysim(mat,KK=k1)
  scam<- mat*mask
  scam <- normalize(scam)
  
  lcem <- lce(scam,k=k2)
  lcem <- (lcem+t(lcem))/2
  lcem <- normalize(lcem)
  
  Dj <- modjaccard(lcem,TRUE)
  DJ <- (Dj +t(Dj))/2
  
  result <- DJ
  result <- (result +t(result))/2
  
  
  
  return(result)
  
}


CISNF <- function(mat,k1=20,k2=20){
  #CISNF implementation
  mask <- .kneighbysim(mat,KK=k1)
  
  scam<- mat*mask
  scam <- normalize(scam)
  
  lcem <- lce(scam,k=k2)
  lcem <- (lcem +t(lcem))/2
  
  modj <- modjaccard(lcem,parallel = TRUE)
  
  
  result <- 1-modj
  
  result <- (result + t(result))/2
  result<-normalize(result)
  
  return(result)
  
}