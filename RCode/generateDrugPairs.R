###############################################################################################################
## Function intersects the benchmark data and the layers of the network and then generates drug pairs for the 
## final (intersection) subset of drugs
##
## input: 
##     benchDat
##     benchDat
##     strcAff
##     sensAff
##     pertAff
##     integration  
##
## output: 
##     list of drug pairs from benchmark, structure, perturbation, sensitivtiy and the integration of all
##
## 
###############################################################################################################

generateDrugPairs <- function(benchDat, strcAff, sensAff, pertAff, integration, iorio, iskar, superpred, drugerank, dnf = NULL) {

    
    ## intersection of the bechmark set and the network layers
    intx <- intersect(colnames(benchDat), colnames(strcAff))
    dataStr <- strcAff[intx,intx]
    dataStr[upper.tri(dataStr, diag=TRUE)] <- NA
    dataSens <- sensAff[intx,intx]
    dataSens[upper.tri(dataSens, diag=TRUE)] <- NA
    dataPert <- pertAff[intx,intx]
    dataPert[upper.tri(dataPert, diag=TRUE)] <- NA
    dataIntegrAll <- integration[intx,intx]
    dataIntegrAll[upper.tri(dataIntegrAll, diag=TRUE)] <- NA
    
    # Add bench mark for dnf
    if(is.null(dnf)==FALSE){
      LDNFS <- length(dnf)
      
      dataDnfAll <- dnf
      
      for(i in 1:LDNFS){
        mati <- dnf[[i]]
        mati <- mati[intx,intx]
        mati[upper.tri(mati, diag=TRUE)] <- NA
        dataDnfAll[[i]] <- mati
      }
      
    }
    
    ind <- match(colnames(benchDat), colnames(iorio))
    if (!all(colnames(iorio)[ind] == colnames(benchDat))) { stop("error!")}
    data <- iorio[ind,ind]
    # #data <- iorio
    if (!all(colnames(data) == colnames(benchDat))){ stop("error!")}
    data[upper.tri(data, diag=TRUE)] <- NA
    
    indx <- match(colnames(benchDat), colnames(iskar))
    if (!all(colnames(iskar)[indx] == colnames(benchDat))) { stop("error!")}
    datax <- iskar[indx,indx]
    # #datax <- iskar
    if (!all(colnames(datax) == colnames(benchDat))){ stop("error!")}
    datax[upper.tri(datax, diag=TRUE)] <- NA
    
    
    if (class(superpred)==class(iorio)) {
        ind2 <- match(colnames(benchDat), colnames(superpred))
        if (!all(colnames(superpred)[ind2] == colnames(benchDat))) { stop("error!")}
        data2 <- superpred[ind2,ind2]
        #data2 <- SuperPredsimil
        all(colnames(data2) == colnames(benchDat))
        data2[upper.tri(data2, diag=TRUE)] <- NA
    }
    
    
    if (class(drugerank)==class(iorio)) {
        ind2 <- match(colnames(benchDat), colnames(drugerank))
        if (!all(colnames(drugerank)[ind2] == colnames(benchDat))) { stop("error!")}
        data2 <- drugerank[ind2,ind2]
        #data2 <- SuperPredsimil
        all(colnames(data2) == colnames(benchDat))
        data2[upper.tri(data2, diag=TRUE)] <- NA
    }
    
    
    ##create pairs of drugs from benchmark 1 if same drug set from GMT and 0 otherwise
    benchDat[upper.tri(benchDat, diag=TRUE)] <- NA
    benchPairs <- melt(benchDat)
    benchPairs <- na.omit(benchPairs)
    colnames(benchPairs)[3] <- "bench"
    
    ## create pairs of drugs from structure
    strcPairs <- melt(dataStr)
    strcPairs <- na.omit(strcPairs)
    colnames(strcPairs)[3] <- "obs.str"
    
    ## create pairs of drugs from sensitivty
    sensPairs <- melt(dataSens)
    sensPairs <- na.omit(sensPairs)
    colnames(sensPairs)[3] <- "obs.sens"
    
    ## create pairs of drugs from perturbation
    pertPairs <- melt(dataPert)
    pertPairs <- na.omit(pertPairs)
    colnames(pertPairs)[3] <- "obs.pert"
    
    ## create pairs of drugs from combination
    integrPairs <- melt(dataIntegrAll)
    integrPairs <- na.omit(integrPairs)
    colnames(integrPairs)[3] <- "obs.combiall"
    
    # ###
    iorioPairs <- melt(data)
    iorioPairs <- na.omit(iorioPairs)
    colnames(iorioPairs)[3] <- "obs.iorio"
    # #
    #
    iskarPairs <- melt(datax)
    iskarPairs <- na.omit(iskarPairs)
    colnames(iskarPairs)[3] <- "obs.iskar"
    #
    
    ## create pairs of drugs from dnf
    if(is.null(dnf)==FALSE){
      integrDNFPairs <- list()
      for(i in 1:LDNFS){
        integrDNFPairs[[i]]<- melt(dataDnfAll[[i]] )
        integrDNFPairs[[i]]<- na.omit(integrDNFPairs[[i]])
        colnames(integrDNFPairs[[i]])[3] <- paste("obs.dnfall",i,sep='.')
      }
      # integrDNFPairs <- melt(dataIntegrAll)
      # integrDNFPairs <- na.omit(integrDNFPairs)
      # colnames(integrDNFPairs)[3] <- "obs.dnfall"
    }

    
    if (class(superpred)==class(iorio)) {
        superPairs <- melt(data2)
        superPairs <- na.omit(superPairs)
        colnames(superPairs)[3] <- "obs.superPred"
        
        res <- list(benchPairs=benchPairs, strcPairs=strcPairs, sensPairs=sensPairs, pertPairs=pertPairs, integrPairs=integrPairs, iorioPairs=iorioPairs, iskarPairs=iskarPairs, superPairs=superPairs, drugePairs=NULL)
        
        if(is.null(dnf)==FALSE){
          for(i in 1:LDNFS){
            pairname <- paste("dnf",i,"Pairs",sep='')
            res[[pairname]] <- integrDNFPairs[[i]]
          }
        }
        return(res)
    }
    
    
    if (class(drugerank)==class(iorio)) {
        drugePairs <- melt(data2)
        drugePairs <- na.omit(drugePairs)
        colnames(drugePairs)[3] <- "obs.drugerank"
        res <- list(benchPairs=benchPairs, strcPairs=strcPairs, sensPairs=sensPairs, pertPairs=pertPairs, integrPairs=integrPairs, iorioPairs=iorioPairs, iskarPairs=iskarPairs, superPairs=NULL, drugePairs=drugePairs)
        if(is.null(dnf)==FALSE){
          for(i in 1:LDNFS){
            pairname <- paste("dnf",i,"Pairs",sep='')
            res[[pairname]] <- integrDNFPairs[[i]]
          }
        }
        return(res)
        
    }
}

###############################################################################################################
## Function intersects the benchmark data and the layers of the network and then generates drug pairs for the 
## final (intersection) subset of drugs
##
## input: 
##     benchDat
##     benchDat
##     strcAff
##     sensAff
##     pertAff
##     integration  
##
## output: 
##     list of drug pairs from benchmark, structure, perturbation, sensitivtiy and the integration of all
##
## 
###############################################################################################################
generateIntDrugPairs <- function(benchDat, integration) {
  
  
  ## intersection of the bechmark set and the network layers
  intx <- intersect(colnames(benchDat), colnames(integration))
  dataIntegrAll <- integration[intx,intx]
  dataIntegrAll[upper.tri(dataIntegrAll, diag=TRUE)] <- NA

  
  ##create pairs of drugs from benchmark 1 if same drug set from GMT and 0 otherwise
  benchDat[upper.tri(benchDat, diag=TRUE)] <- NA
  benchPairs <- melt(benchDat)
  benchPairs <- na.omit(benchPairs)
  colnames(benchPairs)[3] <- "bench"
  
  ## create pairs of drugs from combination
  integrPairs <- melt(dataIntegrAll)
  integrPairs <- na.omit(integrPairs)
  colnames(integrPairs)[3] <- "obs.combiall"
  
  res <- list(benchPairs=benchPairs,integrPairs=integrPairs)
  
  return(res)
  
}