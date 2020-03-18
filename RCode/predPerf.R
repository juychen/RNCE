

predPerf <- function(layerPair, benchPair, plotType="ROC",cut=FALSE,
                     criterion= c("CB","MaxSumNPVPPV","MCT","ValuePPV","MinPvalue","MaxSpSe")
                     ) {
  #Max prosPSE youden roc01 kappa too much FP, ValuePPV PROC01 too much FN,efficicient MCT same as cb, Minimax both high
  
  
  if (plotType == "ROC") {
    pred <- prediction(layerPair, benchPair)
    perf <- performance(pred,"tpr","fpr")
    auc <- performance(pred,"auc")
    f <- performance(pred,"f")
    
    f <- performance(pred,"f")
    
    
    auc <- unlist(slot(auc, "y.values"))
    
    df_bench <- data.frame(layerPair,benchPair)
    colnames(df_bench) <- c("layerPair", "benchPair")
    
    if(cut!=FALSE){
    optimal.cut<- optimal.cutpoints(X="layerPair",status = "benchPair",tag.healthy = 0,methods = criterion,data=df_bench)
    print(summary(optimal.cut))
    print(cut)
    
    sink(paste("Output/",cut,".txt",sep = ""))
    print(sum(benchPair))
    print(summary(optimal.cut))
    sink()
    
    cutoffs<-c()
    
    for(cri in criterion){
      coff<-optimal.cut[[cri]]$Global$optimal.cutoff$cutoff
      
      if(is.null(coff)){
        coff<-optimal.cut[[cri]][[1]]$optimal.cutoff$cutoff
      }
      cutoffs<-append(coff,cutoffs)
    }
    
    names(cutoffs)<-criterion
    
    }
    
    else{
      cutoffs <- c()
    }
    return(list(auc=auc, pred=pred, perf=perf, cutoff = cutoffs))
    
  }
  
  ## NOTE  
  if (plotType == "PR") {
       PR <- PRROC::pr.curve(scores.class0=layerPair, weights.class0=benchPair, curve=TRUE, 
                             max.compute=TRUE, min.compute=TRUE, rand.compute=TRUE)
       return(PR)
  }
  

}
