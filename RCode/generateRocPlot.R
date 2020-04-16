###############################################################################################################
## Function generates ROC plots for every single layer separately and for the integration of all layers
## 
## input: 
##     allPairs: list of all pairs obtained for the benchmark, structure, sensitivity, perturbation, and integrative method 
##     d1Name: name of the dataset which is used for sensitivity information, e.g., "nci60" 
##     d2Name: currently must be set to "lincs"
##     benchName : name of the benchmark set, e.g., "stitch", "chembl"
## output: 
##     
##
## 
###############################################################################################################

generateRocPlot <- function(allPairs, d1Name, d2Name="lincs", benchName, datestr="" , dnfnames=list(),cut=FALSE,format="emf") {
  
  if(cut!=FALSE){
    predIntegr <- predPerf(allPairs$integrPairs$obs.combiall, allPairs$benchPairs$bench,cut = paste("cutpoint-integrate",d1Name,benchName,sep = "-"))
  }else{
    predIntegr <- predPerf(allPairs$integrPairs$obs.combiall, allPairs$benchPairs$bench,cut = cut)
  }
  predStrc <- predPerf(allPairs$strcPairs$obs.str, allPairs$benchPairs$bench)
  predSens <- predPerf(allPairs$sensPairs$obs.sens, allPairs$benchPairs$bench)
  predPert <- predPerf(allPairs$pertPairs$obs.pert, allPairs$benchPairs$bench)
  iorio <- predPerf(allPairs$iorioPairs$obs.iorio, allPairs$benchPairs$bench)
  iskar <- predPerf(allPairs$iskarPairs$obs.iskar, allPairs$benchPairs$bench)
  
  if (length(allPairs)>=9 & !is.null(allPairs$superPairs)) {
      predSuper <- predPerf(allPairs$superPairs$obs.superPred, allPairs$benchPairs$bench)
  }
  
  if (length(allPairs)>=9 & !is.null(allPairs$drugePairs)) {
      predDrugE <- predPerf(allPairs$drugePairs$obs.drugerank, allPairs$benchPairs$bench)
  } 
  if(length(allPairs)>9){
    NDNF <- length(allPairs)-9
    predDNFS <- list()
    for(i in 1:NDNF){
      
      if(cut!=FALSE){
        predDNFS[[i]] <- predPerf(allPairs[[9+i]][3], allPairs$benchPairs$bench,cut = paste("cutpoint",dnfnames[[i]],d1Name,benchName,sep = "-"))
      }else{
        predDNFS[[i]] <- predPerf(allPairs[[9+i]][3], allPairs$benchPairs$bench,cut = cut)
        
      }
      
      }
    
  }
  
  #Mod by chjy 20190124 created a dir to store the figures
  # rmchar <- "[\xb5]|[\n]|[,]|[;]|[:]|[-]|[+]|[*]|[%]|[$]|[#]|[{]|[}]|[[]|[]]|[|]|[\\^]|[/]|[\\]|[.]|[_]|[ ]"
  # datestr <- date()
  # mainDir <- getwd()
  #datestr <- gsub(rmchar, "",  datestr)
  subDir<-paste("Output/",datestr,sep = "")
  # options(someUniqueTag.mainDir = mainDir)
  # options(someUniqueTag.subDir = "subDir")
  # 
  # if(!file_test("-d", file.path(mainDir, subDir))){
  #   if(file_test("-f", file.path(mainDir, subDir))) {
  #     stop("Path can't be created because a file with that name already exists.")
  #   } else {
  #     dir.create(file.path(mainDir, subDir))
  #   }
  # }
  # 
  
  
  # changing params for the ROC plot - width, etc
  #filename = paste(getwd(), "/Output/", "ROC_", d1Name, "_", d2Name, "_", benchName, ".pdf", sep="")
  filename = paste("ROC_", d1Name, "_", d2Name, "_", benchName, ".",format, sep="")
  
  if(format=="pdf"){  pdf(file.path(mainDir, subDir, filename), width=2.5, height = 2.5)}
  if(format=="emf"){  emf(file.path(mainDir, subDir, filename), width=2.5, height = 2.5)}
  
  par(mgp = c(1.3,0.4, 0),mar=c(2.5,2.5,0.5,0.5), xaxs = "i", yaxs = "i", cex.axis=0.8, cex.lab=1.1,cex.main=1.1,bty = 'n')
  
  plot.new()
  
  #library(RColorBrewer)
  
  family <- c("black","#d7191c","#41ab5d","#2b83ba", "pink", "purple", "orange", "cyan","magenta","yellow")
  # plotting the ROC curve
  #plot(predIntegr$perf, col="black", lwd=2)
  # plot(predStrc$perf, col="#d7191c", lwd=2,add = TRUE)
  # plot(predSens$perf, col = "#41ab5d", lwd=2,add = TRUE)
  # plot(predPert$perf, col = "#2b83ba", lwd=2,add = TRUE)
  # plot(iorio$perf, col = "pink", lwd=2,add = TRUE)
  # plot(iskar$perf, col = "purple", lwd=2,add = TRUE)
  # 
  # if (length(allPairs)>=9 & !is.null(allPairs$superPairs)) {
  #     plot(predSuper$perf, col = "orange", lwd=2,add = TRUE)
  # }
  # if (length(allPairs)>=9 & !is.null(allPairs$drugePairs)) {
  #     plot(predDrugE$perf, col = "orange", lwd=2,add = TRUE)
  # }
  # if(length(allPairs)>9){
  #   NDNF <- length(allPairs)-9
  #   for(i in 1:NDNF){
  #     plot(predDNFS[[i]]$perf,lwd=2,add = TRUE)
  #     
  #   }
  #   
  # }
  
  rect(par("usr")[1],par("usr")[3],par("usr")[2],par("usr")[4],col = "#ebebeb")
  grid (lty = 1, col = "white")
  
  plot(predIntegr$perf, col=family[1], lwd=2)
  plot(predStrc$perf, col=family[2], lwd=2,add = TRUE)
  plot(predSens$perf, col = family[3], lwd=2,add = TRUE)
  plot(predPert$perf, col = family[4], lwd=2,add = TRUE)
  plot(iorio$perf, col = family[5], lwd=2,add = TRUE)
  plot(iskar$perf, col = family[6], lwd=2,add = TRUE)
  abline(0,1, col = "gray")
  
  if (length(allPairs)>=9 & !is.null(allPairs$superPairs)) {
    plot(predSuper$perf, col = family[7], lwd=2,add = TRUE)
  }
  if (length(allPairs)>=9 & !is.null(allPairs$drugePairs)) {
    plot(predDrugE$perf, col = family[7], lwd=2,add = TRUE)
  }
  if(length(allPairs)>9){
    NDNF <- length(allPairs)-9
    for(i in 1:NDNF){
      plot(predDNFS[[i]]$perf,col=family[7+i],lwd=2,add = TRUE)
      
    }
  }
  
  # Revised
  
  
  aucLegIntegr <- paste(c("RNCE = "), round(predIntegr$auc,3), sep="")
  aucLegStr <- paste(c("Structure = "), round(predStrc$auc,3),sep="")
  aucLegSen <- paste(c("Sensitivity = "), round(predSens$auc,3) , sep="")
  aucLegPer <- paste(c("Perturbation = "), round(predPert$auc,3), sep="")
  aucLegIorio <- paste(c("IorioPGX = "), round(iorio$auc,3), sep="")
  aucLegIskar <- paste(c("Iskar = "), round(iskar$auc,3), sep="")
  
  methods <- c("RNCE = ","Structure = ","Sensitivity = ","Perturbation = ","IorioPGX = ","Iskar = ")
  scores <- c(round(predIntegr$auc,3),round(predStrc$auc,3),round(predSens$auc,3),round(predPert$auc,3),round(iorio$auc,3),round(iskar$auc,3))
  
  
  if(length(allPairs)>9){
    aucLegDNFs <- list()
    for(i in 1:NDNF){
      if(length(dnfnames)<1){
        aucLegDNFs[[i]] <- paste(c("DNF methods = "), round(predDNFS[[i]]$auc,3), sep="")
        methods[[i+4]] <- "DNF methods = "
        scores[[i+4]] <-  round(predDNFS[[i]]$auc,3)
        
      }else{
        aucLegDNFs[[i]] <- paste(dnfnames[[i]]," = ", round(predDNFS[[i]]$auc,3), sep="")
        scores[[i+4]] <-  round(predDNFS[[i]]$auc,3)
        methods[[i+4]] <- paste(dnfnames[[i]]," = ", sep="")
        
      }
    }
    
  }
  
  score.frame <- data.frame("methods"= c(methods,"Rand = "),"scores"=c(scores,0.5),"color" = c(family[1:(4+length(dnfnames))],'grey'))
  
  score.frame <- score.frame[order(score.frame$scores,decreasing=TRUE),]
  
  legendtext <- paste(score.frame$methods,score.frame$scores,sep = "")

  if (length(allPairs)==9 & !is.null(allPairs$superPairs)) {
    
    
    
    aucLegSuper<- paste(c("SuperPred = "), round(predSuper$auc,3), sep="")
    legend(0.5,0.6,
           #c(aucLegIntegr, aucLegStr, aucLegSen, aucLegPer, aucLegIorio, aucLegIskar, aucLegSuper,rand), 
           c(aucLegIntegr, aucLegStr, aucLegSen, aucLegPer, aucLegIorio, aucLegIskar, aucLegSuper,rand), 
           bg="white",
           border="white", cex=0.6,
           box.col = "white",
           fill=c("black","#d7191c","#41ab5d","#2b83ba", "pink", "purple", "orange", "gray")
    )
  }
  if (length(allPairs)==9 & !is.null(allPairs$drugePairs)) {
    aucLegDrugE <- paste(c("DrugERank = "), round(predDrugE$auc,3), sep="")
    legend(0.35,0.6,
           c(aucLegIntegr, aucLegStr, aucLegSen, aucLegPer, aucLegIorio, aucLegIskar, aucLegDrugE,rand),  bg="white",
           border="white", cex=0.6,
           box.col = "white",#
           fill=c("black","#d7191c","#41ab5d","#2b83ba", "pink", "purple", "orange", "gray")
    )
  }
  
  if (length(allPairs)>=9 & !is.null(allPairs$superPairs)) {
    
    methods <- append(methods,"SuperPred = ",6)
    scores <- append(scores, round(predSuper$auc,3),6)
    
    score.frame <- data.frame("methods"= c(methods,"Rand ="),
                              "scores"=c(scores,0.5),
                              "color" = c(family[1:(length(allPairs)-1)],'grey'))
    
    score.frame <- score.frame[order(score.frame$scores,decreasing=TRUE),]
    legendtext <- paste(score.frame$methods,score.frame$scores,sep = "")
    
    
    aucLegSuper<- paste(c("SuperPred = "), round(predSuper$auc,3), sep="")
    # legend(0.35,0.6,
    #        #c(aucLegIntegr, aucLegStr, aucLegSen, aucLegPer, aucLegIorio, aucLegIskar, aucLegSuper,rand), 
    #        c(aucLegIntegr, aucLegStr, aucLegSen, aucLegPer, aucLegIorio, aucLegIskar, aucLegSuper,aucLegDNFs), 
    #        bg="white",
    #        border="white", cex=0.6,
    #        box.col = "white",#fill=c("black","#d7191c","#41ab5d","#2b83ba", "pink", "purple", "orange", "gray")
    #        fill=c(family[1:(length(allPairs)-1)])
    #        )
    
    legend(0.35,0.6,
           #c(aucLegIntegr, aucLegStr, aucLegSen, aucLegPer, aucLegIorio, aucLegIskar, aucLegSuper,rand), 
           legendtext, 
           bg="white",
           border="white", cex=0.6,
           box.col = "white",#fill=c("black","#d7191c","#41ab5d","#2b83ba", "pink", "purple", "orange", "gray")
           fill=as.character(score.frame$color)
    )
  }
  if (length(allPairs)>=9 & !is.null(allPairs$drugePairs)) {
    methods <- append(methods,"DrugERank = ",6)
    scores <- append(scores,round(predDrugE$auc,3),6)
    
    score.frame <- data.frame("methods"= c(methods,"Rand ="),
                              "scores"=c(scores,0.5),
                              "color" = c(family[1:(length(allPairs)-1)],'grey'))
    
    score.frame <- score.frame[order(score.frame$scores,decreasing=TRUE),]
    legendtext <- paste(score.frame$methods,score.frame$scores,sep = "")
    
    
    aucLegDrugE <- paste(c("DrugERank = "), round(predDrugE$auc,3), sep="")
    # legend(0.35,0.6,
    #        c(aucLegIntegr, aucLegStr, aucLegSen, aucLegPer, aucLegIorio, aucLegIskar, aucLegDrugE,aucLegDNFs), bg="white",
    #        border="white", cex=0.6,
    #        box.col = "white",#fill=c("black","#d7191c","#41ab5d","#2b83ba", "pink", "purple", "orange", "gray")
    #        fill=c(family[1:(length(allPairs)-1)])
    #        )
    
    legend(0.35,0.6,
           #c(aucLegIntegr, aucLegStr, aucLegSen, aucLegPer, aucLegIorio, aucLegIskar, aucLegSuper,rand), 
           legendtext, 
           bg="white",
           border="white", cex=0.6,
           box.col = "white",#fill=c("black","#d7191c","#41ab5d","#2b83ba", "pink", "purple", "orange", "gray")
           fill=as.character(score.frame$color)
    )
  }
  
  
  dev.off()
  
  return(append(predDNFS,predIntegr,0))
  
}

customRocPlot <- function(allPairs, d1Name, d2Name="lincs", benchName, datestr="" , dnfnames=list(),cut=FALSE,format="emf") {
  
  if(cut!=FALSE){
    predIntegr <- predPerf(allPairs$integrPairs$obs.combiall, allPairs$benchPairs$bench,cut = paste("cutpoint-integrate",d1Name,benchName,sep = "-"))
  }else{
    predIntegr <- predPerf(allPairs$integrPairs$obs.combiall, allPairs$benchPairs$bench,cut = cut)
  }
  predStrc <- predPerf(allPairs$strcPairs$obs.str, allPairs$benchPairs$bench)
  predSens <- predPerf(allPairs$sensPairs$obs.sens, allPairs$benchPairs$bench)
  predPert <- predPerf(allPairs$pertPairs$obs.pert, allPairs$benchPairs$bench)

  
  if(length(allPairs)>=9){
    NDNF <- length(allPairs)-9
    predDNFS <- list()
    for(i in 1:NDNF){
      
      if(cut!=FALSE){
        predDNFS[[i]] <- predPerf(allPairs[[9+i]][3], allPairs$benchPairs$bench,cut = paste("cutpoint",dnfnames[[i]],d1Name,benchName,sep = "-"))
      }else{
        predDNFS[[i]] <- predPerf(allPairs[[9+i]][3], allPairs$benchPairs$bench,cut = cut)
        
      }
      
    }
    
  }
  

  subDir<-paste("Output/",datestr,sep = "")

  filename = paste("ROC_", d1Name, "_", d2Name, "_", benchName, ".",format, sep="")
  
  if(format=="pdf"){  pdf(file.path(mainDir, subDir, filename), width=2.5, height = 2.5)}
  if(format=="emf"){  emf(file.path(mainDir, subDir, filename), width=2.5, height = 2.5)}
  
  par(mgp = c(1.3,0.4, 0),mar=c(2.5,2.5,0.5,0.5), xaxs = "i", yaxs = "i", cex.axis=0.8, cex.lab=1.1,cex.main=1.1,bty = 'n')
  
  plot.new()

  
  #library(RColorBrewer)
  
  
  family <- c("black","#d7191c","#41ab5d","#2b83ba", "pink", "purple", "orange", "cyan","magenta","yellow")

  
  rect(par("usr")[1],par("usr")[3],par("usr")[2],par("usr")[4],col = "#ebebeb")
  grid (lty = 1, col = "white")
  
  plot(predIntegr$perf, col=family[1], lwd=2)
  plot(predStrc$perf, col=family[2], lwd=2,add = TRUE)
  plot(predSens$perf, col = family[3], lwd=2,add = TRUE)
  plot(predPert$perf, col = family[4], lwd=2,add = TRUE)

  
  if(length(allPairs)>9){
    NDNF <- length(allPairs)-9
    for(i in 1:NDNF){
      plot(predDNFS[[i]]$perf,col=family[4+i],lwd=2,add = TRUE)
      
    }
  }
  
  aucLegIntegr <- paste(c("RNCE = "), round(predIntegr$auc,3), sep="")
  aucLegStr <- paste(c("Structure = "), round(predStrc$auc,3),sep="")
  aucLegSen <- paste(c("Sensitivity = "), round(predSens$auc,3) , sep="")
  aucLegPer <- paste(c("Perturbation = "), round(predPert$auc,3), sep="")
  
  methods <- c("RNCE = ","Structure = ","Sensitivity = ","Perturbation = ")
  scores <- c(round(predIntegr$auc,3),round(predStrc$auc,3),round(predSens$auc,3),round(predPert$auc,3))

  if(length(allPairs)>9){
    aucLegDNFs <- list()
    for(i in 1:NDNF){
      if(length(dnfnames)<1){
        aucLegDNFs[[i]] <- paste(c("DNF methods = "), round(predDNFS[[i]]$auc,3), sep="")
        methods[[i+4]] <- "DNF methods = "
        scores[[i+4]] <-  round(predDNFS[[i]]$auc,3)
        
      }else{
        aucLegDNFs[[i]] <- paste(dnfnames[[i]]," = ", round(predDNFS[[i]]$auc,3), sep="")
        scores[[i+4]] <-  round(predDNFS[[i]]$auc,3)
        methods[[i+4]] <- paste(dnfnames[[i]]," = ", sep="")
        
      }
    }
    
  }
  
  score.frame <- data.frame("methods"= c(methods,"Rand = "),"scores"=c(scores,0.5),"color" = c(family[1:(4+length(dnfnames))],'grey'))
  
  score.frame <- score.frame[order(score.frame$scores,decreasing=TRUE),]
  
  legendtext <- paste(score.frame$methods,score.frame$scores,sep = "")

  abline(0,1, col = "gray")
  
  if (length(allPairs)>=9) {
    # legend(0.35,0.5,
    #        c(aucLegIntegr, aucLegStr, aucLegSen, aucLegPer,aucLegDNFs), bg="white",
    #        border="white", cex=0.6,
    #        box.col = "white",#fill=c("black","#d7191c","#41ab5d","#2b83ba", "pink", "purple", "orange", "gray")
    #        fill=c(family[1:(4+length(dnfnames))])
    # )
    
    legend(0.35,0.55,
           #c(aucLegIntegr, aucLegStr, aucLegSen, aucLegPer,aucLegDNFs), 
           legendtext,
           bg="white",
           border="white", cex=0.6,
           box.col = "white",#fill=c("black","#d7191c","#41ab5d","#2b83ba", "pink", "purple", "orange", "gray")
           fill=as.character(score.frame$color)
    )
  }
  

  
  
  
  dev.off()
  
  return(append(predDNFS,predIntegr,0))
  
}