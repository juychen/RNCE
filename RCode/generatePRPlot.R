


generatePRPlot <- function(allPairs, d1Name, d2Name="lincs", benchName,datestr="", dnfnames=list(),format="emf") {
  
  
  predIntegr <- predPerf(allPairs$integrPairs$obs.combiall, allPairs$benchPairs$bench, "PR")
  predStrc <- predPerf(allPairs$strcPairs$obs.str, allPairs$benchPairs$bench, "PR")
  predSens <- predPerf(allPairs$sensPairs$obs.sens, allPairs$benchPairs$bench, "PR")
  predPert <- predPerf(allPairs$pertPairs$obs.pert, allPairs$benchPairs$bench, "PR")
  iorio <- predPerf(allPairs$iorioPairs$obs.iorio, allPairs$benchPairs$bench, "PR")
  iskar <- predPerf(allPairs$iskarPairs$obs.iskar, allPairs$benchPairs$bench, "PR")
  
  if (length(allPairs)>=9 & !is.null(allPairs$superPairs)) {
      super <- predPerf(allPairs$superPairs$obs.superPred, allPairs$benchPairs$bench, "PR")
  }
  if (length(allPairs)>=9 & !is.null(allPairs$drugePairs)) {
      druge <- predPerf(allPairs$drugePairs$obs.drugerank, allPairs$benchPairs$bench, "PR")
  }
  
  if(length(allPairs)>9){
    NDNF <- length(allPairs)-9
    predDNFS <- list()
    for(i in 1:NDNF){
      thisPair<-as.matrix(allPairs[[9+i]][3])
      predDNFS[[i]] <- predPerf(thisPair, allPairs$benchPairs$bench,"PR")
    }
    
  }
  
  #Mod by chjy 20190124 created a dir to store the figures
  # rmchar <- "[\xb5]|[\n]|[,]|[;]|[:]|[-]|[+]|[*]|[%]|[$]|[#]|[{]|[}]|[[]|[]]|[|]|[\\^]|[/]|[\\]|[.]|[_]|[ ]"
  # datestr <- date()
  # mainDir <- getwd()
  # datestr <- gsub(rmchar, "",  datestr)
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
  

  
  # changing params for the ROC plot - width, etc
  #filename = paste(getwd(), "/Output/", "NEW_PR_", d1Name, "_", d2Name, "_", benchName, ".pdf", sep="")
  #pdf(filename, width=2.5, height = 2.55)
  filename = paste("NEW_PR_", d1Name, "_", d2Name, "_", benchName, ".",format, sep="")
  if(format=="pdf"){pdf(file.path(mainDir, subDir, filename), width=2.5, height = 2.55)}
  if(format=="emf"){emf(file.path(mainDir, subDir, filename), width=2.5, height = 2.55)}
  par(mgp = c(1.3,0.4, 0),mar=c(2.5,2.5,0.5,0.5), xaxs = "i", yaxs = "i", cex.axis=0.8, cex.lab=1.1,cex.main=1.1,bty = 'n')
  # plotting the ROC curve
  
  plot.new()
  
  family <- c("black","#d7191c","#41ab5d","#2b83ba", "pink", "purple", "orange", "cyan","magenta","yellow")
  
  
  
  rect(par("usr")[1],par("usr")[3],par("usr")[2],par("usr")[4],col = "#ebebeb")
  grid (lty = 1, col = "white")
  
  aucLegIntegr <- paste( round(predIntegr$auc.integral,3), sep="")
  

  plot(predIntegr, col=family[1], lwd=2 , main="", ann=FALSE)
  plot(predStrc, col=family[2], lwd=2,add = TRUE)
  plot(predSens, col = family[3], lwd=2,add = TRUE)
  plot(predPert, col = family[4], lwd=2,add = TRUE)
  plot(iorio, col = family[5], lwd=2,add = TRUE)
  plot(iskar, col = family[6], lwd=2,add = TRUE)
  plot(predIntegr$rand, col = "gray", lwd=2,add = TRUE)
  
  if (length(allPairs)>=9 & !is.null(allPairs$superPairs)) {
      plot(super, col = family[7], lwd=2,add = TRUE)
  }
  if (length(allPairs)>=9 & !is.null(allPairs$drugePairs)) {
      plot(druge, col = family[7], lwd=2,add = TRUE)
  }
  if(length(allPairs)>9){
    NDNF <- length(allPairs)-9
    for(i in 1:NDNF){
      plot(predDNFS[[i]],col=family[7+i],lwd=2,add = TRUE)
      
    }
  }
  
  # Label the x and y axes with dark green text
  title(xlab="Recall")
  title(ylab="Precision")

  #aucLegIntegr <- paste( round(predIntegr$auc.integral,3), sep="")
  aucLegStr <- paste( round(predStrc$auc.integral,3),sep="")
  aucLegSen <- paste( round(predSens$auc.integral,3) , sep="")
  aucLegPer <- paste( round(predPert$auc.integral,3), sep="")
  aucLegIorio <- paste( round(iorio$auc.integral,3), sep="")
  aucLegIskar <- paste( round(iskar$auc.integral,3), sep="")
  
  if(length(allPairs)>9){
    aucLegDNFs <- list()
    for(i in 1:NDNF){
      if(length(dnfnames)<1){
        aucLegDNFs[[i]] <- paste( round(predDNFS[[i]]$auc.integral,3), sep="")
        
      }else{
        aucLegDNFs[[i]] <- paste(round(predDNFS[[i]]$auc.integral,3), sep="")
        
      }
    }
  }
  rand <- paste(round(predIntegr$rand$auc.integral,3), sep="")
  
  
  if (length(allPairs)==9 & !is.null(allPairs$superPairs)) {
    aucLegSuper<- paste(round(super$auc.integral,3), sep="")
    legend(0.75,1,c(aucLegIntegr, aucLegStr, aucLegSen, aucLegPer, aucLegIorio, aucLegIskar, aucLegSuper, rand), border="white", cex=0.6,
           box.col = "white",bg="white",
           fill=c("black","#d7191c","#41ab5d","#2b83ba", "pink", "purple", "orange", "gray")
    )
  }
  
  if (length(allPairs)==9 & !is.null(allPairs$drugePairs)) {
    aucLegDruge <- paste( round(druge$auc.integral,3), sep="")
    legend(0.75,1,c(aucLegIntegr, aucLegStr, aucLegSen, aucLegPer, aucLegIorio, aucLegIskar, aucLegDruge, rand), border="white", cex=0.6,
           box.col = "white",bg="white",
           fill=c("black","#d7191c","#41ab5d","#2b83ba", "pink", "purple", "orange", "gray")
    )
  }
  if (length(allPairs)>9 & !is.null(allPairs$superPairs)) {
    aucLegSuper<- paste(round(super$auc.integral,3), sep="")
    legend(0.75,1,c(aucLegIntegr, aucLegStr, aucLegSen, aucLegPer, aucLegIorio, aucLegIskar, aucLegSuper,aucLegDNFs, rand), border="white", cex=0.6,
           box.col = "white",bg="white",
           fill=c(family[1:(length(allPairs)-2)], "gray")
           #fill=c("black","#d7191c","#41ab5d","#2b83ba", "pink", "purple", "orange", "gray")
    )
  }
  
  if (length(allPairs)>9 & !is.null(allPairs$drugePairs)) {
    aucLegDruge <- paste( round(druge$auc.integral,3), sep="")
    legend(0.75,1,c(aucLegIntegr, aucLegStr, aucLegSen, aucLegPer, aucLegIorio, aucLegIskar, aucLegDruge,aucLegDNFs, rand), border="white", cex=0.6,
           box.col = "white",bg="white",
           fill=c(family[1:(length(allPairs)-2)], "gray")
           #fill=c("black","#d7191c","#41ab5d","#2b83ba", "pink", "purple", "orange", "gray")
    )
  }
  dev.off()
  
  aucLegIntegr <- paste(c("RNCE = "), round(predIntegr$auc.integral,3), sep="")
  aucLegStr <- paste(c("Structure = "), round(predStrc$auc.integral,3),sep="")
  aucLegSen <- paste(c("Sensitivity = "), round(predSens$auc.integral,3) , sep="")
  aucLegPer <- paste(c("Perturbation = "), round(predPert$auc.integral,3), sep="")
  aucLegIorio <- paste(c("IorioPGX = "), round(iorio$auc.integral,3), sep="")
  aucLegIskar <- paste(c("Iskar = "), round(iskar$auc.integral,3), sep="")
  
  if(length(allPairs)>9){
    aucLegDNFs <- list()
    for(i in 1:NDNF){
      if(length(dnfnames)<1){
        aucLegDNFs[[i]] <- paste(c("DNF methods = "), round(predDNFS[[i]]$auc.integral,3), sep="")
        
      }else{
        aucLegDNFs[[i]] <- paste(dnfnames[[i]]," = ", round(predDNFS[[i]]$auc.integral,3), sep="")
        
      }
    }
    
  }
  
  
  filename = paste("NEW_PR_LEGEND", d1Name, "_", d2Name, "_", benchName, ".",format, sep="")
  
  
  if(format=="pdf"){pdf(file.path(mainDir, subDir, filename), width=2.5, height = 2.55)}
  if(format=="emf"){emf(file.path(mainDir, subDir, filename), width=2.5, height = 2.55)}
  par(mgp = c(1.3,0.4, 0),mar=c(2.5,2.5,0.5,0.5), xaxs = "i", yaxs = "i", cex.axis=0.8, cex.lab=1.1,cex.main=0.001,bty = 'n')
  # plotting the ROC curve
  
  plot.new()
  
  rand <- paste(c("Rand = "),round(predIntegr$rand$auc.integral,3), sep="")
  
  if (length(allPairs)==9 & !is.null(allPairs$superPairs)) {
    aucLegSuper<- paste(c("SuperPred = "), round(super$auc.integral,3), sep="")
    legend(0.35,1,c(aucLegIntegr, aucLegStr, aucLegSen, aucLegPer, aucLegIorio, aucLegIskar, aucLegSuper, rand), border="white", cex=0.6,
           box.col = "white",bg="white",
           fill=c("black","#d7191c","#41ab5d","#2b83ba", "pink", "purple", "orange", "gray")
    )
  }
  
  if (length(allPairs)==9 & !is.null(allPairs$drugePairs)) {
    aucLegDruge <- paste(c("DrugERank = "), round(druge$auc.integral,3), sep="")
    legend(0.35,1,c(aucLegIntegr, aucLegStr, aucLegSen, aucLegPer, aucLegIorio, aucLegIskar, aucLegDruge, rand), border="white", cex=0.6,
           box.col = "white",bg="white",
           fill=c("black","#d7191c","#41ab5d","#2b83ba", "pink", "purple", "orange", "gray")
    )
  }
  if (length(allPairs)>9 & !is.null(allPairs$superPairs)) {
    aucLegSuper <- paste(c("SuperPred = "), round(super$auc.integral,3), sep="")
    legend(0.35,1,c(aucLegIntegr, aucLegStr, aucLegSen, aucLegPer, aucLegIorio, aucLegIskar, aucLegSuper,aucLegDNFs, rand), border="white", cex=0.6,
           box.col = "white",bg="white",
           fill=c(family[1:(length(allPairs)-2)], "gray")
           #fill=c("black","#d7191c","#41ab5d","#2b83ba", "pink", "purple", "orange", "gray")
           )
  }
  
  if (length(allPairs)>9 & !is.null(allPairs$drugePairs)) {
    aucLegSuper<- paste(c("DrugERank = "), round(druge$auc.integral,3), sep="")
    legend(0.35,1,c(aucLegIntegr, aucLegStr, aucLegSen, aucLegPer, aucLegIorio, aucLegIskar, aucLegSuper,aucLegDNFs, rand), border="white", cex=0.6,
           box.col = "white",bg="white",
           fill=c(family[1:(length(allPairs)-2)], "gray")
           #fill=c("black","#d7191c","#41ab5d","#2b83ba", "pink", "purple", "orange", "gray")
           )
  }
  
  dev.off()
  
  
}

customPRPlot <- function(allPairs, d1Name, d2Name="lincs", benchName,datestr="", dnfnames=list(),format="emf") {
  
  
  predIntegr <- predPerf(allPairs$integrPairs$obs.combiall, allPairs$benchPairs$bench, "PR")
  predStrc <- predPerf(allPairs$strcPairs$obs.str, allPairs$benchPairs$bench, "PR")
  predSens <- predPerf(allPairs$sensPairs$obs.sens, allPairs$benchPairs$bench, "PR")
  predPert <- predPerf(allPairs$pertPairs$obs.pert, allPairs$benchPairs$bench, "PR")
  
  
  if(length(allPairs)>9){
    NDNF <- length(allPairs)-9
    predDNFS <- list()
    for(i in 1:NDNF){
      thisPair<-as.matrix(allPairs[[9+i]][3])
      predDNFS[[i]] <- predPerf(thisPair, allPairs$benchPairs$bench,"PR")
    }
    
  }
  

  subDir<-paste("Output/",datestr,sep = "")
  
  
  filename = paste("NEW_PR_", d1Name, "_", d2Name, "_", benchName, ".",format, sep="")
  if(format=="pdf"){pdf(file.path(mainDir, subDir, filename), width=2.5, height = 2.5)}
  if(format=="emf"){emf(file.path(mainDir, subDir, filename), width=2.5, height = 2.5)}
  par(mgp = c(1.3,0.4, 0),mar=c(2.5,2.5,0.5,0.5), xaxs = "i", yaxs = "i", cex.axis=0.8, cex.lab=1.1,cex.main=0.001,bty = 'n')
  # plotting the ROC curve
  
  plot.new()
  family <- c("black","#d7191c","#41ab5d","#2b83ba", "pink", "purple", "orange", "cyan","magenta","yellow")
  
  
  
  
  rect(par("usr")[1],par("usr")[3],par("usr")[2],par("usr")[4],col = "#ebebeb")
  grid (lty = 1, col = "white")
  
  
  aucLegIntegr <- paste( round(predIntegr$auc.integral,3), sep="")
  
  plot(predIntegr, col=family[1], lwd=2,main="",sub="", ann=FALSE)
  plot(predStrc, col=family[2], lwd=2,add = TRUE)
  plot(predSens, col = family[3], lwd=2,add = TRUE)
  plot(predPert, col = family[4], lwd=2,add = TRUE)
  plot(predIntegr$rand, col = "gray", lwd=2,add = TRUE)
  if(length(allPairs)>9){
    NDNF <- length(allPairs)-9
    for(i in 1:NDNF){
      plot(predDNFS[[i]],col=family[4+i],lwd=2,add = TRUE)
      
    }
  }
  
  # Label the x and y axes with dark green text
  title(xlab="Recall")
  title(ylab="Precision")
  
  
  #aucLegIntegr <- paste(round(predIntegr$auc.integral,3), sep="")
  aucLegStr <- paste(round(predStrc$auc.integral,3),sep="")
  aucLegSen <- paste(round(predSens$auc.integral,3) , sep="")
  aucLegPer <- paste(round(predPert$auc.integral,3), sep="")


  if(length(allPairs)>9){
    aucLegDNFs <- list()
    for(i in 1:NDNF){
      if(length(dnfnames)<1){
        aucLegDNFs[[i]] <- paste(round(predDNFS[[i]]$auc.integral,3), sep="")
        
      }else{
        aucLegDNFs[[i]] <- paste(round(predDNFS[[i]]$auc.integral,3), sep="")
        
      }
    }
  }
  
  if (length(allPairs)>=9) {
    par(mgp = c(1.3,0.4, 0),mar=c(2.5,2.5,0.5,0.5), xaxs = "i", yaxs = "i", cex.axis=0.8, cex.lab=1.1,cex.main=0.001,bty = 'n')
    
    plot.new()
    legend(0.75,1,
           c(aucLegIntegr, aucLegStr, aucLegSen, aucLegPer,aucLegDNFs), bg="white",
           border="white", cex=0.6,
           box.col = "white",#fill=c("black","#d7191c","#41ab5d","#2b83ba", "pink", "purple", "orange", "gray")
           fill=c(family[1:(4+length(dnfnames))])
    )
  }
  
  
  dev.off()
  
  
  aucLegIntegr <- paste(c("RNCE = "), round(predIntegr$auc.integral,3), sep="")
  aucLegStr <- paste(c("Structure = "), round(predStrc$auc.integral,3),sep="")
  aucLegSen <- paste(c("Sensitivity = "), round(predSens$auc.integral,3) , sep="")
  aucLegPer <- paste(c("Perturbation = "), round(predPert$auc.integral,3), sep="")
  
  
  
  if(length(allPairs)>9){
    aucLegDNFs <- list()
    for(i in 1:NDNF){
      if(length(dnfnames)<1){
        aucLegDNFs[[i]] <- paste(c("DNF methods = "), round(predDNFS[[i]]$auc.integral,3), sep="")
        
      }else{
        aucLegDNFs[[i]] <- paste(dnfnames[[i]]," = ", round(predDNFS[[i]]$auc.integral,3), sep="")
        
      }
    }
    
  }
  
  
  rand <- paste(c("Rand = "), round(predIntegr$rand$auc.integral,3), sep="")
  
  filename = paste("NEW_PR_LEGEND", d1Name, "_", d2Name, "_", benchName, ".",format, sep="")
  if(format=="pdf"){pdf(file.path(mainDir, subDir, filename), width=2.5, height = 2.5)}
  if(format=="emf"){emf(file.path(mainDir, subDir, filename), width=2.5, height = 2.5)}
  if (length(allPairs)>=9) {
    par(mgp = c(1.3,0.4, 0),mar=c(2.5,2.5,0.5,0.5), xaxs = "i", yaxs = "i", cex.axis=0.8, cex.lab=1.1,cex.main=0.001,bty = 'n')
    
    plot.new()
    legend(0.35,1,
           c(aucLegIntegr, aucLegStr, aucLegSen, aucLegPer,aucLegDNFs), bg="white",
           border="white", cex=0.6,
           box.col = "white",#fill=c("black","#d7191c","#41ab5d","#2b83ba", "pink", "purple", "orange", "gray")
           fill=c(family[1:(4+length(dnfnames))])
    )
  }
  
  
  
  dev.off()
  
  
}
