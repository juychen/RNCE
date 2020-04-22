library(dplyr)
library(RecordLinkage)

matchstring <- function(string, stringVector){
  
  distance <- levenshteinSim(string, stringVector);
  matcn<-stringVector[distance == max(distance)]
  
  #matcn <-paste(matcn,collapse = "; ")
  return(list(score= max(distance),bestmatch=matcn[1]))
}
matchFDAdrugsbysim <- function(df.druginfo,bench='ctrp',thresdhold=1){
  
  if(bench=="nci60"){
    threshold<-0.85
  }
  
  
  df.fda <-read.csv("Data/Products.txt",sep = "\t")
  drugname.fda<-as.character(unique(df.fda$DrugName))
  
  
  drug.names <- as.character(df.druginfo$MOLECULE_NAME)
  temp.bestmatch <- c()
  temp.score <-c()
  
  for (d in drug.names){
    temp.result<-matchstring(d,drugname.fda)
    temp.bestmatch <- append(temp.bestmatch,temp.result$bestmatch)
    temp.score <- append(temp.score,temp.result$score)
    
  }
  
  
  df.match <- data.frame(drug=drug.names,match=temp.bestmatch,score=temp.score)
  df.match <- df.match[order(df.match$score,decreasing = TRUE),]
  
  df.result <- df.match[,c(1,3)] 
  
  df.result[df.result$score>=thresdhold,2] <- TRUE
  df.result[df.result$score<thresdhold,2] <- FALSE
  
  colnames(df.result) <- c("MOLECULE_NAME","FDA_APPROVED")
  df.result$FDA_APPROVED <- as.boolean(df.result$FDA_APPROVED)
  
  return(df.result)
  
}

bestmatch.ctrp <- c() 
scores.ctrp <- c()

df.fda <-read.csv("Data/Products.txt",sep = "\t")
drugname.fda<-as.character(unique(df.fda$DrugName))

df.dinfo.ctrp <-read.csv("Data/druginfo_ctrp.csv")
drugname.ctrp <- as.character(df.dinfo.ctrp$MOLECULE_NAME)


df.dinfo.nci60 <-read.csv("Data/druginfo_nci60.csv")
drugname.nci60 <- as.character(df.dinfo.nci60$MOLECULE_NAME)

flag.ctrp <- matchFDAdrugsbysim(df.dinfo.ctrp,bench = "ctrp")
flag.nci60 <- matchFDAdrugsbysim(df.dinfo.nci60,bench = "nci60")

ctrp.res<-merge(df.dinfo.ctrp,flag.ctrp)
nci60.res<-merge(df.dinfo.nci60,flag.nci60)

write.csv(ctrp.res[,c(1,3,4,5)],"Data/ctrpfulldruginfo.csv")
write.csv(nci60.res[,c(1,3,4,5)],"Data/nci60fulldruginfo.csv")
