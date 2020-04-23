library(netcom)  
library(apcluster)
library(xlsx)
library(igraph)
library(mclust)
library(netbiov)
library(devEMF)
library(ggplot2)
library(RColorBrewer)

source('RCode/utilities.R')
source('RCode/affinityClustering.R')

# Deinfine unexpected charaters
badchars <- "[\xb5]|[\n]|[,]|[;]|[:]|[-]|[+]|[*]|[%]|[$]|[#]|[{]|[}]|[[]|[]]|[|]|[\\^]|[/]|[\\]|[.]|[_]|[ ]"


# Load data files
load("Output/dataBench-ctrpv2.RData") ##load "drug-target bench for ctrpv2" here ...
dtgbench.ctrp <- dataBench
load("Output/ATCBench-ctrpv2.RData") ##load "nci60" results here ... 
ATCbench.ctrp <- dataBench3
load("Output/dataBench-nci60.RData") ##load "iskar" results here ... (see iskar.R)
dtgbench.nci60 <- dataBench1.5
load("Output/ATCBench-nci60.RData") ##load "iskar" results here ... (see iskar.R)
ATCbench.nci60 <- dataBench3

name.integrate.ctrp<- "Data/integrationSimilarityMat-ctrpv2.RData"
name.integrate.nc60<- "Data/integrationSimilarityMat-nci60.RData"

name.integrate.ctrp.dnf<- "Data/ctrpv2-DNFIntegrated.RData"
name.integrate.nc60.dnf<- "Data/nci60-DNFIntegrated.RData"

name.save.ctrp<- gsub(badchars,"",name.integrate.ctrp)
name.save.nc60<- gsub(badchars,"",name.integrate.nc60)

name.save.ctrp.dnf<- gsub(badchars,"",name.integrate.ctrp.dnf)
name.save.nc60.dnf<- gsub(badchars,"",name.integrate.nc60.dnf)

# Funtion to assign names to the matrix
submatbyname<-function(mat,name){
  rmat <- mat[name,name]
  colnames(rmat)<-rownames(rmat) <- name
  return(rmat)
}  

# Get our data and comparing data
load(name.integrate.ctrp)
ctrp <- integrtStrctSensPert
diag(ctrp) <- 1

load(name.integrate.nc60)
nc60 <- integrtStrctSensPert
diag(nc60) <- 1

load(name.integrate.ctrp.dnf)
ctrp.dnf <- integrtStrctSensPert
diag(ctrp.dnf) <- 1

load(name.integrate.nc60.dnf)
nc60.dnf <- integrtStrctSensPert
diag(nc60.dnf) <- 1


# Asssign names to the similarity matrix
name.ctrp.atc <- rownames(ATCbench.ctrp)
name.ctrp.dtg <- rownames(dtgbench.ctrp)
name.nci60.atc <- rownames(ATCbench.nci60)
name.nci60.dtg <- rownames(dtgbench.nci60)


# Select the benchmarking subsets ON CTRP AND NCI60
ctrp.sub.dtg <- submatbyname(ctrp,name.ctrp.dtg)
ctrp.sub.atc <- submatbyname(ctrp,name.ctrp.atc)

nc60.sub.dtg <- submatbyname(nc60,name.nci60.dtg)
nc60.sub.atc <- submatbyname(nc60,name.nci60.atc)

ctrp.subdnf.dtg <- submatbyname(ctrp.dnf,name.ctrp.dtg)
ctrp.subdnf.atc <- submatbyname(ctrp.dnf,name.ctrp.atc)

nc60.subdnf.dtg <- submatbyname(nc60.dnf,name.nci60.dtg)
nc60.subdnf.atc <- submatbyname(nc60.dnf,name.nci60.atc)


# Apply the louvain on the benchmark ctrp dataset
loulabel.dtgbench.ctrp<-affinClustering(dtgbench.ctrp,method="louvain")
loulabel.atcbench.ctrp<-affinClustering(ATCbench.ctrp,method="louvain")


# Apply the louvain on the benchmark nci60 dataset
loulabel.dtgbench.nci60<-affinClustering(dtgbench.nci60,method="louvain")
loulabel.atcbench.nci60<-affinClustering(ATCbench.nci60,method="louvain")

# Apply the spectral clustering on the four datasets
speclabel.ctrp.dtg<-affinClustering(ctrp.sub.dtg,method="spectral", K=max(loulabel.dtgbench.ctrp))
speclabel.ctrp.atc<-affinClustering(ctrp.sub.atc,method="spectral", K=max(loulabel.atcbench.ctrp))
speclabel.nc60.dtg<-affinClustering(nc60.sub.dtg,method="spectral", K=max(loulabel.dtgbench.nci60))
speclabel.nc60.atc<-affinClustering(nc60.sub.atc,method="spectral", K=max(loulabel.atcbench.nci60))

# apcomb.ctrp<-affinClustering(ctrp,method="spectral", K=35)
# apcomb.nc60<-affinClustering(nc60,method="spectral", K=33)

ctrpcom<-load("Data/CTRPv2commSpec.Rdata")
apcomb.ctrp<-get(ctrpcom)
nc60com<-load("Data/NCI60commSpec.RData")
apcomb.nc60<-get(nc60com)


commualign <- function(integrate,bench,apcomb,threshold = 2,benchname="target",savename="",overwrite=TRUE,novel=TRUE,plot=FALSE){
  # Write the community infor and test alignment scores
  score.align<-c(0)
  score.count<-c(0)
  size.intersect<-c(0)
  
  names.novel<-list()
  saveFileName <- paste("output/commnuity",benchname,savename,".xlsx",sep='_')
  
  comm.targets <- list()
  
  # Remove file if exsist and want to overwrite, then manipulate exsisting result
  if(file.exists(saveFileName) & (overwrite == TRUE) ){
    file.remove(saveFileName)
  }
  
  # Deal with the clustering results
  
  clusterNames<-labeltoCluster(apcomb)$nameList
  exeamplars<-exeamplarCluster(clusterNames,integrate)
  # Comments by jy in 20200116 to change apcomb to the labels
  #for(i in 1:length(apcomb@clusters)){
  names(clusterNames)<-exeamplars
  
  sortedCluster <- clusterNames[order(names(clusterNames))]
  colnumber <- (c(1:length(sortedCluster)))
  names(colnumber) <-names(sortedCluster)
  
  
  # Load drug targets information
  
  if(benchname=="target"){
    
    ctrpDrTargs <- read.csv("./Data/CTRPv2_drugtarget.csv", stringsAsFactors = FALSE) # 481 drugs x 3 descriptions
    ctrpDrTargs$compound_name <- toupper(ctrpDrTargs$compound_name)
    ctrpDrTargs$compound_name <- gsub(badchars,"",ctrpDrTargs$compound_name)
  }else{
    
    badchars <- "[\xb5]|[\n]|[,]|[;]|[:]|[-]|[+]|[*]|[%]|[$]|[#]|[{]|[}]|[[]|[]]|[|]|[\\^]|[/]|[\\]|[.]|[_]|[ ]"
    chembl_ATC <- read.delim("Data/chembl_drugtargets-16_5-10-02.txt", stringsAsFactor=F, na.strings=c("", "NA")) #10506 drugs!
    sep = " \\("
    chembl_ATC$MOLECULE_NAME <- toupper(unlist(lapply(chembl_ATC$MOLECULE_NAME, function(x) strsplit(x, sep)[[1]][1])))
    
  }
  
  # Load drug info for ctrp and nci60
  
  if(dim(integrate)[1]==239){
    df.druginfo<-read.csv("Data/ctrpfulldruginfo.csv")
    df.druginfo<-df.druginfo[,c(2,3,4,5)]
  }else{
    df.druginfo<-read.csv("Data/nci60fulldruginfo.csv")
    df.druginfo<-df.druginfo[,c(2,3,4,5)]


  }
  
  drugs.FDAap <- unlist(df.druginfo[df.druginfo$FDA_APPROVED == TRUE,1])
  
  # Set plot colors
  
  if(plot==TRUE){
    ncolors <- max(apcomb)
    qual_col_pals <- brewer.pal.info[brewer.pal.info$category == 'qual',]
    col_vector <- unlist(mapply(brewer.pal, qual_col_pals$maxcolors, rownames(qual_col_pals)))
    #pie(rep(1,ncolors), col=sample(col_vector, ncolors))
    mycols<-col_vector[1:ncolors]
    
  }
  
  # Process clusters and plot clusters
  for(i in 1:max(apcomb)){
    
    names.taxonomy<-names(apcomb[apcomb==i]) 
    
    names.bench<-rownames(bench)
    it<-intersect(names.bench,names.taxonomy)
    names.novel[[i]] <- setdiff(names.taxonomy,it)
    if(length(it)>=threshold){
      
      if(plot==TRUE){
        
        col.node <- mycols[colnumber[exeamplars[[i]]]]
        if ( ! file.exists("Output/communities")) {
          dir.create("Output/communities")
        }
        
        emf(paste("output/communities/community",savename,i,".emf",sep='_'),height=1.5,width=1.9)
        par(mar=c(0,0,0,0))
        
        
        adj.integarte<-0.01/integrate[names.taxonomy,names.taxonomy]
        diag(adj.integarte)<-0
        
        color.vlabs <- rep("black",length(it))
        
        color.vlabs[it %in% drugs.FDAap]<-"red"
        
        
        g.commun <- graph.adjacency(adj.integarte,weighted=T,mode="undirected")
        # plot( mst(g.commun),layout=layout_with_gem,
        #      vertex.size=7,vertex.color = "lightgrey",vertex.label.cex=0.6, vertex.label.dist=2)
        mst.plot( (g.commun),v.lab.col=color.vlabs,vertex.color =col.node,
                  v.size=10,mst.edge.col="grey",colors="gray97",layout=layout_with_gem,bg=NA,v.lab=TRUE,v.lab.cex=0.65,lab.dist=0.5)
        
        dev.off()
        

        
      }
      print(paste("community", i, ":"))
      print(it)
      print((sum(bench[it,it])-length(it))/(length(it)*length(it)-length(it)))
      ali<-align(integrate[it,it], bench[it,it]
                 )
      
      score.align[[i]] <- (ali$score)
      score.count[[i]] <- ((sum(bench[it,it])-length(it))/(length(it)*length(it)-length(it)))
      size.intersect[[i]] <- length(it)
      
      if(benchname=="target"){
        
        community_info <- getdts(it,ctrpDrTargs,savename)
        target_genes <- capture.output(cat(community_info$TARGET_NAME, sep=";")) 
        comm.targets[i] <- unique(strsplit(target_genes,";")[1])
        
        ## Aggragerate rows for uniport and atc
        if(dim(integrate)[1]==238){
          
          print("")
          agg_info <- aggregate(community_info,by=list(community_info$MOLECULE_NAME),FUN=function(x){return(paste(x,collapse = ";"))})
          
          community_info <- data.frame(agg_info[,c(1,3)])
          
          colnames(community_info)[1]<-"MOLECULE_NAME"
        }

      }else{
        
        result <- getats(it,chembl_ATC)
        community_info <-result$agg
        target_genes <- capture.output(cat(community_info$TARGET_NAME, sep=";")) 
        comm.targets[i] <- unique(strsplit(target_genes,";")[1])

      }
      
      if((novel == TRUE) & 
         length(setdiff(names.taxonomy,it)>=0) #&
         #score.count[[i]]>0.5
         
         ){
        
        if(ncol(community_info)==3){
          bind.names<-cbind(names.novel[[i]],names.novel[[i]],names.novel[[i]])
          bind.names[,c(2,3)]<-"-"
          
        }else{
          bind.names<-cbind(names.novel[[i]],names.novel[[i]])
          bind.names[,c(2)]<-"-"
          
        }                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                  
        
        colnames(bind.names)<-colnames(community_info)
        community_info<-rbind(community_info,bind.names)
        
        # if(dim(integrate)==239){
        #   df.druginfo<-read.csv("Data/druginfo_ctrp.csv")
        #   df.druginfo<-df.druginfo[,c(2,3,4)]
        # }else{
        #   df.druginfo<-read.csv("Data/druginfo_nci60.csv")
        #   df.druginfo<-df.druginfo[,c(2,3,4)]
        #   
        #   
        # }
        df.merge<-merge(community_info,df.druginfo)
      }
      
      # Comment by jy to change the exeamplars to the first element
      if ((!file.exists(saveFileName))){
        write.xlsx(data.frame(df.merge), file=saveFileName, 
                   #sheetName=paste(names(apcomb@exemplars[i]),"-community",i,sep=''),
                   #sheetName=paste(names(apcomb[apcomb==i])[1],"-community",i,sep=''),
                   sheetName=paste(exeamplars[[i]],"-community",i,sep=''),
                   row.names=FALSE)
        
      }else if(overwrite==TRUE){
        write.xlsx(data.frame(df.merge), file=saveFileName, 
                   #sheetName=paste(names(apcomb@exemplars[i]),"-community",i,sep=''), 
                   sheetName=paste(exeamplars[[i]],"-community",i,sep=''),
                   append=TRUE, row.names=FALSE)
      }
    }
    
  }
  
  print(mean(na.omit(score.align)))
  
  return(list(align=score.align,count=score.count,novel=names.novel,intsize= size.intersect,info=comm.targets))
}

getdts <- function(cdrugs,ctrpDrTargs,source = "ctrp"){
  badchars <- "[\xb5]|[\n]|[,]|[;]|[:]|[-]|[+]|[*]|[%]|[$]|[#]|[{]|[}]|[[]|[]]|[|]|[\\^]|[/]|[\\]|[.]|[_]|[ ]"
  
  if(length(grep("ctrp",source))!=0){
  
  # ctrpDrTargs <- read.csv("./Data/CTRPv2_drugtarget.csv", stringsAsFactors = FALSE) # 481 drugs x 3 descriptions
  # ctrpDrTargs$compound_name <- toupper(ctrpDrTargs$compound_name)
  # ctrpDrTargs$compound_name <- gsub(badchars,"",ctrpDrTargs$compound_name)
  ## 
  DTargs <- ctrpDrTargs[ctrpDrTargs$compound_name %in% cdrugs,,drop=F] 
  DTargs <- DTargs[,c(1,2,3)] 
  colnames(DTargs) <- c("MOLECULE_NAME","TARGET_NAME","ACTIVITY")
  }else{
    
    ## read CHEMBL drug file downloaded from Chembl website
    chemblDrTargs <- read.csv("./Data/uniprot links.csv", stringsAsFactor=F, na.strings=c("", "NA")) #2043 entries
    drgTargets <- chemblDrTargs[,c("Name", "UniProt.Name")]
    colnames(drgTargets) <-  c("MOLECULE_NAME","TARGET_NAME")
    ## remove badchar from Chembl and common drug file + capitalize 
    drgTargets[,1] <- gsub(badchars, "",  drgTargets[,1])
    drgTargets[,1] <- toupper(drgTargets[,1])
    DTargs <- drgTargets[drgTargets[,1] %in% cdrugs,]
  }
  
  return(DTargs)
}



getats <- function(cdrugs,chembl_ATC){
  # badchars <- "[\xb5]|[\n]|[,]|[;]|[:]|[-]|[+]|[*]|[%]|[$]|[#]|[{]|[}]|[[]|[]]|[|]|[\\^]|[/]|[\\]|[.]|[_]|[ ]"
  # 
  # chembl_ATC <- read.delim("Data/chembl_drugtargets-16_5-10-02.txt", stringsAsFactor=F, na.strings=c("", "NA")) #10506 drugs!
  # ## remove unused strings from synonyms or drug names such as ()
  # sep = " \\("
  # #chembl_ATC$MOLECULE_NAME <- toupper(unlist(lapply(chembl_ATC$MOLECULE_NAME, function(x) strsplit(x, sep)[[1]][1])))
  # 
  # 
  # chembl_ATC$MOLECULE_NAME <- toupper(unlist(lapply(chembl_ATC$MOLECULE_NAME, function(x) strsplit(x, sep)[[1]][1])))
  
  atclist <- chembl_ATC[chembl_ATC$MOLECULE_NAME %in% cdrugs,c("MOLECULE_NAME",
                                                               #"MECHANISM_OF_ACTION",
                                                               "TARGET_NAME","ATC_CODE")]
  
  aggatc <- aggregate(atclist,by=list(atclist$MOLECULE_NAME),FUN=function(x){return(paste(x,collapse = ";"))})
  
  aggatc <- data.frame(aggatc[,c(1,3,4)])
  
  colnames(aggatc)[1]<-"MOLECULE_NAME"

  return(list(agg=aggatc,raw=atclist))
}

# Edit by jy to change apcomb to labels
community_query <- function(drug,apcomb){
  
  #for(i in 1:length(apcomb@clusters)){
    
  for(i in 1:max(apcomb)){
      
    members<-names(apcomb[apcomb==i]) 
    #members <- names(apcomb@clusters[[i]])
    
    
    if((drug %in% members)==TRUE){
      print(i)
      
      return(members)
    }
  }
}

novel.filter <- function(r1,r2,thres=3){
  
  return(union(unlist(r2$novel[(r2$intsize)>=thres]),unlist(r1$novel[(r1$intsize)>=thres])))
}

r1<-commualign(ctrp,dtgbench.ctrp,apcomb.ctrp,threshold = 2,benchname = 'target',savename = name.save.ctrp,plot = TRUE)
r2<-commualign(ctrp,ATCbench.ctrp,apcomb.ctrp,threshold = 2,benchname = 'atc',savename = name.save.ctrp)
r3<-commualign(nc60,dtgbench.nci60,apcomb.nc60,threshold = 2,benchname = 'target',savename = name.save.nc60,plot = TRUE)
r4<-commualign(nc60,ATCbench.nci60,apcomb.nc60,threshold = 2,benchname = 'atc',savename = name.save.nc60)

#r5<-commualign(ctrp,dtgbench.ctrp,louvain.ctrp,threshold = 2,benchname = 'target',savename = name.save.ctrp)

#r5<-commualign(ctrp.dnf,dtgbench.ctrp,apcomb.ctrp.dnf,threshold = 2,benchname = 'target',savename = name.save.ctrp.dnf)
#r6<-commualign(ctrp.dnf,ATCbench.ctrp,apcomb.ctrp.dnf,threshold = 2,benchname = 'atc',savename = name.save.ctrp.dnf)
#r7<-commualign(nc60.dnf,dtgbench.nci60,apcomb.nc60.dnf,threshold = 2,benchname = 'target',savename = name.save.nc60.dnf)
#r8<-commualign(nc60.dnf,ATCbench.nci60,apcomb.nc60.dnf,threshold = 2,benchname = 'atc',savename = name.save.nc60.dnf)

names.benchctrp<-union(rownames(dtgbench.ctrp),rownames(ATCbench.ctrp))
names.benchnc60<-union(rownames(dtgbench.nci60),rownames(ATCbench.nci60))
names.benchall <- union(names.benchctrp,names.benchnc60)


union.novel.ctrp <- novel.filter(r1,r2,3)
union.novel.nc60 <- novel.filter(r3,r4,3)

novel.all <-union( 
  union.novel.ctrp,union.novel.nc60
)


novel.ctrp <- setdiff(union.novel.ctrp,names.benchall)
novel.nc60 <- setdiff(union.novel.nc60,names.benchall)

cv.ctrp <- intersect(union.novel.ctrp,names.benchnc60)
cv.nc60 <- intersect(union.novel.nc60,names.benchctrp)

members.ctrp <- list()
members.nc60 <- list()
members.cv.ctrp <- list()
members.cv.nc60 <- list()

for(i in 1:length(novel.ctrp)){
  members.ctrp[[i]] <- community_query(novel.ctrp[[i]],apcomb = apcomb.ctrp)
}
for(i in 1:length(novel.nc60)){
  members.nc60[[i]] <- community_query(novel.nc60[[i]],apcomb = apcomb.nc60)
}


for(i in 1:length(cv.ctrp)){
  print(cv.ctrp[[i]])
  members.cv.ctrp[[i]] <- community_query(cv.ctrp[[i]],apcomb = apcomb.ctrp)
  members.cv.ctrp[[i]] <- community_query(cv.ctrp[[i]],apcomb = apcomb.nc60)
  print("-----")
}
for(i in 1:length(cv.nc60)){
  print(cv.nc60[[i]])
  members.cv.nc60[[i]] <- community_query(cv.nc60[[i]],apcomb = apcomb.nc60)
  members.cv.nc60[[i]] <- community_query(cv.nc60[[i]],apcomb = apcomb.ctrp)
  print("-----")
  
}