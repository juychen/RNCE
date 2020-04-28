
rm(list=ls())
library(devEMF)


# Deena M.A. Gendoo
# Determining the Correlation between different types of data (single-layer drug taxonomy)
# Do this by computing the spearman correlation between all pairs of similarity matrices
# September 22, 2016
################################################################################################
################################################################################################
badchars <- "[\xb5]|[\n]|[,]|[;]|[:]|[-]|[+]|[*]|[%]|[$]|[#]|[{]|[}]|[[]|[]]|[|]|[\\^]|[/]|[\\]|[.]|[_]|[ ]"

## creating the output directory if not exists
if ( ! file.exists("Output")) {
  dir.create("Output")
}
#datestr <- date()
datestr <- Sys.Date()
mainDir <- getwd()
datestr <- gsub(badchars, "",  datestr)
subDir<-paste("Output/",datestr,sep = "")
options(someUniqueTag.mainDir = mainDir)
options(someUniqueTag.subDir = "subDir")

if(!file_test("-d", file.path(mainDir, subDir))){
  if(file_test("-f", file.path(mainDir, subDir))) {
    stop("Path can't be created because a file with that name already exists.")
  } else {
    dir.create(file.path(mainDir, subDir))
  }
}


bp<-matrix(ncol = 1,nrow = 6)
bp[1]<-bp[2]<-bp[3]<-bp[4]<-bp[5]<-bp[6]<-0.15

#CTRPv2-based DNF
################################
#load Structure Affinity Matrix
load("Data/strcSimilarityMat-ctrpv2.RData")
#load Sensitivity Affinity Matrix
load("Data/sensSimilarityMat-ctrpv2.RData")
#load Perturbation Affinity Matrix
load("Data/pertSimilarityMat-ctrpv2.RData")
#load DNF itself!
load("Data/integrationSimilarityMat-ctrpv2.RData")
temp<-integrtStrctSensPert
#load RERANK itself!
load("Data/rerkSimilarityMat-CTRP.Rdata")
load("Data/ctrpv2-DNFIntegrated.RData")
DNF<-integrtStrctSensPert
integrtStrctSensPert<-temp


#compute correlations between pairs of data types
#Single Layers against themselves
# Structure-Sensitivity
Rerank.Vs.Sens<-cor.test(re_ranked,sensAffMat,method = "spearman",alternative = "two.sided",exact=FALSE)
# Structure - Perturbation
Struct.Vs.Rerank<-cor.test(strcAffMat,re_ranked,method = "spearman",alternative = "two.sided",exact=FALSE)
# Perturbation - Sensitivity
Pert.Vs.Rerank<-cor.test(pertAffMat,re_ranked,method = "spearman",alternative = "two.sided",exact=FALSE)
# Single layers against DNF
# DNF- Pert
integrate.Vs.Rerank<-cor.test(integrtStrctSensPert,re_ranked,method = "spearman",alternative = "two.sided",exact=FALSE)

Pert.Vs.intg<-cor.test(pertAffMat,integrtStrctSensPert,method = "spearman",alternative = "two.sided",exact=FALSE)
Strc.Vs.intg<-cor.test(strcAffMat,integrtStrctSensPert,method = "spearman",alternative = "two.sided",exact=FALSE)
Sens.Vs.intg<-cor.test(sensAffMat,integrtStrctSensPert,method = "spearman",alternative = "two.sided",exact=FALSE)

Pert.Vs.DNF<-cor.test(pertAffMat,DNF,method = "spearman",alternative = "two.sided",exact=FALSE)
Strc.Vs.DNF<-cor.test(strcAffMat,DNF,method = "spearman",alternative = "two.sided",exact=FALSE)
Sens.Vs.DNF<-cor.test(sensAffMat,DNF,method = "spearman",alternative = "two.sided",exact=FALSE)

Pert.Vs.Strc<-cor.test(pertAffMat,strcAffMat,method = "spearman",alternative = "two.sided",exact=FALSE)
Strc.Vs.Sens<-cor.test(strcAffMat,sensAffMat,method = "spearman",alternative = "two.sided",exact=FALSE)
Sens.Vs.Pert<-cor.test(sensAffMat,pertAffMat,method = "spearman",alternative = "two.sided",exact=FALSE)

Rerank.vs.DNF <- cor.test(re_ranked,DNF,method = "spearman",alternative = "two.sided",exact=FALSE)
integrate.vs.DNF <- cor.test(re_ranked,integrtStrctSensPert,method = "spearman",alternative = "two.sided",exact=FALSE)

# 
# CTRPv2_Cors<-c(Struct.Vs.Rerank$estimate,Rerank.Vs.Sens$estimate,Pert.Vs.Rerank$estimate,NA,
#                integrate.Vs.Rerank$estimate,integrate.Vs.Rerank$estimate,integrate.Vs.Rerank$estimate)


CTRPv2_Cors<-c(
  Strc.Vs.DNF$estimate,Sens.Vs.DNF$estimate,Pert.Vs.DNF$estimate,
  
  Struct.Vs.Rerank$estimate,Rerank.Vs.Sens$estimate,Pert.Vs.Rerank$estimate,
              Strc.Vs.intg$estimate,Sens.Vs.intg$estimate,Pert.Vs.intg$estimate,
              integrate.Vs.Rerank$estimate,Rerank.vs.DNF$estimate,integrate.vs.DNF$estimate)

# names(CTRPv2_Cors)<-c("Context & Struct","Context & Sens","Context & Pert","", 
#                       "Inte & Pert","DNF & Sens", "DNF & Struct")

CTRPv2_inputCors<-c(
  Pert.Vs.Strc$estimate,Strc.Vs.Sens$estimate,Sens.Vs.Pert$estimate)

names(CTRPv2_inputCors)<-c(
  "Pert & Struct","Struct & Sens","Sens & Pert"
  )



names(CTRPv2_Cors)<-c(
  "Initial & Struct","Initial & Sens","Initial & Pert",
  
  "KRNN-Context & Struct","Context & Sens","Context & Pert",
                      "RNCE & Struct","RNCE & Sens","RNCE & Pert",
                      "Context & RNCE","Context & Initial","Initial & RNCE" )

CTRPv2_PVals<-c(prettyNum(Struct.Vs.Rerank$p.value,digits=3),prettyNum(Rerank.Vs.Sens$p.value,digits=3),prettyNum(Pert.Vs.Rerank$p.value,digits=3),
                prettyNum(integrate.Vs.Rerank$p.value,digits=3),prettyNum(integrate.Vs.Rerank$p.value,digits=3),prettyNum(integrate.Vs.Rerank$p.value,digits=3))

barplot(CTRPv2_Cors,col = "#d9ef8b",las=2,ylim=c(0,1))
text(bp,labels = CTRPv2_PVals)


#NCI60-based DNF
################################
#load Structure Affinity Matrix
load("Data/strcSimilarityMat-nci60.RData")
#load Sensitivity Affinity Matrix
load("Data/sensSimilarityMat-nci60.RData")
#load Perturbation Affinity Matrix
load("Data/pertSimilarityMat-nci60.RData")
#load DNF itself!
load("Data/integrationSimilarityMat-nci60.RData")
#load RERANK itself!
load("Data/rerkSimilarityMat-nci60.Rdata")
temp<-integrtStrctSensPert
#load RERANK itself!
load("Data/nci60-DNFIntegrated.RData")
DNF<-integrtStrctSensPert
integrtStrctSensPert<-temp


#compute correlations between pairs of data types
#Single Layers against themselves
# Structure-Sensitivity
Rerank.Vs.Sens<-cor.test(re_ranked,sensAffMat,method = "spearman",alternative = "two.sided",exact=FALSE)
# Structure - Perturbation
Struct.Vs.Rerank<-cor.test(strcAffMat,re_ranked,method = "spearman",alternative = "two.sided",exact=FALSE)
# Perturbation - Sensitivity
Pert.Vs.Rerank<-cor.test(pertAffMat,re_ranked,method = "spearman",alternative = "two.sided",exact=FALSE)
# Single layers against DNF
# DNF- Pert
integrate.Vs.Rerank<-cor.test(integrtStrctSensPert,re_ranked,method = "spearman",alternative = "two.sided",exact=FALSE)

Pert.Vs.intg<-cor.test(pertAffMat,integrtStrctSensPert,method = "spearman",alternative = "two.sided",exact=FALSE)
Strc.Vs.intg<-cor.test(strcAffMat,integrtStrctSensPert,method = "spearman",alternative = "two.sided",exact=FALSE)
Sens.Vs.intg<-cor.test(sensAffMat,integrtStrctSensPert,method = "spearman",alternative = "two.sided",exact=FALSE)


Pert.Vs.DNF<-cor.test(pertAffMat,DNF,method = "spearman",alternative = "two.sided",exact=FALSE)
Strc.Vs.DNF<-cor.test(strcAffMat,DNF,method = "spearman",alternative = "two.sided",exact=FALSE)
Sens.Vs.DNF<-cor.test(sensAffMat,DNF,method = "spearman",alternative = "two.sided",exact=FALSE)

Pert.Vs.Strc<-cor.test(pertAffMat,strcAffMat,method = "spearman",alternative = "two.sided",exact=FALSE)
Strc.Vs.Sens<-cor.test(strcAffMat,sensAffMat,method = "spearman",alternative = "two.sided",exact=FALSE)
Sens.Vs.Pert<-cor.test(sensAffMat,pertAffMat,method = "spearman",alternative = "two.sided",exact=FALSE)


Rerank.vs.DNF <- cor.test(re_ranked,DNF,method = "spearman",alternative = "two.sided",exact=FALSE)
integrate.vs.DNF <- cor.test(re_ranked,integrtStrctSensPert,method = "spearman",alternative = "two.sided",exact=FALSE)

NCI60_Cors<-c(
  Strc.Vs.DNF$estimate,Sens.Vs.DNF$estimate,Pert.Vs.DNF$estimate,
  Struct.Vs.Rerank$estimate,Rerank.Vs.Sens$estimate,Pert.Vs.Rerank$estimate,
              Strc.Vs.intg$estimate,Sens.Vs.intg$estimate,Pert.Vs.intg$estimate,
              integrate.Vs.Rerank$estimate,Rerank.vs.DNF$estimate,integrate.vs.DNF$estimate)

names(NCI60_Cors)<-c(
  "Initial & Struct","Initial & Sens","Initial & Pert",
  "Context & Struct","Context & Sens","Context & Pert",
                     "RNCE & Struct","RNCE & Sens","RNCE & Pert",
                     "Context & RNCE","Context & Initial","Initial & RNCE" )

NCI60_PVals<-c(prettyNum(Struct.Vs.Rerank$p.value,digits=3),prettyNum(Rerank.Vs.Sens$p.value,digits=3),prettyNum(Pert.Vs.Rerank$p.value,digits=3),
               prettyNum(integrate.Vs.Rerank$p.value,digits=3),prettyNum(integrate.Vs.Rerank$p.value,digits=3),prettyNum(integrate.Vs.Rerank$p.value,digits=3))
NCI60_inputCors<-c(
  Pert.Vs.Strc$estimate,Strc.Vs.Sens$estimate,Sens.Vs.Pert$estimate)

names(NCI60_inputCors)<-c(
  "Pert & Struct","Struct & Sens","Sens & Pert"
)

barplot(NCI60_Cors,col = "#fee08b",las=2,ylim=c(0,1))
text(bp,labels = NCI60_PVals)


# Plot all Correlations For Both DNF taxonomies
emf(paste(subDir,"/CorrelationBetweenSingleLayersALL.emf",sep = ""),width = 3,height = 4)
par(mfrow=c(2,1),mai=c(0.4,0.5,0.1,0.1),bty="n",cex=0.8)

plot.new()
par(fig=c(0,1,0.6,1),bty="n")
rect(par("usr")[1],par("usr")[3],par("usr")[2],par("usr")[4],col = "#ebebeb")

grid (lty = 1, col = "white")
barplot(CTRPv2_Cors,col = "#ffa07a",las=2,main=NA, xaxt='n',ylim=c(0,1),new=TRUE)
#text(bp,labels = CTRPv2_PVals)


par(fig=c(0,1,0.25,0.65),bty="n")
rect(par("usr")[1],par("usr")[3],par("usr")[2],par("usr")[4],col = "#ebebeb")

grid (lty = 1, col = "white")

barplot(NCI60_Cors,col = "#48d1cc",las=2,main=NA,ylim=c(0,1),add=TRUE)

# Label the x and y axes with dark green text

#text(bp,labels = NCI60_PVals,las=2)
dev.off()

# Plot all Correlations For Both DNF taxonomies
emf(paste(subDir,"/CorrelationBetweeninput.emf",sep = ""),width = 3,height = 4)
par(mfrow=c(2,1),mai=c(0.4,0.5,0.1,0.1),bty="n",cex=0.8)

plot.new()
par(fig=c(0,1,0.6,1),bty="n")
rect(par("usr")[1],par("usr")[3],par("usr")[2],par("usr")[4],col = "#ebebeb")

grid (lty = 1, col = "white")
barplot(c(CTRPv2_inputCors,CTRPv2_Cors),col = "#ffa07a",las=2, xaxt='n',main=NA,ylim=c(0,1),new=TRUE)
#text(bp,labels = CTRPv2_PVals)


par(fig=c(0,1,0.25,0.65),bty="n")
rect(par("usr")[1],par("usr")[3],par("usr")[2],par("usr")[4],col = "#ebebeb")

grid (lty = 1, col = "white")

barplot(c(NCI60_inputCors,NCI60_Cors),col = "#48d1cc",las=2,main=NA,ylim=c(0,1),add=TRUE)

# Label the x and y axes with dark green text

#text(bp,labels = NCI60_PVals,las=2)
dev.off()


# Plot all Correlations For Both DNF taxonomies
emf(paste(subDir,"/CorrelationBetweeninput2.emf",sep = ""),width = 1.5,height = 4)
par(mfrow=c(2,1),mai=c(0.4,0.5,0.1,0.1),bty="n",cex=0.8)

plot.new()
par(fig=c(0,1,0.6,1),bty="n")
rect(par("usr")[1],par("usr")[3],par("usr")[2],par("usr")[4],col = "#ebebeb")

grid (lty = 1, col = "white")
barplot(CTRPv2_inputCors,col = "#ffa07a",las=2,main=NA, xaxt='n',ylim=c(0,1),new=TRUE)
#text(bp,labels = CTRPv2_PVals)


par(fig=c(0,1,0.25,0.65),bty="n")
rect(par("usr")[1],par("usr")[3],par("usr")[2],par("usr")[4],col = "#ebebeb")

grid (lty = 1, col = "white")

barplot(NCI60_inputCors,col = "#48d1cc",las=2,main=NA,ylim=c(0,1),add=TRUE)

# Label the x and y axes with dark green text

#text(bp,labels = NCI60_PVals,las=2)
dev.off()


