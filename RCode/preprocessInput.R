###############################################################################################################
## Function reads in the input datasets, preprocess and finds the set of drugs intersecting between the two sets
##
## input: 
##     dname: name ("character") of the first dataset, e.g., "nci60" , "ctrpv2"   
##     d2: currently, must be set to "lincs"  
## output: 
##     a list ("character") or "dataframe" (if dname==nci60) containing the intersection of the two datasets 			
##
## 
###############################################################################################################



preprocessInput <- function(dname , d2="lincs", venn_plot = FALSE) {

   ## Read LINCS metadata file
   if (d2 == "lincs") {  
      lincs <- read.csv("Data/LINCS.csv",stringsAsFactors = FALSE) # 20326 drugs x 28 descriptions
      ## capitalize + remove badchars from lincs
      lincs$pert_iname <- toupper(lincs$pert_iname)
      lincs$pert_iname <- gsub(badchars,"",lincs$pert_iname)
   } else { stop("ERROR!")}

   if (dname == "ctrpv2") {
      ## Read ctrpv2 metadata file
      ctrpv2 <- read.csv("Data/CTRPv2_drugtarget.csv",stringsAsFactors = FALSE) # 481 drugs x 28 descriptions
      ctrpv2$compound_name <- toupper(ctrpv2$compound_name)
      ctrpv2$compound_name <- gsub(badchars,"",ctrpv2$compound_name)
      #intersect LINCS chemicals from LINCS with ctrpv2 compounds
      intrsctLincsCtrpv2 <- intersect(lincs$pert_iname, ctrpv2$compound_name) # 239 drugs/481 total
      lincsInters <- lincs[lincs$pert_iname %in% intrsctLincsCtrpv2,,drop=F] ## NOTE: 318 drugs found in lincs with the same names and different ids... 
      lincsInters <- lincsInters[!duplicated(lincsInters$pert_iname),,drop=F] ## 239 drugs
      intrsct <- lincsInters
      
      # Plot venn diagram to visualize the data
      flog.threshold(ERROR)
      venn.plot <- venn.diagram(
         x = list(
            "L1000" = lincs$pert_iname,
            "CTRP" =  ctrpv2$compound_name
         ),
         filename = NULL, 
         scaled = FALSE,
         fill = c("blue","red"),
         cat.fontfamily = rep("sans", 2),
         cat.post = 0,
         cat.dist = 0.03,
         #cat.default.pos ="outer",  
         # main = "Venn Diagram for Sensitivity Layer and Perturbation Layer",
         # main.fontfamily ="sans",
         # sub = "L1000 data and CTRPv2 data",
         # sub.fontfamily ="sans",
         
         main.cex = 2,
         sub.cex = 1,
         fontfamily =  "sans"
         
      )
      pdf(file="Output/Venn_ctrp_l1000_.pdf")
      grid.draw(venn.plot)
      dev.off()
      emf(file="Output/Venn_ctrp_l1000_.emf")
      grid.draw(venn.plot)
      dev.off()
      
   } else if (dname== "gdsc") {
           x <- read.csv("GDSC1000_druginfo.csv")
           x$DRUG.NAME <- toupper(x$DRUG.NAME)
           x$DRUG.NAME <- gsub(badchars, "", x$DRUG.NAME)
           length(unique(x$DRUG.NAME))
           #[1] 251
     
           length(intersect(x$DRUG.NAME,lincs$pert_iname))
           ii <- intersect(x$DRUG.NAME,lincs$pert_iname)
           #[1] 148
     
           lincsInters <- lincs[lincs$pert_iname %in% ii,,drop=F] ## NOTE: 318 drugs found in lincs with the same names and different ids... 
     
           # ## generate the drug pert object here for now ...
           load("Data/L1000_compound_signatures.RData")
           pert <- L1000_compounds.perturbation
           ## subset according to the pert_ids above
           xx <- colnames(pert[,,1])
           ind <- which(xx %in% lincsInters$pert_id)
           all(xx[ind] %in% lincsInters$pert_id)
           
           l1000.drug.signaturesGDSC <- pert[,ind,]
           save(l1000.drug.signaturesGDSC, file="Data/DNF_pert_all-gdsc.RData")
           intrsct <- lincsInters
           
     
   }
 else if (dname == "nci60") { 
          nci <- read.csv("Data/DTP_NCI60_ZSCORE.csv",stringsAsFactors = FALSE, na.strings="na") # 20861 drugs x 66 descriptions
          nci <- nci[!nci$Drug.name=="-",,drop=F]
          nci$Drug.name <- toupper(nci$Drug.name)
          nci$Drug.name  <- gsub(badchars,"", nci$Drug.name)
          
          ## intersect LINCS cpds with NCI60 compounds by names or smiles to maximize number of drugs
          intrscLincsNciDrugNamLst <- intersect(lincs$pert_iname, nci$Drug.name) # 192, all unique by default
          intrscLincsNciSmilLst <- intersect(lincs$canonical_smiles, nci$SMILES..d.) # 106
          
          ## keep only the intersecting drugs by names and match the names
          lincsIntrsNciDrugNam <-lincs[match(intrscLincsNciDrugNamLst,lincs$pert_iname),] ## 192 X 28
          nciIntrsLincsDrugNam <- nci[match(intrscLincsNciDrugNamLst,nci$Drug.name),] ## 192 X 66
          
          ## keep only the intersecting drugs by smiles and match smiles
          lincsIntrsNciSmil <- lincs[match(intrscLincsNciSmilLst, lincs$canonical_smiles),] ## 106 X 28 
          nciIntrsLincsSmil <- nci[match(intrscLincsNciSmilLst,nci$SMILES..d.),] ## 106 X 66
          
          ### match between nci and lincs and assign same name
          nciIntrsLincsSmil <- nciIntrsLincsSmil[match(lincsIntrsNciSmil$canonical_smiles, nciIntrsLincsSmil$SMILES..d.), ] # 106 X 66
          nciIntrsLincsSmil$Drug.name <- lincsIntrsNciSmil$pert_iname

                    
          ######## rbind lincs files/rbind nci60 files which keeps drugs in common = 238
          lincsboth <- rbind(lincsIntrsNciDrugNam, lincsIntrsNciSmil)
          lincsboth <- lincsboth[!duplicated(lincsboth$pert_iname),,drop=F]
          
          nciboth <- rbind(nciIntrsLincsDrugNam, nciIntrsLincsSmil)
          nciboth <- nciboth[!duplicated(nciboth$Drug.name),,drop=F]
          

          intrsct <- list(lincsboth=lincsboth, nciboth=nciboth)
          
          # Plot venn diagram to visualize the data
          ## gen veen
          flog.threshold(ERROR)
          venn.plot<-draw.pairwise.venn(area1=length(lincs$canonical_smiles),
                                        area2=length(nci$Drug.name),
                                        cross.area=dim(lincsboth)[1],
                                        scaled=FALSE,
                                        fill = c("blue", "red"),
                                        category = c("L1000","NCI60"),
                                        cat.fontfamily = rep("sans", 2),
                                        cat.post = 0,
                                        cat.dist = 0.03,
                                        #cat.default.pos ="text",  
                                        main = "Complex Venn Diagram",
                                        sub = "Featuring: rotation and external lines",
                                        main.cex = 2,
                                        sub.cex = 1,
                                        fontfamily =  "sans"
          )
          pdf(file="Output/venn_nci60_L1000.pdf")
          grid.draw(venn.plot)
          dev.off()
          
          emf(file="Output/venn_nci60_L1000.emf")
          grid.draw(venn.plot)
          dev.off()
          
 }
  
  
  return(intrsct) 

}