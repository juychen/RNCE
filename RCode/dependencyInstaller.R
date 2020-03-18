
if (!requireNamespace("BiocManager", quietly = TRUE))
install.packages("BiocManager")
BiocManager::install()
# 
list.of.CRAN.packages <- c("knitr","reshape2","doParallel",
                           "Hmisc", "apcluster", "rcdk", "fingerprint", 
                           "SNFtool", "ROCR", "proxy", "PRROC", "Hmisc","doParallel")
new.CRAN.packages <- list.of.CRAN.packages[!(list.of.CRAN.packages %in% installed.packages()[,"Package"])]
if(length(new.packages)) install.packages(new.CRAN.packages,
                                          , dependencies = TRUE, clean = TRUE, Ncpus = 2, verbose = TRUE, quiet = TRUE)
# 
# 
list.of.bioC.packages <- c("R6","PharmacoGx", "annotate", "org.Hs.eg.db", "survcomp","Biostrings","kebabs","netbiov","survcomp")
# new.bioC.packages <- list.of.bioC.packages[!(list.of.bioC.packages %in% installed.packages()[,"Package"])]
# biocLite(pkgs = new.bioC.packages, Ncpus = 2, ask = FALSE)
# 
# dependencies <- c(list.of.CRAN.packages, list.of.bioC.packages)
# 
# message("=== package install success status ===")
# 
# for(package in dependencies){
#   message(paste(package, library(package, character.only = TRUE, quietly = TRAUE, logical.return = TRUE, verbose = FALSE), sep = ": "))
# }
BiocManager::install(list.of.bioC.packages)

BiocManager::install("PharmacoGx")
BiocManager::install("R6")