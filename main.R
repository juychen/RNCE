
## The code assumes that the working directory is "RNCE" folder
rm(list=ls())

# #Generate comparing networks
source("RCode/generatecompare.R")
#
# #Test paramerers
source("RCode/paratest-ctrpv-lincs.R")
source("RCode/paratest-nci60-lincs.R")

# #Generate networks and generated ROC PRC plots
source("RCode/main-ctrpv-lincs.R")
source("RCode/main-nci60-lincs.R")
#
# #Test ARI of community detection and test the number of communities
source("RCode/communityARI.R")
#
# #Generate community plots and community targets
source("RCode/main-network-generation-ctrpv2.R")
source("RCode/main-network-generation-nci60.R")
source("RCode/communityTargets.R")
source("RCode/CommunityEnrichmentNovel.R")
#
# #Analysis the correlation
source("RCode/CorrelationAnalysis.R")
#
# #Unit test
source("RCode/main-ctrpv-layertest.R")
source("RCode/main-nci60-layertest.R")
