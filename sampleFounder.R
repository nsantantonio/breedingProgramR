library(AlphaSimR)
# get parent directory
parDir <- getwd()
# get functions
source(paste0(parDir, "/alphaTools.R"))

founderFile <- paste0(parDir, "/founderPop/testAlphaSimR1000SegSite.RData")
load(founderFile)

nFounder <- 100
founderSamples <- sampleFounderPop(founderPop, size = nFounder, n = 100, seed = 12345)
save(founderSamples, file = paste0(parDir, "/founderPop/founderSamples_nInd", nFounder, "_1000SegSite.RData"))

nFounder <- 10
founderSamples <- sampleFounderPop(founderPop, size = nFounder, n = 100, seed = 12345)
save(founderSamples, file = paste0(parDir, "/founderPop/founderSamples_nInd", nFounder, "_1000SegSite.RData"))
