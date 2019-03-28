library(AlphaSimR)

# Creating Founder Haplotypes
founderPop = runMacs(nInd=1000, nChr=10, segSites=1000, inbred = TRUE, species = "MAIZE")
save(founderPop, file = "testAlphaSimR1000SegSite.RData")


founderPop = runMacs(nInd=1000, nChr=10, segSites=c(10*(9:1), 5), inbred = TRUE, species = "MAIZE")
save(founderPop, file = "testAlphaSimR455SegSite.RData")

# founderPop = runMacs(nInd=1000, nChr=10, segSites=c(4*(9:1), 2), inbred = TRUE, species = "MAIZE")
# summary(founderPop)



load("testAlphaSimR.RData")
