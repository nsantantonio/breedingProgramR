simName <- "testSim"
# Should phenotypes from the founding population be used to update GS model?
# Seems like this would typically be true, but for now I am going to omit them
# such that they are onyl used for the first year of GS.
seed <- 12345 # can be null to randomly select
includeFounderPheno <- FALSE 


# define founder pop
nFounder <- 100
nNuclear <- 50
nChrom <- 10
nLoci <- 100

# define snp/qtl loci
nM <- floor(c(10*(9:1), 4) / 4)
nQTL <- c(2*(9:1), 1)

# traits, heritability etc
Vg <- 1
Vgxe <- 1

founderh2 <- 0.3

h2 <- c(0.15, 0.3, 0.3, 0.3, 0.3)
# h2 <- 0.3
# h2hr <- h2 / 2

# generations, selection intensity
nYr <- 20
nTrial <- 5

# nDHfam <- 50 
# DHfamSize <- 20
# selOutOfHR <- 500
# selOutOf1plot <- 100
# selOutOf3plot <- 30
# selVariety <- 3

nDHfam <- 10 
DHfamSize <- 20
# selOutOfHR <- 50
# selOutOf1plot <- 10
# selOutOf3plot <- 5
# selVariety <- 2


# nLocPrelim <- 1
# nLocAdv <- 3
# nLocMET <- 10

trialReps <- c(1, 1, 3, 3, 3)
trialLocs <- c(1, 1, 3, 10, 10)
selectTrials <- c(50, 10, 5, 5, 2)

# nRepsPerLocPrelim <- 1
# nRepsPerLocAdv <- 3
# nRepsPerLocMET <- 3

GScylcePerYr <- 3

simpleFounder <- TRUE
RGSCselect <- TRUE

# returnVDPtoRGSC <- c("headrow", "prelim")
returnVDPtoRGSC <- c(2, 3)

# only include the last 5 generations for prediction
lgen <- 5
# lgen <- nYr

# other prams.
burnInYears <- 0
nGen <- 20
# generation <- 1:nGen
RGSCintensity <- 0.2
reps <- 10

# h2toVe <- function(h2, Vg = 1) Vg * (1-h2) / h2

# RGSCuse <- if(RGSCselect) "ebv" else  "rand"

# # create simple founder pop
# if(simpleFounder) {
# 	founderPop <- quickHaplo(nFounder, nChrom, nLoci, genLen = 1, ploidy = 2L, inbred = FALSE)
# } else {
# 	load("testAlphaSimR1000SegSite.RData")
# }


# # Setting Simulation Parameters
# SP = SimParam$new(founderPop)
# # SP$nThreads <- 40

# # C0 <- founderPop

# # add trait 
# SP$addTraitA(nQtlPerChr=nQTL, mean = 0, var = Vg)
# # SP$addTraitAG(nQtlPerChr=nQTL, mean = 0, var = Vg, varGxE = Vgxe, varE = h2toVe(h2), corGxE = )

# # add error variance as function of heritability
# SP$setVarE(h2 = h2)

# # add genotype platform
# SP$addSnpChip(nM)

# loci <- pullLoci(SP)
# simParam <- SP