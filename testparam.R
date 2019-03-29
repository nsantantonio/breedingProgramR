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

# h2
founderh2 <- 0.3
h2 <- c(0.15, 0.3, 0.3, 0.3, 0.3)

# generations, selection intensity
nYr <- 20
nDHfam <- 10 
DHfamSize <- 100

# selectTrials <- c(50, 10, 5, 5, 2)
selectTrials <- c(0.50, 0.2, 0.2, 0.5, 0.3)
nDH <- nDHfam * DHfamSize
if(all(selectTrials < 1)) cat("n selected =", nDH * cumprod(selectTrials), "\n")


trialReps <- c(1, 1, 3, 3, 3)
trialLocs <- c(1, 1, 3, 10, 10)


GScylcePerYr <- 3

simpleFounder <- TRUE
RGSCselect <- TRUE

# returnVDPtoRGSC <- c("headrow", "prelim")
returnVDPtoRGSC <- c(2, 3)

# only include the last 5 generations for prediction
lgen <- 3
# lgen <- nYr

# other prams.
burnInYears <- 0
nGen <- 20

RGSCintensity <- 0.2
reps <- 10
