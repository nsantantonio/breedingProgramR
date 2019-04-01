simName <- "testSim"
# Should phenotypes from the founding population be used to update GS model?
# Seems like this would typically be true, but for now I am going to omit them
# such that they are onyl used for the first year of GS.
seed <- 12345 # can be null to randomly select
# includeFounderPheno <- FALSE 

RGSCselect <- TRUE
selF2 <- FALSE
nF2 <- 100

# define founder pop
simpleFounder <- FALSE


nFounder <- 10
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
h2 <- c(0.1, 0.3, 0.3, 0.3, 0.3)

# generations, selection intensity
nYr <- 20
nFam <- 10 
famSize <- 100

# selectTrials <- c(50, 10, 5, 5, 2)
selectTrials <- c(0.50, 0.2, 0.2, 0.5, 0.3)
nDH <- nDHfam * DHfamSize
if(all(selectTrials < 1)) cat("n selected =", nDH * cumprod(selectTrials), "\n")


trialReps <- c(1, 1, 3, 3, 3)
trialLocs <- c(1, 1, 3, 10, 10)


# with 2 cycles per year, this is more realitic, and allows for direct 
#comparison with F2 sleection (i.e. 1 gen of selfing after cross, then 
# sample few F2 families with high expected quantiles)
GScylcePerYr <- 2 

# returnVDPtoRGSC <- c("headrow", "prelim")
returnVDPtoRGSC <- c(2, 3)

# only include the last 5 generations for prediction
lgen <- 3
# lgen <- nYr

# other prams.
# burnInYears <- 0
nGen <- 20

RGSCintensity <- 0.2

# total number of replication
reps <- 10
