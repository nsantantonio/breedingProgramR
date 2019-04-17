# get parent directory
parDir <- getwd()
# get functions
source(paste0(parDir, "/alphaTools.R"))

# set default arguments
defArgs <- list(
# simulation parameters
	seed = 12345,
	nThreads = 11, 
	simName = "simNoTitle",
	simFunc = "fitFuncSingleVariate.R", 
	nYr = 20,
	reps = 10,
	lgen = 5,
	useTrue = FALSE,
	traditional = FALSE, # this selects out of VDP as parents, no RGSC
	nSimCrosses = 10, 
# founder parameters
	founderRData = "founderPop/testAlphaSimR1000SegSite.RData",
	founderh2 = 0.3,
	simpleFounder = FALSE,
	founderBurnIn = 3,
# selection parameters
	selectRGSC = 0.2,
	nProgenyPerCrossIn = 1,
	nProgenyPerCrossOut = 1,
	selectIn = "ebv", # ebv, rand
	selectOut = "ebv", # ebv, var, exp?
	selectVDP = "pheno", # ebv, pheno
	returnVDPcrit = "pheno", # ebv?
	selFuncOut = truncSel, # truncSel, expDist, simDHdist
	selFuncIn = expDistPairs, # truncCross, expDistPairs, simDHdistPairs
	withinFamInt = 1, # none, 
	setXint = NULL, # note that x is the cdf of a normal 
	skip = NULL,
# family parameters
	nFounder = 100,
	nNuclear = 100,
	nFam = 10,
	famSize = 20,
	ssd = FALSE,
	selF2 = FALSE,
	nF2 = 1,
# genetic parameters
	# kinship = "SNP",
	Vg = 1,
	updateVg = FALSE,
	# Vgxe = 1,
	h2 = c(0.1, 0.3, 0.3, 0.3, 0.3),
	selectTrials = c(0.5, 0.5, 0.4, 0.3, 1/3),
	trialReps = c(1, 1, 2, 3, 3),
	trialLocs = c(1, 1, 2, 5, 5),
	cyclePerYr = 2,
	returnVDPtoRGSC = c(0, 0, 0, 0, 0, 0), # default to rep(0, nTrial)?
# chromosome parameters
	nChrom = 10,
	nLoci = 100,
	nM = 100, # floor(c(10*(9:1), 4) / 4),
	nQTL = 100 # c(2*(9:1), 1),
)


# expDist and truncCross seemed to do well earlier. Need to rerun and check. 



defArgs <- getComArgs(defArgs)
attach(defArgs)

# load libraries
library(AlphaSimR)
library(doMC)
library(txtplot) 
# source(paramFile)

print(getwd())

source(simFunc)


nFam * famSize * cumprod(selectTrials)


# make directory to store all results
simDir <- paste0(parDir, "/results/", simName) 
system(paste0("simdir=", simDir, "
if [ ! -d $simdir ]; then
	mkdir -p $simdir
else 
	echo 'Directory already exists! Any results therein will be overwritten without further notice...'
fi"))


# set seed
if(!is.null(seed)) set.seed(seed)


# create simple founder pop, either simple or from loaded MaCS sim
if(simpleFounder) {
	founderPop <- quickHaplo(nFounder, nChrom, nLoci, genLen = 1, ploidy = 2L, inbred = FALSE)
} else {
	load(founderRData)
	founderPop <- founderPop[sample(1:nInd(founderPop), nFounder)]
}

# Setting Simulation Parameters
SP = SimParam$new(founderPop)
# SP$nThreads <- 40

# add trait 
SP$addTraitA(nQtlPerChr=nQTL, mean = 0, var = Vg)
# SP$addTraitAG(nQtlPerChr=nQTL, mean = 0, var = Vg, varGxE = Vgxe, varE = h2toVe(h2), corGxE = )

# add error variance as function of heritability
SP$setVarE(h2 = founderh2)

# add genotype platform
SP$addSnpChip(nM)

loci <- pullLoci(SP)
Reduce(intersect, loci)

# why does expDistPairs fail even when w = 1??!!!!
run1 <- sim(founderPop, simParam = SP, paramL = defArgs, returnFunc = getPopStats, w = 0)
run1$gv
variety <- nFam * famSize * cumprod(selectTrials)
nv <- variety[length(variety)]
tapply(run1$vy, rep(1:nYr, each = nv), mean)

run1$Vg


if(system("hostname", intern = TRUE) == "Bender") {
	setMKLthreads(1)
	registerDoMC(nThreads)
} else {
	registerDoMC(nThreads)
}

# simrun <- foreach(r = 1:reps) %do% sim(founderPop, simParam = SP, paramL = defArgs, returnFunc = getPopStats)

simrun <- foreach(r = 1:reps) %dopar% sim(founderPop, simParam = SP, paramL = defArgs, returnFunc = getPopStats)
save(simrun, SP, file = paste0("results/", simName, "/", simName, ".RData"))
