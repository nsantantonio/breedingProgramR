# get parent directory
parDir <- getwd()
# get functions
source(paste0(parDir, "/alphaTools.R"))

# set default arguments
defArgs <- list(
# simulation parameters
	seed = 12345,
	nThreads = 11, 
	projName = NULL, 
	simName = "simNoTitle",
	simFunc = "fitFuncSingleVariate.R", 
	nYr = 20,
	reps = 10,
	lgen = 5,
	useTrue = FALSE,
	traditional = FALSE, # this selects out of VDP as parents, no RGSC
	# nSimCrosses = 10, 
# founder parameters
	founderRData = "founderPop/testAlphaSimR1000SegSite.RData",
	founderh2 = 0.3,
	simpleFounder = FALSE,
	founderBurnIn = 1,
# selection parameters
	selectRGSC = 0.2,
	nProgenyPerCrossIn = 1,
	nProgenyPerCrossOut = 1,
	selectIn = "ebv", # ebv, rand
	selectOut = "ebv", # ebv, var, exp?
	selectVDP = "pheno", # ebv, pheno
	returnVDPcrit = "pheno", # ebv?
	selFuncOut = NULL, # truncSel, expDist, simDHdist
	selFuncIn = simDHdistPairs, # truncCross, expDistPairs, simDHdistPairs, maxVar
	withinFamInt = 1, #  
	setXint = NULL, # note that x is the cdf of a normal 
	skip = NULL,
	# weight = 0.5,
# family parameters
	nFounder = 10,
	nNuclear = 100,
	nFam = 10,
	famSize = 50,
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
# other args - need to be able to add additional things here that 
)

stdArgNames <- names(defArgs)
defArgs <- getComArgs(defArgs)
altArgs <- names(defArgs)[!names(defArgs) %in% stdArgNames] 
attach(defArgs)

# load libraries
library(AlphaSimR)
library(doMC)
library(txtplot) 

# source simulation function
source(simFunc)

# print selection sizes
nFam * famSize * c(1, cumprod(selectTrials))

if(!is.null(projName)) projName <- paste0(projName, "/")
# make directory to store all results
simDir <- paste0(parDir, "/results/", projName, simName) 
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
# Reduce(intersect, loci)

testRun <- FALSE
if(testRun){
	defArgs[["maxCrossPerParent"]] <- 1
	defArgs[["weight"]] <- 0.2
	altArgs <- c(altArgs, "maxCrossPerParent", "weight")
	run1 <- do.call(sim, c(list(founderPop = founderPop, simParam = SP, paramL = defArgs, returnFunc = getPopStats), defArgs[altArgs]))
	run1$gv
	run1$vy
}

if(system("hostname", intern = TRUE) == "Bender") {
	setMKLthreads(1)
	registerDoMC(nThreads)
} else {
	registerDoMC(nThreads)
}

# simrun <- foreach(r = 1:reps) %dopar% sim(founderPop, simParam = SP, paramL = defArgs, returnFunc = getPopStats)
simrun <- foreach(r = 1:reps) %dopar% do.call(sim, c(list(founderPop = founderPop, simParam = SP, paramL = defArgs, returnFunc = getPopStats), defArgs[altArgs]))
save(simrun, SP, file = paste0(simDir, "/", simName, ".RData"))
