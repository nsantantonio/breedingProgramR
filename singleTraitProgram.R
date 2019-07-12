# get parent directory
parDir <- getwd()
# get functions
source(paste0(parDir, "/alphaTools.R"))

# set default arguments
defArgs <- list(
# simulation parameters
	seed = 12345,
	nThreads = 30, 
	projName = "testProj", 
	simName = "testSim",
	simFunc = "fitFuncSingleVariate.R", 
	maxIter = 1000L,
	reps = 10,
	nFounderPops = 10,
	lgen = 5,
	useTruth = 0, # make this an integer, 1 = use QTL, 2 = use true effects?
	traditional = FALSE, # this selects out of VDP as parents, no RGSC
# founder parameters
	founderRData = "founderPop/testAlphaSimR1000SegSite.RData",
	founderSamples = "founderPop/founderSamples_nInd100_1000SegSite.RData",
	founderh2 = 0.3,
	simpleFounder = FALSE,
	founderBurnIn = 1,
	founderReps = 1,
	founderKeep = 5,
# selection parameters
	selectRGSC = 0.2,
	pullCycle = NULL, 
	nProgenyPerCrossIn = 1,
	nProgenyPerCrossOut = 1,
	useIn = "ebv", # ebv, rand
	useOut = "ebv", # ebv, var, exp?
	useInbreed = "ebv", # ebv, var, exp?
	useVDP = "pheno", # ebv, pheno
	returnVDPcrit = "pheno", # ebv?
	selFuncOut = NULL, # truncSel, expDist, simDHdist, solqpOut
	selFuncIn = NULL, # truncCross, expDistPairs, simDHdistPairs, maxVar, ACquant solqp
	inbreedFunc = NULL, # 
	withinFamInt = 1, #  
	setXint = NULL, # note that x is the cdf of a normal 
	skip = NULL,
# family parameters
	nFounder = 100,
	nNuclear = 100,
	nFam = 10,
	famSize = 50,
	ssd = FALSE,
	selF2 = FALSE,
	nF2 = 1,
# genetic parameters
	Vg = 1,
	updateVg = FALSE,
	h2 = c(0.3, 0.3, 0.3, 0.3),
# program parameters
	nYr = 30,
	selectTrials = c(0.5, 0.4, 0.3, 1/3),
	trialReps = c(1, 2, 3, 3),
	trialLocs = c(1, 2, 5, 5),
	cyclePerYr = 3,
	returnVDPtoRGSC = c(0, 0, 0, 0, 0), # default to rep(0, nTrial)?
	nGenOut = NULL,
	nGenInbr = NULL,
	phenoRGSC = 0,
# chromosome parameters
	nChrom = 10,
	nLoci = 100,
	nM = 100, # floor(c(10*(9:1), 4) / 4),
	nQTL = 100 # c(2*(9:1), 1),
)

type.convert(list(A = "TRUE", B = "10"))

stdArgNames <- names(defArgs)
defArgs <- getComArgs(defArgs)
altArgs <- names(defArgs)[!names(defArgs) %in% stdArgNames] 
# attach(defArgs)

# load libraries
library(AlphaSimR)
library(doMC)
library(txtplot) 

# source simulation function
source(defArgs$simFunc)

# print selection sizes
# nFam * famSize * c(1, cumprod(selectTrials))

if(!is.null(defArgs$projName)) defArgs$projName <- paste0(defArgs$projName, "/")
# make directory to store all results
simDir <- paste0(parDir, "/results/", defArgs$projName, defArgs$simName) 
system(paste0("simdir=", simDir, "
if [ ! -d $simdir ]; then
	mkdir -p $simdir
else 
	echo 'Directory already exists! Any results therein will be overwritten without further notice...'
fi"))


# set seed
if (!is.null(defArgs$seed)) set.seed(defArgs$seed)


# create simple founder pop, either simple or from loaded MaCS sim
if (defArgs$simpleFounder) {
	founderPop <- quickHaplo(nFounder, nChrom, nLoci, genLen = 1, ploidy = 2L, inbred = FALSE)
} else {
	load(defArgs$founderRData)
	if(!is.null(defArgs$founderSamples)) {
		load(defArgs$founderSamples)
		defArgs$founderSamples <- founderSamples
	}
	### NOTE, I MOVED THE FOUNDER SAMPLING INSIDE THE SIM FUNCTION SO THAT REPS WOULD BE ACROSS DIFFERENT FOUNDERS. (I.E. NO "FOUNDER" EFFECT)  
	# if (is.null(founderSamples)) founderPop <- founderPop[sample(1:nInd(founderPop), nFounder)]
	# sampleFounderPop(founderPop, size = defArgs$nFounder, n = 100)
}

# Setting Simulation Parameters
SP = SimParam$new(founderPop)
# SP$nThreads <- 40

# add trait 
SP$addTraitA(nQtlPerChr=defArgs$nQTL, mean = 0, var = defArgs$Vg)
# SP$addTraitAG(nQtlPerChr=nQTL, mean = 0, var = Vg, varGxE = Vgxe, varE = h2toVe(h2), corGxE = )

# add error variance as function of heritability
SP$setVarE(h2 = defArgs$founderh2)

# add genotype platform
SP$addSnpChip(defArgs$nM)

loci <- pullLoci(SP)
# Reduce(intersect, loci)


if(system("hostname", intern = TRUE) == "Bender") {
	setMKLthreads(1)
	registerDoMC(defArgs$nThreads)
} else {
	registerDoMC(defArgs$nThreads)
}

simrun <- foreach(k = 1:defArgs$nFounderPops) %:% foreach(r = 1:defArgs$reps, .errorhandling='pass') %dopar% do.call(sim, c(list(k = k, founderPop = founderPop, simParam = SP, paramL = defArgs, returnFunc = getPopStats), defArgs[altArgs]))
msg(0, "saving results in", paste0(simDir, "/", defArgs$simName, ".RData"))
save(simrun, SP, file = paste0(simDir, "/", defArgs$simName, ".RData"))
