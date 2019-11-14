# load breedingProgramR
library(breedingProgramR)
library(AlphaSimR)
# library(doMC)

# get parent directory
parDir <- getwd()
# set default arguments
userArgs <- list(
# simulation parameters
	maxIter = 1000L,
	reps = 1,
	nFounderPops = 1,
	lgen = 4,
	useTruth = 0, # make this an integer, 1 = use QTL, 2 = use true effects
	traditional = TRUE, # this selects out of VDP as parents, no RGSC
# founder parameters
	founderRData = "founderPop/testAlphaSimR1000SegSite.RData",
	founderSamples = "founderPop/founderSamples_nInd100_1000SegSite.RData",
	founderh2 = 0.3,
	founderBurnIn = 1,
	founderReps = 1,
	founderKeep = 4,
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
	nYr = 10,
	selectTrials = c(0.5, 0.4, 0.3, 1/3),
	trialReps = c(1, 2, 3, 3),
	trialLocs = c(1, 2, 5, 5),
	cyclePerYr = 3,
	returnVDPtoRGSC = c(0, 0, 0, 0, 0), # default to rep(0, nTrial)?
	nGenOut = NULL,
	nGenInbr = NULL,
	phenoRGSC = 0,
	separateTrain = FALSE,
# chromosome parameters
	nChrom = 10,
	nLoci = 1000,
	nM = 100, # floor(c(10*(9:1), 4) / 4),
	nQTL = 100 # c(2*(9:1), 1),
)

if(!interactive()){
	userArgs <- getComArgs(userArgs)
}

# create simple founder pop, either simple or from loaded MaCS sim
founderPop <- quickHaplo(userArgs$nFounder, userArgs$nChrom, userArgs$nLoci, genLen = 1, ploidy = 2L, inbred = FALSE)

# Setting Simulation Parameters
SP = SimParam$new(founderPop)
# SP$nThreads <- 40

# add trait 
SP$addTraitA(nQtlPerChr=userArgs$nQTL, mean = 0, var = userArgs$Vg)

# add error variance as function of heritability
SP$setVarE(h2 = userArgs$founderh2)

# add genotype platform
SP$addSnpChip(userArgs$nM)

loci <- pullLoci(SP)
# Reduce(intersect, loci)


trad <- simSingleTraitInbred(k = k, founderPop = founderPop, simParam = SP, paramL = userArgs, returnFunc = getPopStats)
names(trad)

plot(1:userArgs$nYr, trad$gvVDP$variety, type = "l")

userArgs
trad <- simSingleTraitInbred(k = k, founderPop = founderPop, simParam = SP, paramL = userArgs, returnFunc = getPopStats)

if(system("hostname", intern = TRUE) == "Bender") {
	print("Bender is great!")
	if(userArgs$nThreads > 1) setMKLthreads(1)
	registerDoMC(userArgs$nThreads)
} else {
	registerDoMC(userArgs$nThreads)
}

simrun <- foreach(k = 1:userArgs$nFounderPops) %:% foreach(r = 1:userArgs$reps, .errorhandling='pass') %dopar% do.call(simSingleTraitInbred, c(list(k = k, founderPop = founderPop, simParam = SP, paramL = userArgs, returnFunc = getPopStats), userArgs[altArgs]))
msg(0, "saving results in", paste0(simDir, userArgs$simName, ".RData"))
save(simrun, SP, file = paste0(simDir,  userArgs$simName, ".RData"))
