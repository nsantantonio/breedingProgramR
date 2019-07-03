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
	nFounderPops = 5,
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
	pullCycle = 1, 
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

testRun <- FALSE
if(testRun){
	# defArgs[["weight"]] <- 0.2
	# altArgs <- c(altArgs, "maxCrossPerParent", "weight")
	# defArgs[["GSfunc"]] <- RRBLUP2
	# defArgs[["maxIter"]] <- 500L
	defArgs0 <- defArgs
	run0 <- do.call(sim, c(list(founderPop = founderPop, simParam = SP, paramL = defArgs0, returnFunc = getPopStats), defArgs0[altArgs]))
	run0$gvRGSC
	# run0$gvVDP
	run0$varMean
	run0$VgRGSC

	defArgs0.1 <- defArgs
	defArgs0.1[["selFuncIn"]] <- truncSel 
	altArgs <- c("selFuncIn")

	run0.1 <- do.call(sim, c(list(founderPop = founderPop, simParam = SP, paramL = defArgs0, returnFunc = getPopStats), defArgs0.1[altArgs]))
	cor(run0$gvRGSC, run0.1$gvRGSC)
	# run0$gvVDP
	run0$varMean
	run0$VgRGSC

	# defArgs[["selFuncOut"]] <- solqpOut 
	defArgs[["selFuncIn"]] <- solqp 
	defArgs[["fthresh"]] <- 0.01
	# defArgs[["fthreshOut"]] <- 0.2
	altArgs <- c(altArgs, "selFuncIn","fthresh")
	# altArgs <- altArgs[[1]]

	run1 <- do.call(sim, c(list(founderPop = founderPop, simParam = SP, paramL = defArgs, returnFunc = getPopStats), defArgs[altArgs]))
	run1$gvRGSC
	run1$varMean
	run1$VgRGSC


	altArgs <- NULL
	defArgs[["selFuncOut"]] <- solqpOut 
	defArgs[["selFuncIn"]] <- solqp 
	defArgs[["pullCycle"]] <- 1
	defArgs[["fthresh"]] <- 0.01
	defArgs[["fthreshOut"]] <- 0.2
	altArgs <- c(altArgs, "selFuncOut","selFuncIn","pullCycle","fthresh","fthreshOut")
	# altArgs <- altArgs[[1]]

	run2 <- do.call(sim, c(list(founderPop = founderPop, simParam = SP, paramL = defArgs, returnFunc = getPopStats), defArgs[altArgs]))
	run2$gvRGSC
	run2$varMean
	run2$VgRGSC

	pdf("checkQPvsTrunc005In02Out.pdf")
	plot(1:defArgs$nYr, run0$varMean, type = "l", ylim = c(4, 35), main = "Variety Means")
 	lines(1:defArgs$nYr, run1$varMean, lty = 2)
	legend("topleft", lty = c(1, 2), legend = c("truncation", "quadProg"))
	
	plot(seq(0, defArgs$nYr, 1/3), run0$gvRGSC, type = "l", ylim = c(0, 35), main = "Recurrent Population Mean")
 	lines(seq(0, defArgs$nYr, 1/3), run1$gvRGSC, lty = 2)
	legend("topleft", lty = c(1, 2), legend = c("truncation", "quadProg"))
	
	plot(seq(0, defArgs$nYr, 1/3), run0$sdRGSC, type = "l", ylim = c(0, 1.2), main = "Recurrent Population Variance")
 	lines(seq(0, defArgs$nYr, 1/3), run1$sdRGSC, lty = 2)
	legend("topleft", lty = c(1, 2), legend = c("truncation", "quadProg"))
	
	dev.off()


	simL <- list()
	for(i in c(0.005, 0.01, 0.05, 0.1, 0.5)) {
	}
	gvRGSC <- lapply(simL, "[[", "gvRGSC")
	varMean <- lapply(simL, "[[", "varMean")
	Vg <- lapply(simL, "[[", "VgRGSC")
	sdRGSC <- lapply(simL, "[[", "sdRGSC")

	names(simL[[1]])


	defArgs[["fthresh"]] <- 0.1
	altArgs <- c(altArgs, "fthresh")
	altArgs <- altArgs[[1]]


	run2 <- do.call(sim, c(list(founderPop = founderPop, simParam = SP, paramL = defArgs, returnFunc = getPopStats), defArgs[altArgs]))
	run2$gvRGSC
	# run2$gvVDP
	run2$varMean

	pdf("qpSolByF_Variety.pdf")
	plot(1:defArgs$nYr, run0$varMean, type = "l", ylim = c(4, 30), main = "Variety Means")
	for(i in 1:length(varMean)) lines(1:defArgs$nYr, varMean[[i]], lty = i + 1)
	legend("topleft", lty = c(1, 1:length(varMean) + 1), legend = c("truncation", names(varMean)))
	dev.off()


	pdf("qpSolByF_RGSC.pdf")
	plot(seq(0, defArgs$nYr, 0.5), run0$gvRGSC, type = "l", ylim = c(0, 30), main = "Recurrent Population Mean")
	for(i in 1:length(varMean)) lines(seq(0, defArgs$nYr, 0.5), gvRGSC[[i]], lty = i + 1)
	legend("topleft", lty = c(1, 1:length(varMean) + 1), legend = c("truncation", names(varMean)))
	dev.off()

	pdf("qpSolByF_RGSC_Var.pdf")
	plot(seq(0, defArgs$nYr, 0.5), run0$VgRGSC, type = "l", ylim = c(0, 1.2), main = "Recurrent Population Variance")
	for(i in 1:length(varMean)) lines(seq(0, defArgs$nYr, 0.5), Vg[[i]], lty = i + 1)
	legend("topright", lty = c(1, 1:length(varMean) + 1), legend = c("truncation", names(varMean)))
	dev.off()

	pdf("qpSolByF_RGSC_Sd.pdf")
	plot(seq(0, defArgs$nYr, 0.5), run0$sdRGSC, type = "l", ylim = c(0, 1.2), main = "Recurrent Population Variance")
	for(i in 1:length(varMean)) lines(seq(0, defArgs$nYr, 0.5), sdRGSC[[i]], lty = i + 1)
	legend("topright", lty = c(1, 1:length(varMean) + 1), legend = c("truncation", names(varMean)))
	dev.off()


	reps <- 2
	setMKLthreads(10)
	registerDoMC(2)
	simrun <- foreach(r = 1:reps) %dopar% do.call(sim, c(list(r = r, founderPop = founderPop, simParam = SP, paramL = defArgs, returnFunc = getPopStats), defArgs[altArgs]))
}

if(system("hostname", intern = TRUE) == "Bender") {
	setMKLthreads(1)
	registerDoMC(defArgs$nThreads)
} else {
	registerDoMC(defArgs$nThreads)
}


# defArgs[["selFuncOut"]] <- solqpOut 
# defArgs[["selFuncIn"]] <- solqp 
# defArgs[["pullCycle"]] <- 1
# defArgs[["fthresh"]] <- 0.01
# defArgs[["fthreshOut"]] <- 0.2
# altArgs <- c(altArgs, "selFuncOut","selFuncIn","pullCycle","fthresh","fthreshOut")

# simrun <- foreach(r = 1:reps) %dopar% sim(founderPop, simParam = SP, paramL = defArgs, returnFunc = getPopStats)
# simrun <- foreach(k = 1:defArgs$nFounderPops) %:% foreach(r = 1:defArgs$reps) %do% do.call(sim, c(list(k = k, founderPop = founderPop, simParam = SP, paramL = defArgs, returnFunc = getPopStats), defArgs[altArgs]))
simrun <- foreach(k = 1:defArgs$nFounderPops) %:% foreach(r = 1:defArgs$reps, .errorhandling='pass') %dopar% do.call(sim, c(list(k = k, founderPop = founderPop, simParam = SP, paramL = defArgs, returnFunc = getPopStats), defArgs[altArgs]))
msg(0, "saving results in", paste0(simDir, "/", defArgs$simName, ".RData"))
save(simrun, SP, file = paste0(simDir, "/", defArgs$simName, ".RData"))
