#' simSingleTraitInbred function
#'
#' function to (do something)
#'
#' @param founderPop and object of Pop-class
#' @param paramL parameter list. See details for parameter values. 
#' @param simParam an object of an object of SimParam. Default is current global-env SP (see AlphaSimR manual for details)
#' @param returnFunc function, typically used to extract summary statistics from breeding program populations. Default is identity, which returns a list of populations in RCRS and VDP. 
#' @param verbose logical. Should information abou the breeding population be printed to stdout? Default is TRUE
#' @param checkParam logical. Should the parameter list be checked prior to running simulation? Default is TRUE
#' @param GSfunc function to estimate marker effects, must return an object of RRsol-class, SCAsol-class, or GCAsol-class. Default is NULL, which will use RRBLUP() until the number of individuals exceeds the number of markers * switchGSfunc, and will switch to RRBLUP2() for computational speed. 
#' @param switchGSfunc integer. Used to determine when to switch from RRBLUP to RRBLUP2 when the number of individuals exceeds number of markers * switchGSfunc. Default is 4.
#' @param k integer. Index of paramL$founderSamples to be used. Only used when paramL$founderSamples is a listof length > 1. 
#' @return [value]
#' @details
#' parameter list must contain the following elements to specify the breeding program. 
#'
#' paramL <- list(
#'
#' maxIter = 1000L, # integer. maximum iterations for GSfunc
#'
#' lgen = 4L, # integer. number of last generations to use for training GS model, greatly speeds up simulation to only use phenotypes from the last few generations to train the model
#'
#' useTruth = 0L,  # logical or integer. If set to 0, breeding program uses markers and phenotypes for training/selection. If set to 1, simulation uses QTL instead of markers, and if set to 2, simulation uses QTL and true breeding values instead of phenotypes for selection
#'
#' traditional = FALSE,  # logical or integer. A traditional program has no RCRS population, but instead selects the best individuals from the trial indicated by the argument (e.g. traditional = 2 selects new cadidates from the 2nd trial), and crosses them to an elite population that is updated with the best lines at the end of each year. 
#'
#' intAcross = 1, # numeric. proportion of families to select from when traditional > 0. Valid between 0 and 1.  
#'
#' intWithin = 0.2, # numeric. proportion of individuals within families to select from when traditional > 0. Valid between 0 and 1.
#'
#' founderSamples = NULL, # integer vector. Which individuals in the founder population should be sampled. Used when repeating simulations with different starting populations
#'
#' founderh2 = 0.3, # numeric. heritability of phenotypes of founder population to be used for intial selection from the founder population. Valid between 0 and 1.  
#'
#' founderBurnIn = 1L, # interger. Number of generations of random mating within founder population before breeding program begins.
#'
#' founderReps = 1L, # integer. number of replicates for founder population phenotypes
#'
#' founderKeep = 4L, # integer. Number of generations to include founder population phenotypes for training model.
#'
#' selectRCRS = 0.3, # numeric. selection intensity within recurrent population, only used for truncation selection. Valid between 0 and 1.  
#'
#' pullCycle = NULL,  # integer. Cycle (within generation) at which to branch material from RCRS, valid between 0 and cyclePerYr, where lower integers branch earlier.
#'
#' nProgenyPerCrossIn = 1L, # Number of progeny per cross within RCRS (check this!)
#'
#' nProgenyPerCrossOut = 1L, # Number of progeny per cross of material exiting RCRS for VDP
#'
#' useIn = "ebv",  # selection criterion for within RCRS. Typically wont deviate from 'ebv'
#'
#' useOut = "ebv",  # selection criterion for pulling material out of RCRS. Typically wont deviate from 'ebv'
#'
#' useInbreed = "ebv",  # selection criterion for lines within family during inbreeding. Only used when withinFamInt < 1. Typically wont deviate from 'ebv'
#'
#' useVDP = "pheno",  # selection criterion for advancing lines through the VDP. Typically 'pheno', but can also be 'ebv' if the user decides to use both phenotypes and genotypes for selection
#'
#' returnVDPcrit = "pheno",  # criterion used to determine best lines for return of material from VDP into RCRCS. Typically 'pheno', but can also be 'ebv'.
#'
#' selFuncOut = NULL,  # function that selects individuals out of RCRS. Takes in and outputs an argument of Pop-class.
#'
#' selFuncIn = NULL,  # function that selects individuals within the RCRS. Takes in and outputs an argument of Pop-class.
#'
#' selFuncVDP = NULL,  # function that selects individuals within the VDP. Takes in and outputs an argument of Pop-class.
#'
#' inbreedFunc = NULL,  # function that selects individuals during the inbreeding stage. Only used when withinFamInt < 1. Takes in and outputs an argument of Pop-class.
#'
#' withinFamInt = 1,  # intensity of selection within family before lines are entered into the VDP. forces famSize to be famSize / withinFamInt before the selection happens. 
#'
#' setXint = NULL,  # Set intensity used for quantile selection. [explain further]
#'
#' skip = NULL, # integer vector of length up to length(selectTrials) specifying which trials should be skipped, i.e. use 'ebv' for selection instead of phenotyping.
#'
#' nFounder = 100, # integer. Number of founders to be used.  
#'
#' nRCRS = 100, # integer. RCRS population size. [should modify to let change through time?]
#'
#' nFam = 10, # integer. number of families evaluated per year
#'
#' famSize = 50, # integer. family size. total size of first trial is then nFam * famSize
#'
#' ssd = FALSE, # logical. Should single seed decent be used? If FALSE, doubled haploids are produced. If TRUE, the population is selfed at every generation to produce one distinct individual 
#'
#' selF2 = FALSE, # logical. Should selection within the F2 generation be performed. [not sure this works anymore] 
#'
#' nF2 = 1, # integer. number of F2 individuals to be selected from if selF2 = TRUE. 
#'
#' Vg = 1, # numeric. Genetic variance of founder population.
#'
#' updateVg = FALSE, # logical. should Vg be held constant? If FALSE, heritability will change with decreasing Vg. If TRUE, heritability will be updated with current Vg (note this is not realistic). 
#'
#' h2 = c(0.3, 0.3, 0.3, 0.3), # numeric vector, length of number of trials, with plot level heritability of each trial. Valid between 0 and 1.
#'
#' nYr = 30, # integer. Number of years to run the breeding program. Phenotype trials will continue to nYr + number of trials
#'
#' selectTrials = c(0.5, 0.2, 0.5, 0.4), # numeric. Selection intensity within each trial. Valid between 0 and 1. Alternatively, positive integers can be provided to indicate the number of lines selected at each stage.
#'
#' trialReps = c(1, 2, 3, 3), # integer vector of length trials. Number of reps at each trial. 
#'
#' trialLocs = c(1, 2, 5, 5), # integer vector of length trials. Number of locations for each trial. 
#'
#' cyclePerYr = 3, # integer. number of cycles per year in RCRS. 
#'
#' returnVDPtoRCRS = c(0, 0, 0, 0, 0),  # numeric. proportion of the best lines from each trial to be returned to the RCRS 
#'
#' nGenOut = NULL, # integer. number of generations prior to current RCRS cycle for selecting lines for VDP. default is cyclePerYr - pullCycle
#'
#' nGenInbr = NULL, # integer. Number of generations of inbreeding lines destined for VDP. Only valid when inbreedFunc is NULL. Default is cyclePerYr.
#'
#' phenoRCRS = 0, # numeric. proportion of famSize to be used for phenotyping inbred lines out of RCRS population 
#'
#' separateTrain = FALSE # logical. Should only phenotyped individuals from RCRS be used to train model for selection within RCRS? Only used when phenoRCRS > 0.
#'
#')
#' @examples none
#' @export
simSingleTraitInbred <- function(founderPop, paramL, simParam = SP, returnFunc = identity, verbose = TRUE, checkParam = TRUE, GSfunc = NULL, switchGSfunc = 4, k = 1, ...){
# k = 1; paramL = userArgs; simParam <- SP; intAcross = 0.5; intWithin = 0.2; verbose = TRUE; checkParam = FALSE; GSfunc = RRBLUP; nGenOut = NULL; nGenInbr = NULL; returnFunc = getPopStats
	defArgs <- list(
		maxIter = 1000L,
		lgen = 4L,
		useTruth = 0L, 
		traditional = FALSE, 
		intAcross = 1,
		intWithin = 0.2,
		founderSamples = NULL,
		founderh2 = 0.3,
		founderBurnIn = 1L,
		founderReps = 1L,
		founderKeep = 4L,
		selectRCRS = 0.3,
		pullCycle = NULL, 
		nProgenyPerCrossIn = 1L,
		nProgenyPerCrossOut = 1L,
		useIn = "ebv", 
		useOut = "ebv", 
		useInbreed = "ebv", 
		useVDP = "pheno", 
		returnVDPcrit = "pheno", 
		selFuncOut = NULL, 
		selFuncIn = NULL, 
		selFuncVDP = NULL, 
		inbreedFunc = NULL, 
		withinFamInt = 1, 
		setXint = NULL, 
		skip = NULL,
		nFounder = 100,
		nRCRS = 100,
		nFam = 10,
		famSize = 50,
		ssd = FALSE,
		selF2 = FALSE,
		nF2 = 1,
		Vg = 1,
		updateVg = FALSE,
		h2 = c(0.3, 0.3, 0.3, 0.3),
		nYr = 30,
		selectTrials = c(0.5, 0.2, 0.5, 0.4),
		trialReps = c(1, 2, 3, 3),
		trialLocs = c(1, 2, 5, 5),
		cyclePerYr = 3,
		returnVDPtoRCRS = c(0, 0, 0, 0, 0), 
		nGenOut = NULL,
		nGenInbr = NULL,
		phenoRCRS = 0,
		separateTrain = FALSE
	)
	paramL <- argumentChange(defArgs, paramL)

	if (checkParam){
		paramNames <- c("maxIter", "lgen", "useTruth", "traditional", "intAcross", "intWithin", "founderSamples", "founderh2", "founderBurnIn", 
						"founderReps", "founderKeep", "selectRCRS", "pullCycle", "nProgenyPerCrossIn", "nProgenyPerCrossOut", "useIn", "useOut", 
						"useInbreed", "useVDP", "returnVDPcrit", "selFuncOut", "selFuncIn", "selFuncVDP", "inbreedFunc", "withinFamInt", "setXint",
						"skip", "nFounder", "nRCRS", "nFam", "famSize", "ssd", "selF2", "nF2", "Vg", "updateVg", "h2", "nYr", "selectTrials", "trialReps", 
						"trialLocs", "cyclePerYr", "returnVDPtoRCRS", "nGenOut", "nGenInbr", "phenoRCRS", "separateTrain")
		xtraParams <- names(paramL)[!names(paramL) %in% paramNames]
		if (length(xtraParams)) {
			cat("Unrecognized parameters have been provided to the 'paramL' argument. If these were meant to be passed to user defined functions, please provide them directly to simSingleTraitInbred().The following will be ignored:\n", paste(xtraParams, collapse = ", "), "\n")
		}
		if (!all(paramNames %in% names(paramL))) stop("not all parameters in 'paramL'! Please include all these parameters in parameter list:\n", paste0(paramNames, "\n"))
	}
	
	for (p in names(paramL)) assign(p, paramL[[p]])
	if(traditional > 0) msg(1, "Intensity across families:", intAcross, "Intensity within family:", intWithin)
	
	selFuncStop <- c(" is a list, but of the wrong structure. Please provide either a single function, a list of functions for for each year, each cycle, or a nested list of length 'nYr' with each element of length 'cyclePerYr'")
	if(phenoRCRS < 0 | phenoRCRS > 1) stop("'phenoRCRS' must be between 0 and 1 representing the proportion of the first VDP trial dedicated to phenotyping the RCRS!")

	if(is.list(selFuncIn)) {
		if (length(selFuncIn) %in% c(cyclePerYr, nYr, cyclePerYr * nYr)) {
			if(all(sapply(selFuncIn, class) == "list")){
				if (!all(unique(sapply(selFuncIn, length)) == cyclePerYr)) stop(paste0("selFuncIn", selFuncStop))
			}
			if(length(selFuncIn) == cyclePerYr) selFuncIn <- rep(list(selFuncIn), nYr)
			if(length(selFuncIn) == nYr) selFuncIn <- lapply(1:nYr, function(x) rep(list(selFuncIn[[x]]), cyclePerYr))
		} else {
			stop(paste0("selFuncIn", funcListStop))
		}
	} else {
		if(!is.null(selFuncIn)) {
			if(class(selFuncIn) == "function") selFuncIn <- rep(list(rep(list(selFuncIn), cyclePerYr)), nYr) else stop(paste0("selFuncIn", selFuncStop))
		}
	}

	if(is.function(selFuncOut)) selFuncOut <- rep(list(selFuncOut), nYr) else if(is.list(selFuncOut)) {if (length(selFuncOut) != nYr) stop(paste0("selFuncOut", funcListStop))}
	if(is.function(inbreedFunc)) inbreedFunc <- rep(list(inbreedFunc), nYr) else if(is.list(inbreedFunc)) {if (length(inbreedFunc) != nYr) stop(paste0("inbreedFunc", funcListStop))}

	if(!(length(selectRCRS) == 1 | length(selectRCRS) == nYr)) stop("'selectRCRS' must be of length 1 or nYr!")
	if(all(selectRCRS <= 1)) {
		selectRCRS <- nRCRS * selectRCRS
		if(selectRCRS %% 1 != 0) selectRCRS <- round(selectRCRS)
		deltaRCRS <- if(length(selectRCRS) == nYr) TRUE else FALSE
	}

	if(nFounder < nInd(founderPop)) {
		subFounder <-  if (is.null(founderSamples)) sample(1:nInd(founderPop), nFounder) else if(is.list(founderSamples)) founderSamples[[k]] else founderSamples
		founderPop <- founderPop[subFounder]
	}


	if (selF2 & cyclePerYr > 1) warning("Selection on F2 is being performed, and there is more than 1 GS cycle per year. You may want to reduce 'cyclePerYr' to 1")

	if (!all(selectTrials > 0) | (any(selectTrials < 1) & any(selectTrials > 1))) stop("'selectTrials' must have elements between 0 and 1 or positive integers > 0")
	if (!all(returnVDPtoRCRS >= 0) | (any(returnVDPtoRCRS < 1) & any(returnVDPtoRCRS > 1))) stop("'returnVDPtoRCRS' must have elements between 0 and 1 or positive integers")

	nI <- nFam * famSize
	if (all(selectTrials <= 1) & !all(selectTrials == 1)) selectTrials <- nI * cumprod(selectTrials) # note this does not allow all to be exactly 1

	if (any(selectTrials %% 1 != 0)){
		selectTrials <- round(selectTrials)
		actInt <- selectTrials / c(nI, selectTrials[-length(selectTrials)])
		msg(0, "NOTE: Selection intensities have been rounded to the nearest integer:\n", selectTrials, "\nThese correspond to selection intensities of:\n", actInt)
	}
	if(length(returnVDPtoRCRS) != length(selectTrials) + 1) stop("returnVDPtoRCRS must be of length(selectTrials) + 1 (for variety)!")
	if (all(returnVDPtoRCRS <= 1)) returnVDPtoRCRS <- returnVDPtoRCRS * c(nI, selectTrials) 

	if (withinFamInt > 1) withinFamInt <- famSize / ((nI + withinFamInt) / nFam)

	nTrial <- length(selectTrials)
	trials <- c(paste0("trial", 1:nTrial), "variety")
	if(is.function(selFuncVDP)) selFuncVDP <- rep(list(selFuncVDP), nYr + nTrial - 1) else if(is.list(selFuncVDP)) {if (length(selFuncVDP) != nYr + nTrial - 1) stop(paste0("selFuncVDP", funcListStop))}
	if (!is.null(skip)) skip <- trials[skip]

	if(traditional > 0) {
		paramL$cyclePerYr <- cyclePerYr <- 1
	}
	if(useTruth > 0) {
		msg(0, "NOTICE: using QTL as markers.")
		useTrue <- TRUE
	}
	if(useTruth > 1){
		msg(0, "NOTICE: using true marker/breeding values for all selection methods.")
		useOut <- "bv"
		useIn <- "bv"
		useInbreed <- "bv"
		gen0use <- "bv"
		estIntFunc <- bv
		returnVDPcrit <- "bv"
		useVDP <- "bv"
		useGS <- "bv"
		ebv <- bv
		Gvar <- varA
	} else {
		gen0use <- "pheno"
		estIntFunc <- pheno
		useGS <- "pheno"
		pullGenoFunc <- if(useTruth == 1) pullQtlGeno else pullSnpGeno 
		Gvar <- estVg
		useTrue <- FALSE
	}

	# define cycles
	if(is.null(pullCycle)) pullCycle <- cyclePerYr
	cycle <- 1:cyclePerYr

	# ignore GS models if only phenotypes used. 
	noGS <- if(useIn == "pheno" & useOut == "pheno" & withinFamInt == 1) TRUE else FALSE 
	if(noGS) msg(1, "noGS is true!")
	# initialize lists to store populations
	RCRS <- list()
	if(phenoRCRS > 0){
		RCRStoVDP <- list()
		RCRStoVDPmodel <- list()
	} 
	GSmodel <- list()
	predAcc <- list()
	predAcc[["VDP"]] <- list()
	intensity <- list()
	names(trials) <- trials
	VDP <- lapply(trials, function(x) list())
	if(traditional > 0) elite <- list() else elite <- NULL
	if(SP$isTrackPed) ped <- list()

	# initialize nuclear population, train GS model (necessary?), predict ebv (note the ebv's should be bad if markers are not exactly on QTL) 
	RCRS[[gen(0)]] <- newPop(founderPop)
	if(verbose) msg(1, "Founder population has genetic variance of:", popVar(getTrueBV(RCRS[[gen(0)]], simParam = simParam)))
	
	# burn in using random crosses
	printBurnin <- TRUE
	while (nFam > nInd(RCRS[[gen(0)]]) | founderBurnIn > 0) {
		if(nFam > nInd(RCRS[[gen(0)]]) & verbose) msg(1, "nFounder < nFam. Random mating to make nFam parents...")
		if(printBurnin & founderBurnIn > 0 & verbose) msg(1, "Running", founderBurnIn, "burn-in cycle(s) of random mating...")
		RCRS[[gen(0)]] <- selectCross(RCRS[[gen(0)]], nInd = nInd(RCRS[[gen(0)]]), use = "rand", simParam = simParam, nCrosses = max(nFam, nInd(RCRS[[gen(0)]]))) 
		founderBurnIn <- founderBurnIn - 1
		if(printBurnin) printBurnin <- FALSE
	}
	# phenotype founders
	RCRS[[gen(0)]] <- setPheno(RCRS[[gen(0)]], varE = h2toVe(founderh2, Vg), reps = founderReps)

	if(!noGS){	
		# train GS model and set EBV (after selfing if ssd)
		if (selF2) RCRS[[gen(0)]] <- self(RCRS[[gen(0)]], nProgeny = nF2, simParam = simParam)
		GSmodel[[gen(0)]] <- RRBLUP(RCRS[[gen(0)]], traits = 1, use = gen0use, snpChip = 1, simParam = simParam, useQtl = useTrue)
		# GSmodel[[gen(0)]] <- do.call(GSfunc, getArgs(GSfunc, pop = RCRS[[gen(0)]], traits = 1, use = gen0use, snpChip = 1, simParam = simParam, useQtl = useTrue, maxIter = 1000L, ...))
		RCRS[[gen(0)]] <- setEBV(RCRS[[gen(0)]], GSmodel[[gen(0)]], simParam = simParam)
		predAcc[["RCRS"]][[gen(0)]] <- getAcc(pop = RCRS[[gen(0)]], simParam = simParam)		
	}
	
	xInt <- if(is.null(setXint)) 1 - selectTrials[nTrial] / selectTrials[1] else setXint 
	if(verbose) msg(0, "Estimated selection intensity:", round(qnorm(xInt, sd = sqrt(Vg)), 3))

	# txtdensity(initSNPAlFreq)
	# pop <- RCRS[[gen(0)]]; GSfit <- GSmodel[[gen(0)]]
	# pop <- RCRS[[lastRCRSgen]]; GSfit <- GSmodel[[lastGSmodel]]

	pullRCRSgen <- gen(0)
	pullGSmodel <- gen(0)

	if(!exists("nElite")) nElite <- nFam 
	if(traditional > 0) elite[[gen(0)]] <- selectInd(RCRS[[gen(0)]], nInd = min(nInd(RCRS[[gen(0)]]), nElite), use = useIn)
	# sapply(RCRS, function(x){mean(gv(x))})
	# rlapply(VDP, function(x){mean(gv(x))}, level = 2, combine = c)
	# run program for nYr years
	for (i in 1:(nYr + nTrial - 1)) { 
	# for (i in 1:5) { 
		# i = 1
		lastRCRSgen <- names(RCRS)[length(RCRS)]
		lastGSmodel <- if (i <= nYr) gen(i-1) else gen(nYr)
		selectRCRSi <- if(deltaRCRS) selectRCRS[i] else selectRCRS[1]
		# selectTrialsi <- selectTrials
		# nProgenyPerCrossIni <- if(identical(expDistPairs, selFuncIn)) nRCRS / selectRCRSi  * nProgenyPerCrossIn else nProgenyPerCrossIn

		# update elite pop from previous year
		if (traditional > 0){
			if(i == 1){
				elite[[gen(i)]] <- elite[[gen(0)]]
				# if(nElite > nFam) elite[[gen(i)]] <- elite[[gen(i)]][sample(nElite, nFam)]
			} else if(i == 2 & nElite > nFam){
				elite[[gen(i)]] <- elite[[gen(0)]][!elite[[gen(0)]]@id %in% elite[[gen(1)]]@id]
			} else if(i > 1) {
				mostAdvancedTrial <- tail(which(sapply(VDP[trials[-length(trials)]], length) > 0), 1)
				elSel <- if(i <= traditional) VDP[mostAdvancedTrial] else lapply(VDP, tail, 1)[trials[{traditional + 1}:mostAdvancedTrial]] 
				elite[[gen(i)]] <- selectInd(mergePopsRec(elSel), nInd = nElite, use = useIn)
				rm(elSel)
			}
			if(nInd(elite[[gen(i)]]) > nFam) elite[[gen(i)]] <- elite[[gen(i)]][sample(nInd(elite[[gen(i)]]), nFam)]
		}

		if (i <= nYr){ 
			if (verbose) msg(0, "Year: ", i)

			# predict latest RCRS with updated GS model 
			if (!noGS & i > 1) RCRS[[lastRCRSgen]] <- setEBV(RCRS[[lastRCRSgen]], GSmodel[[lastGSmodel]], simParam = simParam)

			# estimate VDP selection intensity from previous generations
			if (i > nTrial) {
				intensity[[gen(i)]] <- estIntensity(VDP, i, nT = nTrial, start = "trial1", end = "variety", estFunc = estIntFunc, Gvar = Gvar)
				if(is.null(setXint)) xInt <- pnorm(intensity[[gen(i)]]) 
				if(verbose) msg(1, "Realized selection intensity:", round(intensity[[gen(i)]], 3))
				if(verbose) msg(1, "best variety:", round(max(gv(VDP[["variety"]][[gen(i - nTrial)]])), 3))
			}
			
			# use selected lines from last trial to make new crosses. Not sure this timeline is realistic...
			if(verbose) msg(1, "Selecting lines for entry into the VDP...")
			if(i == 1 | is.null(selFuncOut)){

				selToP <- if(nFam < nInd(RCRS[[lastRCRSgen]])) selectInd(RCRS[[lastRCRSgen]], nInd = nFam, trait = 1, use = useOut) else RCRS[[lastRCRSgen]]
			} else {
				pullGSmodel <- gen(i - 1) # note, this should actually be defaulting to lastGSmodel
				# select out of RCRS
				# cyclePerYr - pullCycle + delay * cyclePerYr
				if(is.null(nGenOut)) nGenOut <- cyclePerYr - pullCycle # if longer you need to count!
				# GSmodelOut <- if(separateTrain) RCRStoVDPmodel[c(pullGSmodel, lastGSmodel)] else GSmodel[c(pullGSmodel, lastGSmodel)]
				GSmodelOut <- if(separateTrain) RCRStoVDPmodel[[pullGSmodel]] else GSmodel[[lastGSmodel]]
				selToP <- do.call(selFuncOut[[i]], getArgs(selFuncOut[[i]], pop = RCRS[[pullRCRSgen]], GSfit = GSmodelOut, nSel = nFam, nGenOut = nGenOut, nGenThisYr = cyclePerYr - pullCycle, 
												  trait = 1, use = useOut, quant = xInt, nProgeny = nProgenyPerCrossOut, Gvar = Gvar, simParam = simParam, ...))
												  # trait = 1, use = useOut, quant = xInt, nProgeny = nProgenyPerCrossOut, Gvar = Gvar, simParam = simParam, fthreshOut = 0.2))
				# double check this is the correct GS model!!!! Maybe needs to be pullGSmodel???
			}
			# check GP accuracy for material to VDP
			if(!noGS) {
				selToP <- setEBV(selToP, GSmodel[[lastGSmodel]], simParam = simParam) #
				predAcc[["RCRSout"]][[gen(i)]] <- getAcc(pop = selToP, simParam = simParam)
			}

			if(verbose) msg(1, "VDP input Vg:", round(varA(selToP), 6))
			if(verbose) msg(1, "VDP input pop mean:", round(mean(gv(selToP)), 6))
			famSizei <- round(nFam / nInd(selToP) * famSize) 

			# Save some proportion of RCRS to put into VDP to update training models. 
			if(phenoRCRS > 0) {
				RCRStoVDPgen <- lastRCRSgen
				# see how many lines from RCRS you can sample, then add back to famSizei
				newFamSizei <- round(famSizei * (1 - phenoRCRS))
				availPlots <- (nInd(selToP) * famSizei) - (nInd(selToP) * newFamSizei)
				RCRSfamSize <- availPlots / nInd(RCRS[[RCRStoVDPgen]])
				if(RCRSfamSize < 1) {
					nPhRS <- RCRSfamSize * nInd(RCRS[[RCRStoVDPgen]]) 
					RCRSfamSize <- 1
					whichRCRStoVDP <- sample(nInd(RCRS[[RCRStoVDPgen]]), nPhRS)
				} else {
					nPhRS <- nInd(RCRS[[RCRStoVDPgen]]) 
					RCRSfamSize <- floor(RCRSfamSize)
					whichRCRStoVDP <- 1:nInd(RCRS[[RCRStoVDPgen]])
				}
				availPlots <- availPlots - RCRSfamSize * nPhRS
				famSizei <- newFamSizei + floor(availPlots / nInd(selToP))
				
				if(ssd){
					RCRStoVDP[[gen(i)]] <- self(RCRS[[RCRStoVDPgen]][whichRCRStoVDP], nProgeny = RCRSfamSize)
					for(i in 1:(cyclePerYr-1)) RCRStoVDP <- self(RCRS[[RCRStoVDPgen]], nProgeny = 1)
				} else {
					RCRStoVDP[[gen(i)]] <- makeDH(RCRS[[RCRStoVDPgen]][whichRCRStoVDP], nDH = RCRSfamSize)
				}
				if(verbose) {
					msg(1, nInd(RCRStoVDP[[gen(i)]]), "inbred lines from", nPhRS, "RCRS included in VDP for making crosses")
					msg(1, "Family size reduced from", famSize, "to", famSizei)
				}
				# selectTrialsi[1] <-  selectTrialsi[1] # dont update, allows for no selection in first
			}

			if(traditional > 0 & i > 1) {
				# if(verbose) msg(1, nInd(selToP), "lines selected out of VDP", VDPsel, "for making crosses")
				if(verbose) msg(1, nInd(selToP), "crosses made out of VDP", VDPsel)
			} else if(i > 1) {
				if(verbose) msg(1, nInd(selToP), "crosses selected out of", nInd(RCRS[[pullRCRSgen]]), "RCRS population for VDP")
			} else {
				if(verbose) msg(1, nInd(selToP), "crosses selected out of", nInd(RCRS[[gen(0)]]), "founder population for VDP")
			}
			
			# make DH families or self
			if(is.null(inbreedFunc)) {
				nProgPerFam <- famSizei / withinFamInt
				if (nProgPerFam %% 1 != 0) {
					nProgPerFam <- round(nProgPerFam)
					nSelToTrial <- round(nProgPerFam * withinFamInt)
					msg(0, "NOTICE: Selection intensities within familiy have been rounded to the nearest integer resulting in", nSelToTrial, "progeny per family selected from", nProg, "progeny per family")
				}

				# HERE IT IS!? for some reason mean(pheno(makeDH(selToP, nDH = nProgPerFam))) << mean(pheno(selToP))
				VDP[[trials[1]]][[gen(i)]] <- if (ssd) self(selToP, nProgeny = nProgPerFam) else makeDH(selToP, nDH = nProgPerFam)
				
				msg(1, "DH parent mean:", round(mean(pheno(selToP)), 6), "DH mean:", round(mean(pheno(VDP[[trials[1]]][[gen(i)]])), 6))
				msg(1, "DH parent GV mean:", round(mean(gv(selToP)), 6), "DH GV mean:", round(mean(gv(VDP[[trials[1]]][[gen(i)]])), 6))

				#select within family if intensity < 1
				if (withinFamInt < 1) {
					VDP[[trials[1]]][[gen(i)]] <- setEBV(VDP[[trials[1]]][[gen(i)]], GSmodel[[lastGSmodel]], simParam = simParam)
					fams <- split(1:(nProgPerFam * nInd(selToP)), rep(1:nInd(selToP), each = nProgPerFam))
					faml <- list()
					for (j in names(fams)){
						faml[[j]] <- selectInd(VDP[[trials[1]]][[gen(i)]][fams[[j]]], nInd = famSizei, trait = 1, use = useOut) 
					}
					VDP[[trials[1]]][[gen(i)]] <- mergePops(faml)		
				}
			} else {
				if(is.null(nGenInbr)) nGenInbr <- cyclePerYr  
				VDP[[trials[1]]][[gen(i)]] <- do.call(inbreedFunc[[i]], getArgs(inbreedFunc[[i]], pop = selToP, GSfit = GSmodel[[lastGSmodel]], # use last GS model here because selectOut will burn through the rest of the year. Need to make flexible if not. 
					trait = 1, use = useInbreed, int = withinFamInt, ssd = ssd, nProgeny = famSizei, nGenInbr = nGenInbr, simParam = simParam, ...))
					# trait = 1, use = useInbreed, int = withinFamInt, ssd = ssd, nProgeny = famSizei, nGenInbr = nGenInbr, simParam = simParam))
			}

			if(!noGS) VDP[[trials[1]]][[gen(i)]] <- setEBV(VDP[[trials[1]]][[gen(i)]], GSmodel[[lastGSmodel]], simParam = simParam) #
			predAcc[["VDPin"]][[gen(i)]] <- if(noGS) NULL else getAcc(pop = VDP[[trials[1]]][[gen(i)]], simParam = simParam)

			msg(1, "best new line:", round(max(gv(VDP[[trials[1]]][[gen(i)]])), 6))
		}
		# get generation indices
		genI <- tail(1:i, min(nTrial, i))
		genI <- genI[genI <= nYr]
		genBack <- abs(genI - i) + 1
		index = 1:length(genI)

		# phenotype, or predict if skipped
		for (g in index) {
			Vgi <- if (updateVg) varG(VDP[[trials[genBack[g]]]][[gen(genI[g])]])[[1]] else Vg
			if (!trials[genBack[g]] %in% skip) VDP[[trials[genBack[g]]]][[gen(genI[g])]] <- setPheno(VDP[[trials[genBack[g]]]][[gen(genI[g])]], varE = h2toVe(h2[genBack[g]], Vgi), reps = trialReps[genBack[g]] * trialLocs[genBack[g]])
			if(phenoRCRS > 0 & i > 1 & g == 1 & i <= nYr) RCRStoVDP[[gen(i)]] <- setPheno(RCRStoVDP[[gen(i)]], varE = h2toVe(h2[genBack[g]], Vgi), reps = trialReps[genBack[g]] * trialLocs[genBack[g]])
		}

		if (i <= nYr){
			if(SP$isTrackPed) ped[[gen(i)]] <- cbind(VDP[[trials[1]]][[gen(i)]]@mother, VDP[[trials[1]]][[gen(i)]]@father, VDP[[trials[1]]][[gen(i)]]@id)
			# if(traditional) 
			for (j in cycle){
				jp <- which(cycle == j)
				if(traditional > 0) {
					whichTrials <- which(sapply(VDP, length) > 0)
					VDPsel <- if(is.logical(traditional) | max(whichTrials) < traditional) tail(names(VDP)[whichTrials], 1) else paste0("trial", traditional)
					oriGen <- if(VDPsel == "variety") i - length(trials) + 1 else i - as.numeric(gsub("[A-z]", "", VDPsel)) + 1 # this is hacky, but it works...
					# oriGen <- i - min(i, traditional) + 1
					selPop <- VDP[[VDPsel]][[oriGen]] 
					families <- split(VDP[[trials[1]]][[gen(oriGen)]]@id, rep(1:nFam, each = famSizei)) # this assumes famSizei is a single integer, need to allow for diff number ind per fam!!!!!!!!!!!!!!!
					if(verbose) msg(1, "Individuals selected out of", VDPsel, "from generation", oriGen, "for crossing...")
				
				} else {
					RCRS[[gen(j-1)]] <- setEBV(RCRS[[gen(j-1)]], GSmodel[[lastGSmodel]], simParam = simParam)
					selPop <- RCRS[[gen(j-1)]]
					predAcc[["RCRS"]][[gen(j-1)]] <- getAcc(selPop)
				}
				# run GS model to cycle through RCRS for year i
 				if(traditional > 0) {
					RCRS[[gen(j)]] <- do.call(tradSelCross2, getArgs(tradSelCross2, pop = selPop, elitepop = elite[[gen(i)]], families = families, nFam = nFam, famSize = famSizei, use = useIn, trait = 1, simParam = simParam, 
						nCrosses = nFam, nProgeny = nProgenyPerCrossIn, verbose = verbose, ...))
						# nCrosses = nFam, nProgeny = nProgenyPerCrossIn, verbose = verbose))
				} else if(is.null(selFuncIn)){
					RCRS[[gen(j)]] <- selectCross(pop = selPop, nInd = min(selectRCRSi, nInd(selPop)), use = useIn,  trait = 1, simParam = simParam, nCrosses = nRCRS, nProgeny = nProgenyPerCrossIn)
				} else {
					RCRS[[gen(j)]] <- do.call(selFuncIn[[i]][[jp]], getArgs(selFuncIn[[i]][[jp]], nSel = selectRCRSi, pop = selPop, GSfit = GSmodel[[lastGSmodel]],
					trait = 1,  use = useIn, nCrosses = nRCRS, nProgeny = nProgenyPerCrossIn, quant = xInt, verbose = verbose, Gvar = Gvar, simParam = simParam, ...))
					# trait = 1,  use = useIn, nCrosses = nRCRS, nProgeny = nProgenyPerCrossIn, quant = xInt, verbose = verbose, Gvar = Gvar, simParam = simParam, fthresh = 0.01))
					if(nInd(RCRS[[gen(j)]]) != nRCRS) msg(2, "Only", nInd(RCRS[[gen(j)]]), "crosses made in RCRS...")
				}
				# would be good to be able to select within f2 family if f2 > 1
				if (selF2) RCRS[[gen(j)]] <- self(RCRS[[gen(j)]], nProgeny = nF2, simParam = simParam)
				# if(jp == pullCycle) {
				# 	pullRCRSgen <- gen(j) # should I push material out earlier?
				# }
				if(verbose) msg(1, "RCRS Vg:", round(varA(RCRS[[gen(j)]]), 6))
				if(verbose) msg(1, "RCRS Pop Mean:", round(mean(RCRS[[gen(j)]]@gv), 6))
			}
			pullRCRSgen <- if(pullCycle > 0) gen(cycle[[pullCycle]]) else lastRCRSgen
			# update GScycle number
			cycle <- cycle + cyclePerYr
			
			if(!noGS){
				# so this removes the selections from the training population. 
				trnSet <- lapply(VDP[trials[!grepl("variety", trials)]], function(x) x[names(x) %in% gen(max(1, i-max(1, lgen)):i)])
				trnSet <- trnSet[sapply(trnSet, length) > 0]
				names(trnSet)
				# include founders in prediction model
				if(founderKeep >= i) trnSet[["founder"]] <- RCRS[[gen(0)]]

				# Inlcude RCRS individuals into training
				if(phenoRCRS > 0) {
					if(separateTrain) {
						RCRStoVDPtrian <- RCRStoVDP
						# RCRStoVDPtrian <- RCRStoVDP[names(RCRStoVDP) %in% gen(max(1, i-max(1, lgen)):i)]  
						if(founderKeep >= i) RCRStoVDPtrian["founder"] <- RCRS[gen(0)]
						RCRStoVDPtrian <- mergePopsRec(RCRStoVDPtrian) 
						msg(1, "Training set for RCRSout branch has", RCRStoVDPtrian@nInd, "individuals...")	
						if(!is.null(GSfunc)){
							RCRStoVDPmodel[[gen(i)]] <- do.call(GSfunc, getArgs(GSfunc, pop = RCRStoVDPtrian, traits = 1, use = useGS, snpChip = 1, simParam=simParam, useQtl = useTrue, maxIter = maxIter))#, ...))
						} else {
							if(nInd(RCRStoVDPtrian) > simParam$snpChips[[1]]@nLoci * switchGSfunc) {
								RCRStoVDPmodel[[gen(i)]] <- do.call(RRBLUP2, getArgs(RRBLUP2, pop = RCRStoVDPtrian, traits = 1, use = useGS, snpChip = 1, simParam=simParam, useQtl = useTrue, maxIter = maxIter))#, ...))
							} else {
								RCRStoVDPmodel[[gen(i)]] <- do.call(RRBLUP, getArgs(RRBLUP, pop = RCRStoVDPtrian, traits = 1, use = useGS, snpChip = 1, simParam=simParam, useQtl = useTrue, maxIter = maxIter))#, ...))
							}
						}
					} else {
						trnSet[["RCRStoVDP"]] <- RCRStoVDP[names(RCRStoVDP) %in% gen(max(1, i-max(1, lgen)):i)]
					}
				}
				
				# concatenate training set and train GS model
		 		train <- mergePopsRec(trnSet) 
				msg(1, "Training set has", train@nInd, "individuals...")	
				if(genParam(train)$varG == 0) {
					GSmodel[[gen(i)]] <- GSmodel[[gen(i - 1)]]
				} else {
					if(!is.null(GSfunc)){
						GSmodel[[gen(i)]] <- do.call(GSfunc, getArgs(GSfunc, pop = train, traits = 1, use = useGS, snpChip = 1, simParam=simParam, useQtl = useTrue, maxIter = maxIter))#, ...))
					} else {
						if(nInd(train) > simParam$snpChips[[1]]@nLoci * switchGSfunc) {
							GSmodel[[gen(i)]] <- do.call(RRBLUP2, getArgs(RRBLUP2, pop = train, traits = 1, use = useGS, snpChip = 1, simParam=simParam, useQtl = useTrue, maxIter = maxIter))#, ...))
						} else {
							GSmodel[[gen(i)]] <- do.call(RRBLUP, getArgs(RRBLUP, pop = train, traits = 1, use = useGS, snpChip = 1, simParam=simParam, useQtl = useTrue, maxIter = maxIter))#, ...))
						}
					}
				}
			}
		}
		# check if final year
		if (i - nYr == 1) msg(0, "Final year reached, selecting on phenotypes / ebv trained with last year training set ...")
		
		for (g in index) {
			GSmodelVDP <- if (trials[genBack[g]] %in% skip) lastGSmodel else gen(min(i, nYr))
			if(!noGS) {
				VDP[[trials[genBack[g]]]][[gen(genI[g])]] <- setEBV(VDP[[trials[genBack[g]]]][[gen(genI[g])]], GSmodel[[GSmodelVDP]], simParam = simParam)
				predAcc[["VDP"]][[trials[[genBack[g]]]]][gen(genI[g])] <- getAcc(pop = VDP[[trials[genBack[g]]]][[gen(genI[g])]])
			}
			# select based on ebv or phenotype
			sel <- if (trials[genBack[g]] %in% skip) "ebv" else  useVDP
			if (i - genI[g] < nTrial) {
				if(is.null(selFuncVDP)){
					VDP[[trials[genBack[g] + 1]]][[gen(genI[g])]] <- selectInd(VDP[[trials[genBack[g]]]][[gen(genI[g])]], nInd = min(selectTrials[genBack[g]], nInd(VDP[[trials[genBack[g]]]][[gen(genI[g])]])), trait = 1, use = sel, returnPop = TRUE)
				} else {
					VDP[[trials[genBack[g] + 1]]][[gen(genI[g])]] <- do.call(selFuncVDP[[i]], getArgs(selFuncVDP[[i]], nSel = min(selectTrials[genBack[g]], nInd(VDP[[trials[genBack[g]]]][[gen(genI[g])]])), 
																									  pop = VDP[[trials[genBack[g]]]][[gen(genI[g])]], GSfit = GSmodel[[lastGSmodel]], trait = 1, use = sel, 
																									  returnPop = TRUE, verbose = verbose, Gvar = Gvar, simParam = simParam, ...))
				}
			}
			if (ssd) VDP[[trials[genBack[g] + 1]]][[gen(genI[g])]] <- self(VDP[[trials[genBack[g] + 1]]][[gen(genI[g])]])
		}

		if (i <= nYr){
			# return lines from VDP into the RCRS 
			if (any(returnVDPtoRCRS > 0)){
				returnToRCRS <- genBack %in% which(returnVDPtoRCRS > 0)
				if (sum(returnToRCRS) > 0){	
					addToRCRS <- list()
					for (g in index[returnToRCRS]) {
						addToRCRS[[gen(g)]] <- selectInd(VDP[[trials[genBack[g]]]][[gen(genI[g])]], nInd = returnVDPtoRCRS[genBack[g]], trait = 1, use = returnVDPcrit, returnPop = TRUE) 
					}
					if (length(addToRCRS) > 0){
						addToRCRS <- mergePopsRec(addToRCRS)
						RCRS[[gen(cycle[1] - 1)]] <- c(RCRS[[gen(cycle[1] - 1)]], addToRCRS)
					}
				}
			}
		}
		if(verbose) cat("\n")
	}

	rL <- list(SP = SP, paramL = paramL, RCRS = RCRS, VDP = VDP, GSmodel = GSmodel, predAcc = predAcc, selIntensity = intensity)
	do.call(returnFunc, getArgs(returnFunc, resultL = rL, ...))
}
