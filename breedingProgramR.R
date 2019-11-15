getComArgs <- function(defaultArgs = NULL) {
  defaults <- !is.null(defaultArgs)
  args <- commandArgs(TRUE)
  isAssn <- grepl("=", args)
  userArgs <- args[isAssn]
  needEval <- grepl("\\(|\\)|\\:|\\'", userArgs) 
  argSplit <- strsplit(userArgs, "=")
  argList <- lapply(argSplit, "[[", 2)
  names(argList) <- lapply(argSplit, "[[", 1)
  argList[needEval] <- lapply(argList[needEval], function(x) eval(parse(text = x)))
  argList[!needEval] <- lapply(argList[!needEval], function(x) strsplit(x, ",")[[1]])
  argList[!needEval] <- type.convert(argList[!needEval], as.is = TRUE)
  print(argList)
  if (defaults){
    defaultArgs[names(argList)] <- argList
    return(defaultArgs)
  } else {
    return(argList)
  }
}

getArgs <- function(f, ...) { list(...)[names(list(...)) %in% names(formals(f))] }

argumentChange <- function(defaultArgs, userArgs){
	defaultArgs[names(userArgs)] <- userArgs
	return(defaultArgs)
}

msg <- function(n, ...) { cat(do.call(paste, c(rep(list("    "), n), ..., collapse = "", sep = " ")), "\n", sep = "") }


simSingleTraitInbred <- function(founderPop, paramL, simParam = SP, returnFunc = identity, verbose = TRUE, checkParam = FALSE, GSfunc = NULL, switchGSfunc = 4, ...){
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
		inbreedFunc = NULL, 
		withinFamInt = 1, 
		setXint = NULL, 
		skip = NULL,
		nFounder = 100,
		nNuclear = 100,
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
		paramNames <- c("seed", "nThreads", "projName", "simName",  "maxIter", "reps", "nFounderPops", "lgen", "useTruth", "traditional", 
					      "founderh2", "simpleFounder", "founderBurnIn", "founderReps", "founderKeep", "selectRCRS", 
					    "pullCycle", "nProgenyPerCrossIn", "nProgenyPerCrossOut", "useIn", "useOut", "useInbreed", "useVDP", "returnVDPcrit", 
					    "selFuncOut", "selFuncIn", "inbreedFunc", "withinFamInt", "setXint", "skip", "nFounder", "nNuclear", "nFam", "famSize", 
					    "ssd", "selF2", "nF2", "Vg", "updateVg", "h2", "nYr", "selectTrials", "trialReps", "trialLocs", "cyclePerYr", "returnVDPtoRCRS", 
					    "nGenOut", "nGenInbr", "phenoRCRS", "separateTrain", "nChrom", "nLoci", "nM", "nQTL")
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
		selectRCRS <- nNuclear * selectRCRS
		if(selectRCRS %% 1 != 0) selectRCRS <- round(selectRCRS)
		deltaRCRS <- if(length(selectRCRS) == nYr) TRUE else FALSE
	}

	if(nFounder < nInd(founderPop)) {
		subFounder <-  if (is.null(founderSamples)) sample(1:nInd(founderPop), nFounder) else founderSamples[[k]]
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
	if (!is.null(skip)) skip <- trials[skip]

	if(traditional > 0) {
		cyclePerYr <- 1
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
		# nProgenyPerCrossIni <- if(identical(expDistPairs, selFuncIn)) nNuclear / selectRCRSi  * nProgenyPerCrossIn else nProgenyPerCrossIn

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
					nRCRS <- RCRSfamSize * nInd(RCRS[[RCRStoVDPgen]]) 
					RCRSfamSize <- 1
					whichRCRStoVDP <- sample(nInd(RCRS[[RCRStoVDPgen]]), nRCRS)
				} else {
					nRCRS <- nInd(RCRS[[RCRStoVDPgen]]) 
					RCRSfamSize <- floor(RCRSfamSize)
					whichRCRStoVDP <- 1:nInd(RCRS[[RCRStoVDPgen]])
				}
				availPlots <- availPlots - RCRSfamSize * nRCRS
				famSizei <- newFamSizei + floor(availPlots / nInd(selToP))
				
				if(ssd){
					RCRStoVDP[[gen(i)]] <- self(RCRS[[RCRStoVDPgen]][whichRCRStoVDP], nProgeny = RCRSfamSize)
					for(i in 1:(cyclePerYr-1)) RCRStoVDP <- self(RCRS[[RCRStoVDPgen]], nProgeny = 1)
				} else {
					RCRStoVDP[[gen(i)]] <- makeDH(RCRS[[RCRStoVDPgen]][whichRCRStoVDP], nDH = RCRSfamSize)
				}
				if(verbose) {
					msg(1, nInd(RCRStoVDP[[gen(i)]]), "inbred lines from", nRCRS, "RCRS included in VDP for making crosses")
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
					RCRS[[gen(j)]] <- selectCross(pop = selPop, nInd = min(selectRCRSi, nInd(selPop)), use = useIn,  trait = 1, simParam = simParam, nCrosses = nNuclear, nProgeny = nProgenyPerCrossIn)
				} else {
					RCRS[[gen(j)]] <- do.call(selFuncIn[[i]][[jp]], getArgs(selFuncIn[[i]][[jp]], nSel = selectRCRSi, pop = selPop, GSfit = GSmodel[[lastGSmodel]],
					trait = 1,  use = useIn, nCrosses = nNuclear, nProgeny = nProgenyPerCrossIn, quant = xInt, verbose = verbose, Gvar = Gvar, simParam = simParam, ...))
					# trait = 1,  use = useIn, nCrosses = nNuclear, nProgeny = nProgenyPerCrossIn, quant = xInt, verbose = verbose, Gvar = Gvar, simParam = simParam, fthresh = 0.01))
					if(nInd(RCRS[[gen(j)]]) != nNuclear) msg(2, "Only", nInd(RCRS[[gen(j)]]), "crosses made in RCRS...")
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
			if (i - genI[g] < nTrial) VDP[[trials[genBack[g] + 1]]][[gen(genI[g])]] <- selectInd(VDP[[trials[genBack[g]]]][[gen(genI[g])]], nInd = min(selectTrials[genBack[g]], nInd(VDP[[trials[genBack[g]]]][[gen(genI[g])]])), trait = 1, use = sel, returnPop = TRUE)
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


sampleFounderPop <- function(founderPop, size, n, savePop = FALSE, seed = NULL) {
	if(!is.null(seed)) set.seed(seed)
	popSamples <- list()
	for(i in 1:n){
		samplei <- sample(1:nInd(founderPop), nFounder)
		popSamples[[i]] <- if(savePop) founderPop[samplei] else samplei 
	}
	popSamples
}


 # nFam * famSize

inbreedFunc <- function(pop, GSfit, nProgPerFam, ...) {
	pop <- if(ssh) self(pop, nProgeny = nProgPerFam) else makeDH(pop, nDH = nProgPerFam)
	if (withinFamInt < 1) {
		pop <- setEBV(pop, GSmodel[[lastGSmodel]], simParam = simParam)
		fams <- split(1:(nProgPerFam * nFam), rep(1:nFam, each = nProgPerFam))
		faml <- list()
		for (j in names(fams)){
			faml[[j]] <- selectInd(pop[fams[[j]]], nInd = famSize, trait = 1, use = selectOut) 
		}
		pop <- mergePops(faml)		
	}
	# # print mean genotypic value of DH 
	# if (verbose) print(sapply(VDP[[trials[1]]], function(x) mean(gv(x))))
	if(nInd(pop) != nFam * famSize) stop("selToP is wrong...")
}
# popList <- list(1, 2, list(3, 4, list(5, 6, list(7))), list(8, list(9, 10)))

# cR <- function(popList){
#   unlist(lapply(popList, function(x) if (is.list(x)) cR(x) else x), recursive = FALSE)
# }
# cR(popList)

mergePopsRec <- function(popList) { mergePops(lapply(popList, function(x) {if (is.list(x)) mergePopsRec(x) else x})) }

getAF <- function(pop, pullGeno = pullSnpGeno) { colMeans(pullGeno(pop)) / 2 }
getF <- function(pop, pullGeno = pullSnpGeno) {1 + (0.5 - rowMeans(pullGeno(pop) == 1)) * 2} # not sue this is right...
h2toVe <- function(h2, Vg = 1) { Vg * (1-h2) / h2 }
gen <- function(i) { paste0("gen", i) }
logSelInd <- function(pop, sel) { pop@id %in% sel   }
rSel <- function(sel) { Reduce("&", sel) }
maxBv <- function(simParam, traits = 1) { sapply(simParam$traits[traits], function(x) { sum(abs(x@addEff)) }) }
sdUnCor <- function(x) { sqrt(mean(x^2) - mean(x)^2) }
getRRh2 <- function(rrFit) { solve(pop0pred@Vu + pop0pred@Ve) %*% pop0pred@Vu }

tradSelCross2 <- function(pop, elitepop, nCrosses, nFam, famSize, families, use, best = TRUE, nProgeny = 1, intWithin = 0.2, intAcross = 1, equalWeight = FALSE, useFamPrior = FALSE, verbose = FALSE){
	if(is.character(use)) use <- match.fun(use)

	if(any(table(elitepop@id) > 1)) msg(2, "WARNING: some elites repeated!")
	if(!all(pop@id %in% unlist(families))) stop("wrong family information... fix me!")
	famL <- lapply(families, function(x) x[x %in% pop@id])
	repFam <- sapply(famL, length) > 0
	famNum <- sum(repFam)
		
	if(famNum == 1) msg(2, "NOTE: Only one family represented!")
	# if(nInd(elitepop) > nFam) elitepop <- elitepop[sample(nInd(elitepop), nFam)]

	nFamSel <- ceiling(nFam * intAcross)
	if(famNum < nFamSel) msg(2, "NOTE: insufficient families represented to use specified intAcross. Using individuals from", famNum ,"represented families")
	msg(2, "Number of families available for crossing", min(famNum, nFamSel))

	famMeans <- sapply(famL[repFam], function(x) mean(use(pop[x])))
	famSel <- names(sort(famMeans, decreasing = TRUE))[1:min(famNum, nFamSel)]

	candFamL <- famL[famSel]
	inFamNum <- round(famSize * intWithin)

	parentL <- list()
	for(i in famSel) {
		parentL[[i]] <- selectInd(pop[candFamL[[i]]], nInd = min(nInd(pop[candFamL[[i]]]), inFamNum))
		parentL[[i]] <- parentL[[i]][order(use(parentL[[i]]), decreasing = TRUE)]
	}

	candL <- lapply(parentL, function(x) x@id)
	selFamMeans <- sapply(parentL, function(x) mean(use(x)))
	
	parents <- mergePopsRec(parentL)
	parents <- parents[order(use(parents), decreasing = TRUE)]
	parNames <- parents@id
	elNames <- elitepop@id
	
	needToAvoidInbreeding <- any(elNames %in% parNames)
	if(needToAvoidInbreeding) {
		msg(2, "Selection candidates already in elite pop! Avoiding within family crosses...")
	}

	if(nInd(parents) < nCrosses) msg(2, "NOTE: insufficient individuals to make desired number of unique crosses! Some parents will be reused...")

	if(!best){
		fm <- if (useFamPrior) famMeans[famSel] else selFamMeans 
		inFamWeights <- lapply(parentL, function(x) use(x)^2 / sum(use(x)^2))
		w <- if (equalWeight) rep(1/nFamSel, nFamSel) else fm^2 / sum(fm^2)
		whichFams <- sample(famSel, nCrosses, replace = TRUE, prob = w)
	} else {
		whichFams <- sapply(parNames, function(x) names(candL)[sapply(candL, function(xx) x %in% xx)])
	}

	selection <- list()
	for(i in 1:nCrosses){
		p1 <- if(best) parNames[i] else sample(candL[[whichFams[i]]], 1, prob = c(inFamWeights[[whichFams[i]]]))
		notRelated <- !elNames %in% famL[[whichFams[i]]]
		if(sum(notRelated) == 0) {
			nR <- !elitepop@id %in% famL[[whichFams[i]]]
			if(sum(nR) == 0){
				msg(1, "No unrelated selection candidates! Using lines from another family in same generation... ")
				otherFam <- sample(famSel[!famSel %in% whichFams[i]], 1)
				if(best) {
					p2 <- candL[[otherFam]][1]
				} else {
					p2 <- sample(candL[[otherFam]], 1, prob = c(inFamWeights[[otherFam]]))
				}
			} else {
				p2 <- sample(elitepop@id[nR], 1)
			}
		} else {		
			p2cand <- elNames[notRelated]
			elPar <- which(elNames == sample(p2cand, 1))
			p2 <- elNames[elPar]
			elNames <- elNames[-elPar]
		}
		selection[[i]] <- c(p1, p2)
		if(length(p1) > 1 | length(p2) > 1) break
	}

	selection <- do.call(rbind, selection)
	crossFamilies <- lapply(candL, function(x) x[x %in% selection[, 1]])
	if(verbose) msg(2, "Number of lines selected from each family for crosses:", sapply(crossFamilies, length))

	if(nProgeny > 1) selection <- selection[rep(1:nrow(selection), each = nProgeny), ] 
	newpop <- makeCross(mergePopsRec(list(pop, elitepop)), crossPlan = selection) 

	if(verbose) msg(2, "Selection mean:", round(mean(gv(newpop)), 6), "from pop mean:", round(mean(gv(pop)), 6))
	newpop
}


# pop <- selPop; nCrosses = nFam; use = "pheno"; intWithin = 0.5; intAcross = 1; equalWeight = FALSE; useFamPrior = FALSE
tradSelCross <- function(pop, nCrosses, nFam, famSize, families, use, nProgeny = 1, intWithin = 0.2, intAcross = 1, equalWeight = FALSE, useFamPrior = FALSE, verbose = FALSE){
# RCRS[[gen(j)]] <- tradSel(pop = selPop, nInd = min(selectRCRSi, nInd(selPop)), use = useIn,  trait = 1, simParam = simParam, nCrosses = nNuclear, nProgeny = nProgenyPerCrossIn)
	# lpop <- function(x, l, whichElem = NULL) {
	# 	whichElem <- if(is.null(whichElem)) which(sapply(l, function(elem) x %in% elem))
	# 	for(i in whichElem) l[[i]] <- l[[i]][!l[[i]] %in% x]
	# 	l
	# }
	if(is.character(use)) use <- match.fun(use)

	if(!all(pop@id %in% unlist(families))) stop("wrong family information... fix me!")
	famL <- lapply(families, function(x) x[x %in% pop@id])
	# famL[10] <- NULL
	repFam <- sapply(famL, length) > 0
	famNum <- sum(repFam)
	if(famNum == 1){
		msg(2, "NOTE: Only one family represented! Continuing by mating within familiy...")
		newpop <- selectCross(pop, nInd = round(nInd(pop) * intWithin), nCrosses = nCrosses)
	} else {
		nFamSel <- ceiling(nFam * intAcross)
		if(famNum < nFamSel) msg(2, "NOTE: insufficient families represented to use specified intAcross. Using individuals from all", famNum ,"represented families")
		msg(2, "Number of families used for crossing", min(famNum, nFamSel))

		famMeans <- sapply(famL[repFam], function(x) mean(use(pop[x])))
		famSel <- names(sort(famMeans, decreasing = TRUE))[1:min(famNum, nFamSel)]

		candFamL <- famL[famSel]
		inFamNum <- famSize * intWithin

		parentL <- list()
		for(i in famSel) {
			parentL[[i]] <- selectInd(pop[candFamL[[i]]], nInd = min(nInd(pop[candFamL[[i]]]), inFamNum))
			parentL[[i]] <- parentL[[i]][order(use(parentL[[i]]), decreasing = TRUE)]
		}

		selFamMeans <- sapply(parentL, function(x) mean(use(x)))
		
		parents <- mergePopsRec(parentL)
		parents <- parents[order(use(parents), decreasing = TRUE)]
		parNames <- parents@id

		if(choose(nInd(parents), 2) < nCrosses) {
			msg(2, "NOTE: insufficient individuals to make desired number of unique crosses! Some mate pairs will be repeated...")
			nTimes <- ceiling(nCrosses / choose(nInd(parents), 2))
			selL <- list(t(combn(parents@id, 2)))
			selection <- do.call(rbind, rep(selL, nTimes - 1))
			selection <- rbind(selection, selL[[1]][sample(1:nrow(selL[[1]])), ])[1:nCrosses, ]
		} else {
			selection <- list()
			if(intWithin == 1) {
				wp <- 1
				for(i in 1:nCrosses) {
					p1 <- parNames[1]
					diffFam <- !sapply(candFamL, function(f) p1 %in% f)
					cand <- parNames[parNames %in% unlist(candFamL[diffFam])]
					if(length(parNames) <= 1) {
						parNames <- parents@id
						wp <- wp + 1
					}
					p2 <- parNames[parNames %in% unlist(candFamL[diffFam])][[wp]]
					selection[[i]] <- c(p1, p2)
					parNames <- parNames[!parNames %in% c(p1, p2)]
				}
			} else {
				nFamPairs <- choose(nFamSel, 2)
				famPairs <- combn(famSel, 2)
				fm <- if (useFamPrior) famMeans[famSel] else selFamMeans 
				inFamWeights <- lapply(parentL, function(x) use(x)^2 / sum(use(x)^2))
				pbar <- combn(fm, 2, FUN = mean)
				w <- if (equalWeight) rep(1/nFamPairs, nFamPairs) else pbar^2 / sum(pbar^2)
	 			withRep <- if(nCrosses > nFamPairs) TRUE else FALSE
				whichPairs <- sample(1:ncol(famPairs), nCrosses, replace = withRep, prob = w)
				for(i in 1:nCrosses){
					pairi <- famPairs[, whichPairs[i]]
					p1 <- sample(parentL[[pairi[1]]]@id, 1, prob = c(inFamWeights[[pairi[1]]]))
					p2 <- sample(parentL[[pairi[2]]]@id, 1, prob = c(inFamWeights[[pairi[2]]]))
					selection[[i]] <- c(p1, p2)
				}
			}
			if(any(sapply(selection, function(x) length(unique(x))) != 2)) stop("something wrong happened!")
			selection <- do.call(rbind, selection)
		}

		if(nProgeny > 1) selection <- selection[rep(1:nrow(selection), each = nProgeny), ] 
		newpop <- makeCross(pop, crossPlan = selection) 
	}
	if(verbose) msg(2, "Selection mean:", round(mean(gv(newpop)), 6), "from pop mean:", round(mean(gv(pop)), 6))
	newpop
}

getAcc <- function(pop, simParam) {
	popTrue <- gv(pop)
	popPred <- ebv(pop)
	if(nrow(popTrue) == 1 | var(c(popTrue)) == 0 | var(c(popPred)) == 0) NA else cor(popTrue, popPred)
}

dummyFunc <- function(x, retrn) { retrn }

# get locus index
pullLoci <- function(simParam, snpChip = 1, asList = FALSE) {
	splitIndex <- function(lpc, loc) split(loc, rep(1:length(lpc), times = lpc))
	getIndex <- function(sites, nPerChr, loc){
		indexList <- splitIndex(sites, 1:sum(sites))
		locList <- splitIndex(nPerChr, loc)
		unlist(lapply(1:length(locList), function(i) indexList[[i]][locList[[i]]]))
	}
	sites <- SP$segSites
	QTLnPerChr <- simParam$traits[[snpChip]]@lociPerChr
	SNPnPerChr <- simParam$snpChips[[snpChip]]@lociPerChr
	QTLloc <- simParam$traits[[snpChip]]@lociLoc 
	SNPloc <- simParam$snpChips[[snpChip]]@lociLoc
	QTLsites <- getIndex(sites, QTLnPerChr, QTLloc)
	SNPsites <- getIndex(sites, SNPnPerChr, SNPloc)
	intersect(QTLsites, SNPsites)
	list(QTL = QTLsites, SNP = SNPsites)
}
# get all segSites, as pullSegSiteGeno doesnt function when there are diff segSites per chrom. Can this be true?
pullSegSites <- function(pop, returnMatrix = TRUE){
	rawToSum <- function(xk) {
		class(xk) <- "integer"
		rowSums(xk)
	}
	geno <- lapply(pop@geno, function(x) t(apply(x, 3, rawToSum)))
	if (returnMatrix) geno <- do.call(cbind, geno)
	geno
}

rlapply <- function(l, f = identity, level = 1, combine = list, counter = 1, ...){
	args <- list(...)
	if (counter < level){
		do.call(lapply, c(list(X = l, FUN = rlapply, f = f, level = level, combine = combine, counter = counter + 1), args))
	} else {
		result <- do.call(lapply, c(list(X = l, FUN = f), args))
		if (identical(combine, list)) return(result) else return(do.call(combine, result))
	}
}

estIntensity <- function(VDP, i, nT = nTrial, start = "trial1", end = "variety", estFunc = pheno, Gvar = estVg) {
	S <- mean(pheno(VDP[[end]][[gen(i - nT)]])) - mean(pheno(VDP[[start]][[gen(i - nT)]]))
	i <- S / sqrt(Gvar(VDP[[start]][[gen(i - nT)]]))
	if(is.matrix(i) & prod(dim(i)) == 1) i <- i[[1]] else if(is.matrix(i)) msg(2, "intensity has dimensions:", dim(i))
	i
}


dFSel <- function(dF, limit = 1, val = "selCrit", parentCols = c("p1", "p2"), returnPar = TRUE) {
	dFord <- list(fwd = order(dF[[val]]), rev = order(-dF[[val]]))
	dup <- duplicated(dF[[val]])
	dFord <- lapply(dFord, function(x) x[!dup])
	if(!(list({1:nrow(dF)}[!dup]) %in% dFord | list({nrow(dF):1}[!rev(dup)]) %in% dFord)) {
		stop("dF must be sorted in order to select!")
	}
	
	parMat <- as.matrix(dF[, parentCols])
	parents <- sort(unique(c(parMat)))
	parCount <- rep(0, length(parents))
	names(parCount) <- parents
	rows <- NULL

	for (i in 1:nrow(parMat)) {
		parenti <- parMat[i, ]
		if(all(parCount[parenti] < limit)) {
			parCount[parenti] <- parCount[parenti] + 1
			rows <- c(rows, i)
		}
	}
	if(returnPar) parMat[rows, , drop = FALSE] else rows
}


getSel <- function(selCrit, n, high = TRUE, variable = "selCrit", parentCols = c("p1", "p2"), maxP = NULL) {
	len <- if(is.data.frame(selCrit)) nrow(selCrit) else length(selCrit)
	if(len < n) n <- nrow(selCrit)
	if (is.data.frame(selCrit)){
		selCrit <- selCrit[order(selCrit[[variable]], decreasing = high), ]
		if(!is.null(maxP)){
			sel <- dFSel(selCrit, limit = maxP, val = variable, parentCols = parentCols)
			sel <- sel[1:min(nrow(sel), n), , drop = FALSE] 
		} else {
			sel <- as.matrix(selCrit[1:n, parentCols])
		}
	} else {
		sel <- names(sort(selCrit, decreasing = high))[1:n]
	}
	sel
}

getTrueQTLeff <- function(simParam, trait = 1) { simParam$traits[[trait]]@addEff }

getTrueBV <- function(pop, simParam, trait = 1) {
	M <- pullQtlGeno(pop)	
	u <- getTrueQTLeff(simParam)
	M %*% u
} 

getSelfVar <- function(M, u = NULL, fdiff = NULL) {
	if(is.null(u)) u <- rep(1, ncol(M))
	H <- M == 1
	Hu <- H %*% u^2 
	if(!is.null(fdiff)) Hu <- Hu * (1 - 2^(-fdiff))
	# if(ncol(Hu) == 1) Hu <- c(Hu)
	Hu
}

weightedQuantile <- function(mu, sigmasq, quant = 0.9, w = 0.5) {
	msg(2, "            weight parameter has value:", w)
	w * mu + (1-w) * qnorm(quant, sd = sqrt(sigmasq))
}

estVg <- function(pop, GSfit = NULL, GSfunc = RRBLUP) {
	if(is.null(GSfit)) GSfit <- GSfunc(pop)
	af <- getAF(pop)
	Vg <- GSfit@Vu
	Vg <- if(is.matrix(Vg) & prod(dim(Vg)) == 1) Vg[[1]]
	Vg * sum(2 * af * (1-af))
}

geth2 <- function(pop, GSfit) {
	Vg <- estVg(pop, GSfit)
	Vg / sum(Vg, GSfit@Ve)
}
# function to return expected quantiles from sampling DH individuals
# select individuals
simDHdist <- function(nSel, pop, GSfit, retQuant = FALSE, quant = 0.9, nDH = 200, weight = 0.5, nProgeny = 1, returnPop = TRUE, bigmem = FALSE, ...) {
	if(bigmem) {
		DH <- makeDH(pop, nDH = nDH)
		DH <- setEBV(DH, GSfit)
		DHebv <- ebv(DH)
		simDist <- split(DHebv, rep(1:nInd(pop), each =  nDH))
	} else {
		DHdist <- function(i){
			DH <- makeDH(pop[i], nDH = nDH)
			DH <- setEBV(DH, GSfit)
			ebv(DH)
		}
		simDist <- lapply(pop@id, DHdist)
	}
	expVar <- if(retQuant) sapply(simDist, quantile, probs = quant) else weightedQuantile(sapply(simDist, mean), sapply(simDist, var), quant, w = weight)
	names(expVar) <- pop@id
	if(returnPop){	
		selection <- getSel(expVar, n = nSel)
		if(nProgeny > 1) selection <- rep(selection, each = nProgeny)
		return(pop[selection])
	} else {
		return(expVar)
	}
}

# select pairs and cross
# pop <- RCRS[[lastRCRSgen]]; GSfit <- GSmodel[[lastGSmodel]]; use = ebv; nSel = selectRCRSi; nCrosses = nNuclear; maxCrossPerParent = 1; nDH = 200; quant = xInt; retQuant = FALSE;  w = 0.5; weight = 0.5; nSimCrosses = 1
simDHdistPairs <- function(nSel, pop, GSfit, nCrosses, use, retQuant = FALSE, quant = 0.9, nDH = 200, weight = 0.5, maxCrossPerParent = 0, nSimCrosses = 1, nProgeny = 1, verbose = FALSE, ...) {
	n <- nInd(pop)
	if (n < nSel) nSel <-  n
	nCombos <- choose(nSel, 2) 
	nEx <- if(nCombos < nCrosses) ceiling(nCrosses / nCombos) else 1 
	maxP <- if(maxCrossPerParent == 0 | nCombos <  nCrosses) nCrosses else maxCrossPerParent

	if(nSel < n) pop <- truncSel(pop, nSel = nSel, use = use)

	parents <- do.call(rbind, combn(pop@id, 2, simplify = FALSE))
	colnames(parents) <- c("p1", "p2")
	crosses <- rep(1:nrow(parents), each = nSimCrosses)
	popX <- makeCross(pop, parents[crosses, , drop = FALSE])

	if(verbose) msg(2, "simulating distribution of", nDH, "DH for", nSimCrosses, "crosses for each of", nCombos, "parental pairs\n")
	simVar <- do.call(simDHdist, getArgs(simDHdist, nSel = nInd(popX), pop = popX, GSfit = GSfit, retQuant = retQuant, quant = quant, nDH = nDH, weight = weight, returnPop = FALSE, ...))
	if(nSimCrosses > 1) simVar <- tapply(simVar, crosses, mean)

	selCrit <- data.frame(parents, selCrit = simVar)
	lenSel <- 0
	while(lenSel < nCrosses / nEx){
		if(lenSel > 0) {
			msg(2, "Not enough possible crosses with maxP =", maxP, "! Increasing maxP to", maxP + 1, "and retrying...\n")
			maxP <- maxP + 1
		}
		selection <- getSel(selCrit, n = nCrosses, high = TRUE, maxP = maxP)
		lenSel <- nrow(selection)
	}
	if(verbose) print(table(selection))

	if(nEx > 1) selection <- selection[rep(1:nrow(selection), times = nEx)[1:nCrosses], ]
	if(nProgeny > 1) selection <- selection[rep(1:nrow(selection), each = nProgeny), ] 
	makeCross(pop, crossPlan = selection) 
}

# select individuals
expDist <- function(nSel, pop, GSfit, use, quant, returnQuant = TRUE, pullGeno = pullSnpGeno, updateEBV = FALSE, weight = 0.5, nProgeny = 1, ...){
	expVar <- do.call(getSelfVar, getArgs(getSelfVar, M = pullGeno(pop), u = GSfit@markerEff, ...))
	if(returnQuant) {
		if(updateEBV) pop <- setEBV(pop, GSfit)
		parVal <- ebv(pop)
		expVar <- weightedQuantile(mu = parVal, sigmasq = expVar, quant = quant, w = weight)
	} 
	if(ncol(expVar) == 1) expVar <- expVar[, 1]
	selection <- getSel(expVar, n = nSel, high = TRUE)
	if(nProgeny > 1) selection <- rep(selection, each = nProgeny)
	pop[selection]
}

maxVar <- function(pop, GSfit, nSel, nCrosses, use, weightLoci = FALSE, pullGeno = pullSnpGeno, maxCrossPerParent = 0, verbose = FALSE, nProgeny = 1, ...){
	n <- nInd(pop)
	if (n < nSel) nSel <-  n
	nCombos <- choose(nSel, 2)
	nEx <- if(nCombos < nCrosses) ceiling(nCrosses / nCombos) else 1 
	maxP <- if(maxCrossPerParent == 0 | nCombos <  nCrosses) nCrosses else maxCrossPerParent
	
	if(nSel < n) pop <- truncSel(pop, nSel = nSel, use = use)

	M <- pullGeno(pop)
	K <- if(weightLoci) genCov(M, u = c(GSfit@markerEff)) else genCov(M)
	covL <- data.frame(which(lower.tri(K), arr.ind = TRUE), selCrit = K[lower.tri(K)])
	selCrit <- data.frame(covL, p1 = colnames(K)[covL$col], p2 = colnames(K)[covL$row]) 
	
	lenSel <- 0
	while(lenSel < nCrosses / nEx){
		if(lenSel > 0) {
			msg(2, "Not enough possible crosses with maxP =", maxP, "! Increasing maxP to", maxP + 1, "and retrying...\n")
			maxP <- maxP + 1
		}
		selection <- getSel(selCrit, n = nCrosses, high = FALSE, maxP = maxP)
		lenSel <- nrow(selection)
	}
	if(verbose) print(table(selection))
	
	if(nEx > 1) selection <- selection[rep(1:nrow(selection), times = nEx)[1:nCrosses], ]
	if(nProgeny > 1) selection <- selection[rep(1:nrow(selection), each = nProgeny), ] 
	makeCross(pop, crossPlan = selection) 
}

solqp <- function(pop, GSfit, use, nCrosses, simParam, lambda = NULL, fthresh = NULL, gain = NULL, truncqp = NULL, allowSelf = FALSE, weightLoci = FALSE, pullGeno = pullSnpGeno, verbose = FALSE, nProgeny = 1, maxProp = 1, ...){
	suppressMessages(require(LowRankQP))
	inbreedingCoef <- function(cee, Kmat) 1/2 * crossprod(cee, Kmat) %*% cee
	expectedGain <- function(cee, gebvs) crossprod(cee, gebvs )
	betterSample <- function(x, ...) x[sample(length(x), ...)]
	if(is.character(use)) use <- match.fun(use)

	n <- nInd(pop)
	
	M <- pullGeno(pop, simParam = simParam)
	K <- if(weightLoci) genCov(M, u = c(GSfit@markerEff)) else genCov(M)	
	gebvs <- use(pop)
	popTruth <- genParam(pop, simParam = simParam)

	Vg <- popTruth$varG
	if(var(gebvs) == 0){
		msg(2, "No variance in GEBVs! random intermating progressing...\n")
		selection <- randomCross(pop, nFam = nCrosses, nProgeny = nProgeny)
	} else {
		if(is.null(lambda) & is.null(fthresh) & is.null(gain) & is.null(truncqp)) stop("lambda {0, 1}, fthresh {>0}, gain {>0} or truncqp {>0} must be provided a value")
		if(!is.null(gain)) lambda <- 1 else if(is.null(lambda)) lambda <- 0:100 * 0.01 else if (lambda > 1 | lambda < 0) stop("Supplied lambda values must be between 0 and 1.")

		f <- NULL
		g <- NULL
		A <- matrix(1, nrow = 1, ncol = n)
		u <- matrix(1, ncol = 1, nrow = n) * maxProp
		b <- 1
		log <- list()
		cee <- list()
		for(k in 1:length(lambda)){
			H <- 2 * lambda[k] * K # the one half is included in the optimization. So where does the 2 come from? Should be 1?
			d <- if(is.null(gain)) -(1 - lambda[k]) * gebvs else rep(0, length(gebvs))
			if(!is.null(gain)) {
				b <- c(1, gain)
				A <- rbind(A, c(gebvs))
			}
			log[[k]] <- capture.output({solutionqp <- suppressMessages(invisible(LowRankQP(Vmat = H, dvec = d, Amat = A, bvec = b, uvec = u, method = "LU", verbose = FALSE)))})
			cee[[k]] <- solutionqp$alpha
			f <- c(f, c(inbreedingCoef(solutionqp$alpha, K)))
			g <- c(g, c(expectedGain(solutionqp$alpha, gebvs)))
		}
		if (!is.null(fthresh) & is.null(gain)) {
			whichLambda <- if(fthresh < min(f)) which.min(f) else which(f == max(f[f <= fthresh]))
		} else if (!is.null(truncqp)) {
			zero <- 1 / (2 * nCrosses)
			nPars <- sapply(cee, function(x) sum(x >= (zero)))
			whichLambda <- if(truncqp > max(nPars)) which.max(nPars) else which(nPars == min(nPars[nPars >= truncqp]))
			if(length(whichLambda) > 1) whichLambda <- whichLambda[1]
		} else {
			whichLambda <- 1
		}
		msg(2, "lambda:", lambda[whichLambda])
		propPar <- round(cee[[whichLambda]] * 2 * nCrosses)
		rownames(propPar) <- pop@id
		pars <- rep(pop@id, times = propPar)
		if (length(unique(pars)) == 1){
			selection <- do.call(rbind, rep(list(rep(unique(pars), 2)), nCrosses))
		} else {
			parList <- list()
			index <- 1:length(pars)
			for(i in 1:min(nCrosses, floor(length(pars) / 2))) {
				p1 <- betterSample(index, 1)
				samplep2 <- index[pars[index] != pars[p1]]
				p2 <- if(allowSelf) betterSample(index[-p1], 1) else betterSample(samplep2, 1)
				if(length(p2) == 0) next
				if(pars[p1] == pars[p2]) {msg(2, "oops\n"); break}
				crossi <- c(p1, p2)
				parList[[i]] <- pars[crossi]
				index <- index[!index %in% crossi]
			}
			selection <- do.call(rbind, parList)
		}
	}
	if(nProgeny > 1) selection <- selection[rep(1:nrow(selection), each = nProgeny), ] 
	makeCross(pop, crossPlan = selection, simParam = simParam) 
}

solqpOut <- function(pop, GSfit, use, nSel, nProgeny, nGenOut, nGenThisYr, simParam, limitN = 0, lambdaOut = NULL, fthreshOut = NULL, gainOut = NULL, truncqpOut = NULL, verbose = FALSE, ...){
	params <- list(lambdaOut = lambdaOut, fthreshOut = fthreshOut, gainOut = gainOut, truncqpOut = truncqpOut)
	for(i in names(params)){
		if (length(params[[i]]) > 1 & !is.list(params[[i]])) stop(paste0(i, " must be length 1 or a list."))
		if(length(params[[i]]) <= 1) params[[i]] <- rep(list(params[[i]]), nGenOut) else if(length(params[[i]]) != nGenOut) stop(paste0(i, " is the wrong length, must be 1 or cyclePerYr - pullCycle"))
	}
	pop <- setEBV(pop, GSfit)
	if(nGenOut == 0){
		pop <- selectInd(pop, nInd = nSel, use = use)
	} else {
		N <- if(limitN > 1) rep(nSel, nGenOut) else if(limitN > 0) c(rep(nInd(pop), nGenOut - 1), nSel) else rep(nInd(pop), nGenOut)
		gs <- 1
		i <- 0
		if(verbose) msg(2, "solqpOut Initial Pop Mean:", round(mean(pop@gv), 6), "PopVar", round(varA(pop), 6))

		while(i < nGenOut){
			i <- i + 1
			pop <- do.call(solqp, getArgs(solqp, pop = pop, GSfit = GSfit, use = use, nCrosses = N[i], simParam = simParam,
						   nProgeny = nProgeny, lambda = params$lambdaOut[[i]], fthresh = params$fthreshOut[[i]], gain = params$gainOut[[i]], truncqp = params$truncqpOut[[i]]))#, ...))
			pop <- setEBV(pop, GSfit)
			if(verbose) msg(2, "solqpOut Pop Mean:", round(mean(pop@gv), 6), "PopVar", round(varA(pop), 6))
			if(verbose) msg(2, "solqpOut pop size:", nInd(pop))
		}
	}
	if(limitN < 1) pop <- selectInd(pop, nInd = nSel, use = use)
	if(verbose) msg(2, "solqpOut Final Pop Mean:", round(mean(pop@gv), 6), "PopVar", round(varA(pop), 6))
	pop
}

withinFamSel <- function(pop, GSfit, use, int, nProgeny, nGenInbr, nGenThisYr, lambdaInbr = NULL, fthreshInbr = NULL, gainInbr = NULL, truncqpInbr = NULL, ssd, simParam, ...){
	nFam <- nInd(pop)
	nProgPerFam <- nProgeny / int
	if (nProgPerFam %% 1 != 0) {
		nProgPerFam <- round(nProgPerFam)
		nSelToTrial <- round(nProgPerFam * withinFamInt)
		msg(2, "\nNOTE: Selection intensities within familiy have been rounded to the nearest integer resulting in", nSelToTrial, "progeny per family selected from", nProgPerFam, "progeny per family\n")
	}
	if (ssd) {
		i <- 0
		while(i <= nGenInbr) {
			nP <- if(i == 0) nProgPerFam else 1 
			pop <- self(pop, nProgeny = nP)
			i <- i + 1
		}
	} else {
		pop <- makeDH(pop, nDH = nProgPerFam)
	}

	pop <- setEBV(pop, GSfit, simParam = simParam)
	fams <- split(1:(nProgeny * nFam), rep(1:nFam, each = nProgeny))
	faml <- list()
	for (j in names(fams)){
		faml[[j]] <- selectInd(pop[fams[[j]]], nInd = nProgeny, trait = 1, use = use) 
	}
	mergePops(faml)
}

expDistPairs <- function(pop, GSfit, nSel, quant, nCrosses, use, returnQuant = TRUE, weightLoci = FALSE, pullGeno = pullSnpGeno, maxCrossPerParent = 0, Gvar = estVg, weight = 0.5, nProgeny = 1, verbose = FALSE, ...) {
	n <- nInd(pop)
	if (n < nSel) nSel <-  n
	nCombos <- choose(nSel, 2)
	nEx <- if(nCombos < nCrosses) ceiling(nCrosses / nCombos) else 1 
	maxP <- if(maxCrossPerParent == 0 | nCombos <  nCrosses) nCrosses else maxCrossPerParent
	
	if(nSel < n) pop <- truncSel(pop, nSel = nSel, use = use)
	parVal <- ebv(pop)
	rownames(parVal) <- pop@id

	M <- pullGeno(pop)
	K <- if(weightLoci) genCov(M, u = c(GSfit@markerEff)) else genCov(M)

	Vg <- do.call(Gvar, getArgs(Gvar, pop = pop, GSfit = GSfit, ...))
	if (prod(dim(Vg)) > 1) stop("can only handle a single trait!") else Vg <- Vg[[1]]
	parents <- do.call(rbind, combn(pop@id, 2, simplify = FALSE))
	colnames(parents) <- c("p1", "p2")
	pE <- getPopMeanVar(parVal, K, Vg)
	Eq <- weightedQuantile(mu = pE$pbar, sigmasq = pE$crossvar, quant = quant, w = weight)
	# Eq <- w * pE$pbar + (1-w) * qnorm(quant, sd = sqrt(pE$crossvar)) # does this make sense?
	selCrit <- data.frame(parents, selCrit = Eq)

	lenSel <- 0
	while(lenSel < nCrosses / nEx){
		if(lenSel > 0) {
			msg(2, "Not enough possible crosses with maxP =", maxP, "! Increasing maxP to", maxP + 1, "and retrying...\n")
			maxP <- maxP + 1
		}
		selection <- getSel(selCrit, n = nCrosses, high = TRUE, maxP = maxP) # rerun with high = TRUE??
		lenSel <- nrow(selection)
	}
	if(verbose) print(table(selection))
	
	selection <- getSel(selCrit, n = nCrosses)
	if(nEx > 1) selection <- selection[rep(1:nrow(selection), times = nEx)[1:nCrosses], ]
	if(nProgeny > 1) selection <- selection[rep(1:nrow(selection), each = nProgeny), ] 
	makeCross(pop, crossPlan = selection)
}

pherFuncID <- function(x, p) {x}
pherFuncMax <- function(x, p) {(x / max(x))^p}
quant <- function(x, sigma, w) {w * mean(x) + (1 - w) * sigma * sd(x)}
popQuant <- function(sel, ebvs, sigma = 2, w = 0.5) {quant(x = ebvs[sel], sigma = sigma, w = 0.5)}
initFunc <- function(x) {rep(1, length(x))} # can change later 
gsQuant <- function(sel, ebvs, sigma = 2, w = 0.5, nCrosses = 100, nProgenyPerCross = 1) {
	if(!all(ebv(pop[sel]) == ebvs[sel])) stop("ebvs of pop dont match those of ant!")
	simCrosses <- randomCross(pop[sel], nFam = nCrosses, nProgeny = nProgenyPerCross)
	simPop <- makeCross(pop[sel], simCrosses)
	simPop <- setEBV(simPop, GSfit)
	quant(ebv(simPop), sigma = sigma, w = w)
}

acOpt <- function(x, n, targetFunc, pherFunc, xAt0 = FALSE, evapRate = 0.05, nAnts = 500, maxIter = 300, countThresh = 1000, showResult = FALSE, pherPower = 1, dumbAnt = FALSE, plateau = 100, returnStats = FALSE){
	evap <- function(oldP, newP, evapRate) {(1 - evapRate) * oldP + newP }

	if(is.matrix(x)) {if(ncol(x) == 1) x <- c(x) else stop("I cant handle multiple traits yet! dim = ", dim(x))}
	
	if(xAt0){
		minx <- min(x)
		x <- x - minx
	}
	N <- length(x)
	pheromone <- initFunc(x)

	bestAnt <- NULL
	noP = rep(0, N)
	bestestPath = 0
	lastBestPath <- bestestPath
	pathCounter = 0 
	iter = 0
	if(showResult) {
		plotBest <- NULL 
		plotMean <- NULL 
	}

	while(pathCounter <= countThresh & iter < maxIter){
		iter = iter + 1
		ant <- list()
		path <- list()
		antPheromone <- list()
		for(i in 1:nAnts){
			ant[[i]] <- sample(1:N, n, prob = pheromone / sum(pheromone))
			path[[i]] <- targetFunc(sel = ant[[i]], ebvs = x)
			# path[[i]] <- quant(x[ant[[i]]], sigma = 2, w = 0.5)
			Pi <- noP
			Pi[ant[[i]]] <- path[[i]]
			antPheromone[[i]] <- Pi
		}
		bestPath <- which.max(path)
		if (path[[bestPath]] == lastBestPath) {
			pathCounter <- pathCounter + 1
		} else {
			lastBestPath <- path[[bestPath]]
			pathCounter <- 0
		}
		if(path[[bestPath]] >= bestestPath) {
			bestestPath <- path[[bestPath]] 
			bestAnt <- ant[[bestPath]]
		}
		newPheromone <- pherFunc(Reduce("+", antPheromone), p = pherPower)	
		pheromone <- evap(pheromone, newPheromone, evapRate = evapRate)
	}
	return(bestAnt)
}

# pop = RCRS[[lastRCRSgen]]; GSfit = GSmodel[[lastGSmodel]]; acTrunc = 1; evapRate = 0.05; nAnts = 500; pherPower = 1.5; nSel = selectRCRSi; nCrosses = nNuclear; use = ebv; pullGeno = pullSnpGeno; weightLoci = FALSE; maxCrossPerParent = 1; nProgeny = 1
ACquant <- function(pop, GSfit, nSel, nCrosses, use, acTrunc = 1, evapRate = 0.05, nAnts = 500, pherPower = 1.5, verbose = FALSE, nProgeny = 1, ...){
	if(is.character(use)) use <- match.fun(use)
	n <- nInd(pop)
	if (n < nSel) nSel <-  n
	nCombos <- choose(nSel, 2)
	nEx <- if(nCombos < nCrosses) ceiling(nCrosses / nCombos) else 1 
	# maxP <- if(maxCrossPerParent == 0 | nCombos <  nCrosses) nCrosses else maxCrossPerParent
	popVar <- varA(pop)
	popMean <- mean(gv(pop))

	if (acTrunc < 1) pop <- truncSel(pop, nSel = n * acTrunc, use = use)
	if (n == nSel) {
		selectedParents <- 1:n
	} else {
		selectedParents <- acOpt(use(pop), n = nSel, xAt0 = TRUE, targetFunc = popQuant, pherFunc = pherFuncMax, evapRate = evapRate, nAnts = nAnts, pherPower = pherPower)
	}
	selection <- randomCross(pop[selectedParents], nFam = nCrosses, nProgeny = nProgeny)

	if(nEx > 1) selection <- selection[rep(1:nrow(selection), times = nEx)[1:nCrosses], ]
	if(nProgeny > 1) selection <- selection[rep(1:nrow(selection), each = nProgeny), ] 
	newpop <- makeCross(pop, crossPlan = selection) 
	if(verbose){
		msg(2, "Selected population variance diff:", {varA(newpop) - popVar})
		msg(2, "Selected population mean diff:", {mean(gv(newpop)) - popMean})
	}
	newpop
}


randomCross <- function(pop, nFam, nProgeny = 1){ # note this is just a random sampler, to illustrate how one might build a function to pick pairs. 
	allCrosses <- combn(pop@id, 2)
	resample <- if (nFam > ncol(allCrosses)) TRUE else FALSE 
	crosses <- allCrosses[, sample(1:ncol(allCrosses), nFam, replace = resample)]
	if (nProgeny > 1) crosses[, rep(1:nFam, each = nProgeny)]
	t(crosses)
}

selectInd2 <- function(pop, nSel, use, trait = 1, selFunc = identity){
	if(is.character(use)) use <- match.fun(use)
	sel <- use(pop)
	if(ncol(sel) < 1) stop("Something is wrong! I dont appear to have a phenotype/GEBV! perhaps I was never predicted?")
	names(sel) <- pop@id
	selection <- getSel(selFunc(sel), n = nSel)
	pop[selection]
}

truncSel <- function(pop, nSel, use, traits = 1, ...) { do.call(selectInd2, getArgs(selectInd2, pop = pop, nSel = nSel, use = use, trait = traits, ...)) }

truncCross <- function(pop, nSel, nCrosses, use, nProgeny = 1, crossFunc = randomCross, traits = 1, ...) {
	if(nSel < nInd(pop)) selPop <- do.call(selectInd2, getArgs(selectInd2, pop = pop, nSel = nSel, use = use, trait = traits, ...)) else selPop <- pop
	selection <- do.call(crossFunc, getArgs(crossFunc, pop = selPop, nFam = nCrosses, nProgeny = nProgeny, ...))
	makeCross(pop, selection)
}

checkFit <- function(pop){
	GSfit <- GSfunc(pop, traits = 1, use = "pheno", snpChip = 1, simParam = simParam)
	pop <- setEBV(pop, GSfit, simParam = simParam)
	M <- pullSnpGeno(pop)

	msg(2, "alphaSimR RRBLUP iterations:", GSfit@iter)
	
	msg(2, "alphaSimR RR-BLUP Vu:",  GSfit@Vu, ", Ve:", GSfit@Ve)

	require(EMMREML)
	rr <- emmreml(y = pheno(pop), X = matrix(1, nInd(pop), 1), Z = pullSnpGeno(pop), K = diag(simParam$snpChips[[1]]@nLoci))
	af <- getAF(pop)
	msg(2, "emreml RR-BLUP Vu:",  rr$Vu, ", Ve:", rr$Ve)
	msg(2, "Correlation of marker effect estimates:", cor(rr$u, GSfit@markerEff))
	msg(2, "emreml RR-BLUP Vg:", rr$Vu * sum(2 * af * (1-af)))
	M %*% rr$u - ebv(pop) 

	K <- vanRaden1()
	mean(diag(K))
	gblup <- emmreml(y = pheno(pop), X = matrix(1, nInd(pop), 1), Z = diag(nInd(pop)), K = K)
	msg(2, "emreml GBLUP Vg:",  gblup$Vu, ", Ve:", gblup$Ve)

	msg(2, "Prediction accuracy alphaSimR:",  getAcc(pop))
	msg(2, "Prediction accuracy emmreml:",  cor(gblup$u, bv(pop)))
	msg(2, "Correlation of emmreml and alphaSimR genetic effect estimates:",  cor(gblup$u, ebv(pop)))
}

getPopMeanVar <- function(parVal, parCov, Vg){
	if (ncol(parVal) > 1) stop("cannot use more than 1 trait...") 
	pbar <- combn(parVal, 2, mean)
	pCovar <- parCov[lower.tri(parCov)] * Vg
	pVarSum <- combn(diag(parCov) * Vg, 2, sum) 
	crossvar <- pVarSum - 2 * pCovar
	crossvar[crossvar < 0] <- 0
	list(pbar = pbar, crossvar = crossvar, pcovar = pCovar)
}

vanRaden1 <- function(M){
	Z <- scale(M, scale = FALSE)
	p <- attributes(Z)[["scaled:center"]] / 2
	ZZt <- tcrossprod(Z)
	ZZt / (2 * crossprod(p, 1-p)[[1]])
}

wrightsF <- function(M, returnNA = TRUE){
	n <- nrow(M)
	p <- colMeans(M) / 2
	d <- colSums(M == 1) / n # assumes M in 0, 1, 2
	v <- 2 * p *(1-p)
	wF <- 1 - d/v
	wF[is.na(wF)] <- if(returnNA) NA else 0
}


genCov <- function(M, u = NULL, absU = TRUE, sumVar = TRUE, scaleD = TRUE, inclm = TRUE, p = NULL){
	if(is.matrix(u)) u <- c(u)
	Z <- scale(M, scale = FALSE)
	m <- if(inclm) ncol(Z) else 1
	if(is.null(p)) p <- attributes(Z)[["scaled:center"]] / 2
	v <- 2 * p * (1 - p)
	if(is.null(u) & sumVar){
		ZDZt <- tcrossprod(Z) / sum(v)
	} else {
		if(all(u == 1) & length(u) == 1) u <- rep(1, ncol(Z))
		d <- u
		if(absU) d <- abs(d)
		if(!sumVar) {
			seg <- v != 0
			d <- d[seg] / (v*m)[seg]
			Z <- Z[, seg]
		}
		if(scaleD) d <- d / mean(d)
		ZDZt <- tcrossprod(Z %*% diag(d), Z)
		if(sumVar) ZDZt <- ZDZt / sum(v)
	}
	ZDZt
}

genCov2 <- function(M, u = NULL, absD = TRUE, sumVar = TRUE, scaleD = TRUE, inclm = TRUE, p = NULL){
	if(is.matrix(u)) u <- c(u)
	Z <- scale(M, scale = FALSE)
	m <- if(inclm) ncol(Z) else 1
	if(is.null(p)) p <- attributes(Z)[["scaled:center"]] / 2
	v <- 2 * p * (1 - p)
	if(is.null(u) & sumVar){
		ZDZt <- tcrossprod(Z) / sum(v)
	} else {
		if(all(u == 1) & length(u) == 1) u <- rep(1, ncol(Z))
		d <- if (!is.null(p)) 1 / v else u
		if(absD) d <- abs(d)
		if(!sumVar) {
			seg <- v != 0
			d <- d[seg] / (v*m)[seg]
			Z <- Z[, seg]
		}
		if(scaleD) d <- d / mean(d)
		ZDZt <- tcrossprod(Z %*% diag(d), Z)
		if(sumVar) ZDZt <- ZDZt / sum(v)
	}
	ZDZt
}


getTotalIntensity <- function(x) {
    S <- x$vy - x$gv[x$RCRSyr - 1]
	i <- S / x$Vg[x$RCRSyr - 1]
	list(S = S, i = i)
}


getPopStats <- function(resultL, meanVariety = TRUE, verbose = FALSE){
    VDPparam <- rlapply(resultL[["VDP"]], f = genParam, level = 2)
    RCRSparam <- lapply(resultL[["RCRS"]], genParam)

    nYr <- length(resultL[["VDP"]][[1]])
    GScylcePerYr <- (length(resultL[["RCRS"]]) - 1) / nYr 
    yr <- 1:nYr
    Ryr <- yr * GScylcePerYr
    Rcyc <- c(0, 1:(GScylcePerYr * nYr))

    VgRCRS <- sapply(RCRSparam, "[[", "varG")
    gvRCRS <- sapply(RCRSparam, function(x) mean(x$gv_a) + x$gv_mu) # this is correct
    sRCRS <- gvRCRS[-1] - gvRCRS[-length(gvRCRS)]
    iRCRS <- sRCRS / sqrt(VgRCRS[-length(VgRCRS)])

    VgVDP <- rlapply(VDPparam, "[[", i = "varG", level = 2, combine = c)
    gvVDP <- rlapply(VDPparam, function(x) mean(x$gv_a) + x$gv_mu, level = 2, combine = c)
    sVDP <- gvVDP$variety - gvVDP$trial1
    iVDP <- sVDP / sqrt(VgVDP$trial1)
  
  	RyrIndex <- Ryr - GScylcePerYr + 1
  	sTotal <- gvVDP$variety - gvRCRS[RyrIndex]
	iTotal <- sTotal / sqrt(VgRCRS[RyrIndex])

    gvVariety <- lapply(VDPparam[["variety"]], function(x) x$gv_a + x$gv_mu)
    SDgRCRS <- sqrt(VgRCRS)

    varMean <- gvVDP[["variety"]]
	nVariety <- sapply(gvVariety, nrow)
	Yvariety <- unlist(gvVariety)
	Xvariety <- rep(Ryr[1:length(nVariety)], times = nVariety)


    RCRSacc <- resultL$predAcc[["RCRS"]]
    VDPacc <- resultL$predAcc[["VDP"]]
    RCRSoutAcc <- resultL$predAcc[["RCRSout"]]
    VDPinAcc <- resultL$predAcc[["VDPin"]]

    theorMax <- maxBv(resultL$SP)
	
	nVar = unique(nVariety)

    return(list(SP = resultL$SP, paramL = resultL$paramL, Rcyc = Rcyc, varMean = varMean, sdRCRS = SDgRCRS, 
				VgRCRS = VgRCRS, VgVDP = VgVDP, gvRCRS = gvRCRS, gvVDP = gvVDP,
    			sRCRS = sRCRS, iRCRS = iRCRS, sVDP = sVDP, iVDP = iVDP, sTotal = sTotal, iTotal = iTotal, 
    			nVar = nVar, vx = Xvariety, vy = Yvariety, RCRSyr = Ryr, RCRSacc = RCRSacc, VDPacc = VDPacc,
    			RCRSoutAcc = RCRSoutAcc, VDPinAcc = VDPinAcc, theorMax = theorMax))
}


getYrange <- function(simR) { range(c(simR$gv + simR$sdRCRS, simR$gv - simR$sdRCRS, simR$vy)) }

invertList <-  function(ll) {
    nms <- unique(unlist(lapply(ll, function(X) names(X))))
    ll <- lapply(ll, function(X) setNames(X[nms], nms))
    ll <- apply(do.call(rbind, ll), 2, as.list)
    lapply(ll, function(X) X[!sapply(X, is.null)])
}

plotPop <- function(simL, Rgen = RCRSgen, vLine = "none", popcol = "#000000", alpha = "0D", alphaMean = "0D", pch = 1) {
	polycol <- paste0(popcol, alphaMean)
	popcol <- paste0(popcol, alpha)
	xpoly <- c(Rgen, rev(Rgen), Rgen[1])
	ypoly <- c(simL$gvRCRS + simL$sdRCRS, rev(simL$gvRCRS - simL$sdRCRS), simL$gvRCRS[1] + simL$sdRCRS[1])

	polygon(x = xpoly, y = ypoly, col = polycol, border = NA)
	lines(x = Rgen, y = simL$gvRCRS, type = "l", col = popcol, lwd = 2)
	points(simL$vx, simL$vy, col = popcol, pch = pch)
    if (vLine == "linear") {
		abline(with(simL, lm(vy ~ vx)), col = popcol, lwd = 2)
    } else if (vLine %in% c("poly", "curve")){	    	
    	fit <- if(vLine == "poly") with(simL, lm(vy ~ poly(vx, 2))) else with(simL, loess(vy ~ vx))
		smx <- seq(min(simL$vx), max(simL$vx), by = 0.1)
		lines(smx, predict(fit, newdata = data.frame(vx = smx)), type = "l", col = popcol, lwd = 1)
    } else if (vLine == "connect"){
    	lines(simL$vx, simL$vy, col = popcol, lwd = 2)
    }
}


simPlot <- function(popList, cols = "#000000", popLabs = NULL, varLine = "none", meanVariety = TRUE, legendPos = "topleft", plotReps = FALSE, plotVg = TRUE, plotSelInt = TRUE){

    if (length(cols) != length(popList)) stop("cols must be same length as popList!")
    lineCol <- paste0(cols, "FF")
    ptCol <- paste0(cols, "FF")
    polyCol <- paste0(cols, "4D")

    if (is.null(popLabs)) popLabs <- names(popList)

  	nVar <- unique(unlist(rlapply(popList, "[[", i = "nVar", level = 2, combine = c)))
    cyclePerYr <- unique(unlist(rlapply(popList, function(x) x[["paramL"]][["cyclePerYr"]], level = 2, combine = c)))
    if(length(nVar) > 1) warning("number of varieties differ between pops!")
    if(length(nVar) > 1) warning("cycles per year differ between pops!")

    simStats <- lapply(popList, function(x) lapply(x, function(xx) xx[!names(xx) %in% c("SP", "paramL", "VgVDP", "gvVDP")]))

    simStatsInv <- lapply(simStats, invertList)
    simReps <- rlapply(simStatsInv, level = 3, combine = rbind) 
    simAvg <- rlapply(simReps, f = colMeans, level = 2, na.rm = TRUE)

    RCRSyr <- simAvg[[1]][["RCRSyr"]]
    RCRSgen <- simAvg[[1]][["Rcyc"]]
    yr <- RCRSyr / cyclePerYr
    xlims <- range(c(0, RCRSgen))
    ylims <- range(sapply(simReps, getYrange)) * 1.1

    if (meanVariety) {
    	simAvg <- lapply(simAvg, function(x) {x[["vy"]] <- x[["varMean"]]; x[["vx"]] <- x[["RCRSyr"]]; x})
    }
    plot(NA, xlim = xlims, ylim = ylims, xaxt = "n", xlab = "generation", ylab = "standardized genetic value")
    axis(1, at = c(0, RCRSyr), labels = c(0, yr))

    for (i in 1:length(popList)){
    	if(plotReps) invisible(lapply(simStats[[i]], plotPop, Rgen = RCRSgen, popcol = cols[i]))
    	plotPop(simAvg[[i]], popcol = cols[i], alpha = "FF", alphaMean = "4D", Rgen = RCRSgen, vLine = varLine, pch = 16)
    }

    if (length(popList) > 1) {
        legend(legendPos, legend = popLabs, col=cols, lty = 1, lwd = 2, pch = 16)
      } else {
        legend(legendPos, legend = c("RCRS mean", expression(paste('RCRS ', sigma[g])), "Variety mean"), 
         bty = "n",
         col = c(lineCol, "black", ptCol),
         lty = c(1, 0, 0), lwd = c(2, 0, 0),
         pch = c(NA, 22, 16),
         pt.bg = c(NA, polyCol, ptCol),
        )
	}
	if(plotVg){
		ylims2 <- c(0, max(sapply(simAvg, "[[", "VgRCRS")))
		plot(NA, xlim = xlims, ylim = ylims2, xaxt = "n", xlab = "generation", ylab = "Vg", main = "Genetic variance across generations")
	    axis(1, at = c(0, RCRSyr), labels = c(0, yr))
	    for (i in 1:length(simAvg)) {
	    	lines(simAvg[[i]][["Rcyc"]], simAvg[[i]][["VgRCRS"]], type = "l", lwd = 2, lty = 1, col = cols[[i]])
	    }	
	    legend("topright", legend = popLabs, col=cols, lty = 1, lwd = 2, pch = 16)
	}



# THIS NEEDS TO BE CLEANED UP!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	if(plotSelInt){
	
		sVDP <- lapply(simAvg, "[[", "sVDP")
		iVDP <- lapply(simAvg, "[[", "iVDP")
		sRCRS <- lapply(simAvg, "[[", "sRCRS")
		iRCRS <- lapply(simAvg, "[[", "iRCRS")

		ylims3 <- range(unlist(sVDP))
		# ylims4 <- range(unlist(int))
		ylims4 <- c(-2, max(unlist(iVDP)))

		ylims5 <- range(unlist(sRCRS))
		# ylims4 <- range(unlist(int))
		ylims6 <- c(-2, max(unlist(iRCRS)))


		plot(NA, xlim = xlims, ylim = ylims3, xaxt = "n", xlab = "generation", ylab = "S", main = "Selection differential across generations in VDP")
	    axis(1, at = c(0, RCRSyr), labels = c(0, yr))
	    
	    for (i in 1:length(sRCRS)) {
	    	lines(RCRSyr, sVDP[[i]], type = "l", lwd = 2, lty = 1, col = cols[[i]])
	    }	
	    legend("topright", legend = popLabs, col=cols, lty = 1, lwd = 2, pch = 16)

	   	plot(NA, xlim = xlims, ylim = ylims4, xaxt = "n", xlab = "generation", ylab = "S/Vg", main = "Selection intensity across generations in VDP")
	    axis(1, at = c(0, RCRSyr), labels = c(0, yr))
	    
	    for (i in 1:length(iVDP)) {
	    	lines(RCRSyr, iVDP[[i]], type = "l", lwd = 2, lty = 1, col = cols[[i]])
	    }

	    legend("topleft", legend = popLabs, col=cols, lty = 1, lwd = 2, pch = 16)
		plot(NA, xlim = xlims, ylim = ylims5, xaxt = "n", xlab = "generation", ylab = "S", main = "Selection differential across generations in RCRS")
	    axis(1, at = c(0, RCRSyr), labels = c(0, yr))
	    
	    for (i in 1:length(sRCRS)) {
	    	lines(RCRSgen[-1], sRCRS[[i]], type = "l", lwd = 2, lty = 1, col = cols[[i]])
	    }	
	    legend("topright", legend = popLabs, col=cols, lty = 1, lwd = 2, pch = 16)

	   	plot(NA, xlim = xlims, ylim = ylims6, xaxt = "n", xlab = "generation", ylab = "S/Vg", main = "Selection intensity across generations in RCRS")
	    axis(1, at = c(0, RCRSyr), labels = c(0, yr))
	    
	    for (i in 1:length(iRCRS)) {
	    	lines(RCRSgen[-1], iRCRS[[i]], type = "l", lwd = 2, lty = 1, col = cols[[i]])
	    }
	    legend("topleft", legend = popLabs, col=cols, lty = 1, lwd = 2, pch = 16)
	}
}

