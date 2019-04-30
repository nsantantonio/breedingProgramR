
# save x quantile value!!!!!!!!!!!! check
# burn in RGSC random crosses CHECK
# use mu + x sigma inside RGSC too! --> mu or mu + x sigma or x sigma # check?
# add truely traditional VDP -- how to do? CHECK

# x drops as VDP gets small -- CDF of normal distribution. 
# x as function of heritability and VDP size/intensity
# VDP size relative to the heritability for fixed x
# selection efficiency of VDP

# add sparse selection



# paramL = defArgs; simParam <- SP; verbose = TRUE; checkParam = FALSE; GSfunc = RRBLUP
sim <- function(founderPop, paramL, simParam = SP, returnFunc = identity, verbose = TRUE, checkParam = FALSE, GSfunc = RRBLUP, ...){
	# parameter checks and warnings.
	if (checkParam){
		paramNames <- c("founderRData", "simFunc", "nThreads", "simName", "RGSCselect", "selF2", "nF2", 
						"selQuantile", "ssd", "simpleFounder", "nFounder", "nNuclear", "nChrom", "nLoci", 
						"nM", "nQTL", "Vg", "Vgxe", "founderh2", "h2", "nYr", "nFam", "famSize", "selectTrials", 
						"trialReps", "trialLocs", "cyclePerYr", "returnVDPtoRGSC", "lgen", "RGSCintensity", "reps")
		if (!all(paramNames %in% names(paramL))) stop("not all parameters in 'paramL'! Please include all these parameters in parameter list:\n", paste0(paramNames, "\n"))
	}
	for (p in names(paramL)) assign(p, paramL[[p]])
	
	if(!(length(selectRGSC) == 1 | length(selectRGSC) == nYr)) stop("'selectRGSC' must be of length 1 or nYr!")
	if(all(selectRGSC <= 1)) {
		selectRGSC <- nNuclear * selectRGSC
		if(selectRGSC %% 1 != 0) selectRGSC <- round(selectRGSC)
		deltaRGSC <- if(length(selectRGSC) == nYr) TRUE else FALSE
	}

	if (selF2 & cyclePerYr > 1) warning("Selection on F2 is being performed, and there is more than 1 GS cycle per year. You may want to reduce 'cyclePerYr' to 1")
	# if (selF2 & is.null(selFunc)) warning("Selection on F2 is being performed, but not on expected quantiles, you probably want to set selQuantile = TRUE")

	# check selectTrials & nReturnVDPtoRGSC
	if (!all(selectTrials > 0) | (any(selectTrials < 1) & any(selectTrials > 1))) stop("'selectTrials' must have elements between 0 and 1 or positive integers > 0")
	if (!all(returnVDPtoRGSC >= 0) | (any(returnVDPtoRGSC < 1) & any(returnVDPtoRGSC > 1))) stop("'returnVDPtoRGSC' must have elements between 0 and 1 or positive integers")

	# define number of individuals per cycle, and number to select at each stage
	nI <- nFam * famSize
	if (all(selectTrials <= 1) & !all(selectTrials == 1)) selectTrials <- nI * cumprod(selectTrials) # note this does not allow all to be exactly 1

	if (any(selectTrials %% 1 != 0)){
		selectTrials <- round(selectTrials)
		actInt <- selectTrials / c(nI, selectTrials[-length(selectTrials)])
		cat("NOTE: Selection intensities have been rounded to the nearest integer:\n", selectTrials, "\nThese correspond to selection intensities of:\n", actInt, "\n")
	}
	if (all(returnVDPtoRGSC <= 1)) returnVDPtoRGSC <- returnVDPtoRGSC * c(nI, selectTrials) 

	# count and rename trials
	nTrial <- length(selectTrials)
	trials <- c(paste0("trial", 1:nTrial), "variety")
	if (!is.null(skip)) skip <- trials[skip]

	if(traditional) {
		# selectOut <- "pheno" # can still use markers in traditional program...
		cyclePerYr <- 1
	}
	# if use true values
	if(useTrue){
		cat("NOTICE: using true breeding values for all selection methods.")
		selectOut <- "bv"
		selectIn <- "bv"
		gen0use <- "bv"
		estIntFunc <- bv
		returnVDPcrit <- "bv"
		selectVDP <- "bv"
		useGS <- "bv"
		pullGenoFunc <- pullQtlGeno
		ebv <- bv
	} else {
		gen0use <- "pheno"
		estIntFunc <- pheno
		useGS <- "pheno"
		pullGenoFunc <- pullSnpGeno
	}

	# define cycles
	cycle <- 1:cyclePerYr

	# initialize lists to store populations
	RGSC <- list()
	GSmodel <- list()
	predAcc <- list()
	intensity <- list()
	names(trials) <- trials
	VDP <- lapply(trials, function(x) list())
	
	# initialize nuclear population, train GS model (necessary?), predict ebv (note the ebv's should be bad if markers are not exactly on QTL) 
	RGSC[[gen(0)]] <- newPop(founderPop)
	if(verbose) cat("Founder population has genetic variance of:", popVar(getTrueBV(RGSC[[gen(0)]], simParam = simParam)), "\n")
	
	# burn in using random crosses
	printBurnin <- TRUE
	while (nFam > nInd(RGSC[[gen(0)]]) | founderBurnIn > 0) {
		if(nFam > nInd(RGSC[[gen(0)]]) & verbose) cat("nFounder < nFam. Random mating to make nFam parents...\n")
		if(printBurnin & founderBurnIn > 0 & verbose) cat("Running", founderBurnIn, "burn-in cycles of random mating...\n")
		RGSC[[gen(0)]] <- selectCross(RGSC[[gen(0)]], nInd = nInd(RGSC[[gen(0)]]), use = "rand", simParam = simParam, nCrosses = max(nFam, nInd(RGSC[[gen(0)]]))) 
		founderBurnIn <- founderBurnIn - 1
		if(printBurnin) printBurnin <- FALSE
	}
	# train GS model and set EBV (after selfing if ssd)
	GSmodel[[gen(0)]] <- GSfunc(RGSC[[gen(0)]], traits = 1, use = gen0use, snpChip = 1, simParam = simParam, useQTL = useTrue)
	if (selF2) RGSC[[gen(0)]] <- self(RGSC[[gen(0)]], nProgeny = nF2, simParam = simParam)
	RGSC[[gen(0)]] <- setEBV(RGSC[[gen(0)]], GSmodel[[gen(0)]], simParam = simParam)
	
	xInt <- if(is.null(setXint)) 1 - selectTrials[nTrial] / selectTrials[1] else setXint 
	if(verbose) cat("Estimated selection intensity:", qnorm(xInt, sd = sqrt(Vg)), "\n")

	# pop <- RGSC[[gen(0)]]; GSfit <- GSmodel[[gen(0)]]
	# pop <- RGSC[[lastRGSCgen]]; GSfit <- GSmodel[[lastGSmodel]]

	# run program for nYr years
	for (i in 1:(nYr + nTrial - 1)) { 
	# for (i in 1:6) { 
		# i = 1
		lastRGSCgen <- names(RGSC)[length(RGSC)]
		lastGSmodel <- if (i <= nYr) gen(i-1) else gen(nYr)
		selectRGSCi <- if(deltaRGSC) selectRGSC[i] else selectRGSC[1]
		# nProgenyPerCrossIni <- if(identical(expDistPairs, selFuncIn)) nNuclear / selectRGSCi  * nProgenyPerCrossIn else nProgenyPerCrossIn

		if (i <= nYr){ 
			if (verbose) cat("    Year: ", i, "\n")

			# predict latest RGSC with updated GS model 
			if (i > 1) RGSC[[lastRGSCgen]] <- setEBV(RGSC[[lastRGSCgen]], GSmodel[[lastGSmodel]], simParam = simParam)
			
			# estimate VDP selection intensity from previous generations
			if (i > nTrial) {
				intensity[[gen(i)]] <- estIntensity(VDP, i, nT = nTrial, start = "trial1", end = "variety", estFunc = estIntFunc, Gvar = varA)
				if(is.null(setXint)) xInt <- pnorm(intensity[[gen(i)]]) 
				if(verbose) cat("Realized selection intensity:", intensity[[gen(i)]], "\n")
			}
			
			# use selected lines from last trial to make new crosses. Not sure this timeline is realistic...

			if(traditional & i > 1) {
				lastVDPSel <- tail(names(VDP)[sapply(VDP, length) > 0], 1)
				if(verbose) cat("        ", nInd(RGSC[[lastRGSCgen]]), "lines selected out of VDP", lastVDPSel, "for making crosses \n")
			} else {
				if(verbose) cat("        ", nInd(RGSC[[lastRGSCgen]]), "individuals selected out of", nNuclear, " RGSC population \n")
			}
			if(is.null(selFuncOut)){
				selToP <- selectInd(RGSC[[lastRGSCgen]], nInd = nFam, trait = 1, use = selectOut)
			} else {
			# select out of RGSC, on mean, expected quantile, etc...
				selToP <- do.call(selFuncOut, getArgs(selFuncOut, nSel = nFam, pop = RGSC[[lastRGSCgen]], GSfit = GSmodel[[lastGSmodel]], 
												  trait = 1, use = selectOut, quant = xInt, nProgeny = nProgenyPerCrossOut, 
												  # pullGeno = pullGenoFunc, w = weight, ...))
												  pullGeno = pullGenoFunc, w = weight))
			}
			if(nInd(selToP) != nFam) stop("selToP is wrong...")
			
			# (mean(ebv(selToP)) - mean(ebv(RGSC[[lastRGSCgen]]))) / sd(ebv(RGSC[[lastRGSCgen]]))
			# Determine number of individuals to select within familiy (default is all)
			nProgPerFam <- famSize / withinFamInt
			if ((nProgPerFam) %% 1 != 0) {
				nProgPerFam <- round(nProgPerFam)
				nSelToTrial <- round(nProgPerFam * withinFamInt)
				cat("\nNOTE: Selection intensities within familiy have been rounded to the nearest integer resulting in", nSelToTrial, "progeny per family selected from", nProg, "progeny per family\n")
			}

			# make DH families or self
			VDP[[trials[1]]][[gen(i)]] <- if (ssd) self(selToP, nProgeny = nProgPerFam) else makeDH(selToP, nDH = nProgPerFam)

			#select within family if intensity < 1
			if (withinFamInt < 1) {
				VDP[[trials[1]]][[gen(i)]] <- setEBV(VDP[[trials[1]]][[gen(i)]], GSmodel[[lastGSmodel]], simParam = simParam)
				fams <- split(1:(nProgPerFam * nFam), rep(1:nFam, each = nProgPerFam))
				faml <- list()
				for (j in names(fams)){
					faml[[j]] <- selectInd(VDP[[trials[1]]][[gen(i)]][fams[[j]]], nInd = famSize, trait = 1, use = selectOut) 
				}
				VDP[[trials[1]]][[gen(i)]] <- mergePops(faml)		
			}
			# # print mean genotypic value of DH 
			# if (verbose) print(sapply(VDP[[trials[1]]], function(x) mean(gv(x))))
			if(nInd(VDP[[trials[1]]][[gen(i)]]) != nFam * nProgPerFam) stop("selToP is wrong...")
		}

		# get generation indices
		genI <- tail(1:i, min(5, i))
		genI <- genI[genI <= nYr]
		genBack <- abs(genI - i) + 1
		index = 1:length(genI)

		# phenotype, or predict if skipped
		for (g in index) {
			Vgi <- if (updateVg) varG(VDP[[trials[genBack[g]]]][[gen(genI[g])]])[[1]] else Vg
			if (!trials[genBack[g]] %in% skip) VDP[[trials[genBack[g]]]][[gen(genI[g])]] <- setPheno(VDP[[trials[genBack[g]]]][[gen(genI[g])]], varE = h2toVe(h2[genBack[g]], Vgi), reps = trialReps[genBack[g]] * trialLocs[genBack[g]])
		}

		if (i <= nYr){
			# if(traditional) 
			for (j in cycle){
				if(traditional) {
					lastVDPSel <- tail(names(VDP)[sapply(VDP, length) > 0], 1)
					selPop <- VDP[[lastVDPSel]][[length(VDP[[lastVDPSel]])]] 
				} else {
					if (j != cycle[1]) RGSC[[gen(j-1)]] <- setEBV(RGSC[[gen(j-1)]], GSmodel[[lastGSmodel]], simParam = simParam)
					predAcc[["RGSC"]][[gen(j-1)]] <- getAcc(RGSC[[gen(j-1)]])
					selPop <- RGSC[[gen(j-1)]]
				}
			# run GS model to cycle through RGSC for year i
				if(is.null(selFuncIn)){
					RGSC[[gen(j)]] <- selectCross(pop = selPop, nInd = min(selectRGSCi, nInd(selPop)), use = selectIn,  trait = 1, simParam = simParam, nCrosses = nNuclear, nProgeny = nProgenyPerCrossIn)
				} else {
					# print(selectRGSCi)
					# if(identical(expDistPairs, selFuncIn)) nProgenyPerCrossIn <- nNuclear / selectRGSCi  * nProgenyPerCrossIn 
					RGSC[[gen(j)]] <- do.call(selFuncIn, getArgs(selFuncIn, nSel = selectRGSCi, pop = selPop, GSfit = GSmodel[[lastGSmodel]],
					trait = 1,  use = selectIn,  trait = 1, nCrosses = nNuclear, nProgeny = nProgenyPerCrossIn, quant = xInt, verbose = verbose, 
					# pullGeno = pullGenoFunc, w = weight, ...))
					pullGeno = pullGenoFunc, w = weight))
					if(nInd(RGSC[[gen(j)]]) != nNuclear) stop("nNuclear isnt right...")
				}
				# would be good to be able to select within f2 family if f2 > 1
				if (selF2) RGSC[[gen(j)]] <- self(RGSC[[gen(j)]], nProgeny = nF2, simParam = simParam)
			}
			# update GScycle number
			cycle <- cycle + cyclePerYr
			
			# so this removes the selections from the training population. 
			trnSet <- lapply(VDP[trials[!grepl("variety", trials)]], function(x) x[names(x) %in% gen(max(1, i-max(1, lgen)):i)])
			trnSet <- trnSet[sapply(trnSet, length) > 0]

			# concatenate training set and train GS model
	 		train <- mergePopsRec(trnSet) 
			cat("training set has ", train@nInd, "individuals...\n")	
			GSmodel[[gen(i)]] <- GSfunc(train, traits = 1, use = useGS, snpChip = 1, simParam=simParam, useQTL = useTrue)
		} 
		# check if final year
		if (i - nYr == 1) cat("\nFinal year reached, selecting on phenotypes / ebv trained with last year training set ...\n")
		
		for (g in index) {
			# set ebv if using to select / skip generations
			if (selectVDP == "ebv" | !is.null(skip)) {
				GSmodelVDP <- if (trials[genBack[g]] %in% skip) lastGSmodel else gen(min(i, nYr))
				VDP[[trials[genBack[g]]]][[gen(genI[g])]] <- setEBV(VDP[[trials[genBack[g]]]][[gen(genI[g])]], GSmodel[[GSmodelVDP]], simParam = simParam)
				predAcc[[trials[genBack[g]]]][[gen(genI[g])]] <- getAcc(VDP[[trials[genBack[g]]]][[gen(genI[g])]])
			}

			# select based on ebv or phenotype
			sel <- if (trials[genBack[g]] %in% skip) "ebv" else  selectVDP
			if (i - genI[g] < nTrial) VDP[[trials[genBack[g] + 1]]][[gen(genI[g])]] <- selectInd(VDP[[trials[genBack[g]]]][[gen(genI[g])]], nInd = selectTrials[genBack[g]], trait = 1, use = sel, returnPop = TRUE)
			if (ssd) VDP[[trials[genBack[g] + 1]]][[gen(genI[g])]] <- self(VDP[[trials[genBack[g] + 1]]][[gen(genI[g])]])
		}

		if (i <= nYr){
			# return lines from VDP into the RGSC 
			if (any(returnVDPtoRGSC > 0)){
				returnToRGSC <- genBack %in% which(returnVDPtoRGSC > 0)
				if (sum(returnToRGSC) > 0){	
					addToRGSC <- list()
					for (g in index[returnToRGSC]) {
						addToRGSC[[gen(g)]] <- selectInd(VDP[[trials[genBack[g]]]][[gen(genI[g])]], nInd = returnVDPtoRGSC[genBack[g]], trait = 1, use = returnVDPcrit, returnPop = TRUE) 
					}
					if (length(addToRGSC) > 0){
						addToRGSC <- Reduce(c, addToRGSC)
						RGSC[[gen(cycle[1] - 1)]] <- c(RGSC[[gen(cycle[1] - 1)]], addToRGSC)
					}
				}
			}
		}
	}

	rL <- list(SP = SP, paramL = paramL, RGSC = RGSC, VDP = VDP, GSmodel = GSmodel, predAcc = predAcc, selIntensity = intensity)
	do.call(returnFunc, getArgs(returnFunc, resultL = rL, ...))
}
