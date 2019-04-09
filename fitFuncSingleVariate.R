# paramL = defArgs; simParam <- SP; verbose = TRUE; checkParam = FALSE
sim <- function(founderPop, paramL, simParam = SP, returnFunc = identity, verbose = TRUE, checkParam = FALSE){
	# parameter checks and warnings.
	if(checkParam){
		paramNames <- c("founderRData", "simFunc", "nThreads", "simName", "RGSCselect", "selF2", "nF2", 
						"selQuantile", "ssd", "simpleFounder", "nFounder", "nNuclear", "nChrom", "nLoci", 
						"nM", "nQTL", "Vg", "Vgxe", "founderh2", "h2", "nYr", "nFam", "famSize", "selectTrials", 
						"trialReps", "trialLocs", "GScylcePerYr", "returnVDPtoRGSC", "lgen", "nGen", "RGSCintensity", "reps")
		if(!all(paramNames %in% names(paramL))) stop("not all parameters in 'paramL'! Please include all these parameters in parameter list:\n", paste0(paramNames, "\n"))
	}
	for(p in names(paramL)) assign(p, paramL[[p]])
	
	if(selF2 & GScylcePerYr > 1) warning("Selection on F2 is being performed, and there is more than 1 GS cycle per year. You may want to reduce 'GScylcePerYr' to 1")
	# if(selF2 & is.null(selFunc)) warning("Selection on F2 is being performed, but not on expected quantiles, you probably want to set selQuantile = TRUE")

	# check selectTrials & nReturnVDPtoRGSC
	if(!all(selectTrials > 0) | (any(selectTrials < 1) & any(selectTrials > 1))) stop("'selectTrials' must have elements between 0 and 1 or positive integers > 0")
	if(!all(returnVDPtoRGSC >= 0) | (any(returnVDPtoRGSC < 1) & any(returnVDPtoRGSC > 1))) stop("'returnVDPtoRGSC' must have elements between 0 and 1 or positive integers")

	# define number of individuals per cycle, and number to select at each stage
	nI <- nFam * famSize
	if(all(selectTrials <= 1) & !all(selectTrials == 1)) selectTrials <- nI * cumprod(selectTrials) # note this does not allow all to be exactly 1

	if(any(selectTrials %% 1 != 0)){
		selectTrials <- round(selectTrials)
		actInt <- selectTrials / c(nI, selectTrials[-length(selectTrials)])
		cat("NOTE: Selection intensities have been rounded to the nearest integer:\n", selectTrials, "\nThese correspond to selection intensities of:\n", actInt, "\n")
	}
	if(all(returnVDPtoRGSC <= 1)) returnVDPtoRGSC <- returnVDPtoRGSC * c(nI, selectTrials) 

	# count and rename trials
	nTrial <- length(selectTrials)
	trials <- c(paste0("trial", 1:nTrial), "variety")
	if(!is.null(skip)) skip <- trials[skip]

	# define cycles
	GScylce <- 1:GScylcePerYr

	# initialize lists to store populations
	RGSC <- list()
	GSmodel <- list()
	predAcc <- list()
	names(trials) <- trials
	VDP <- lapply(trials, function(x) list())

	# initialize nuclear population, train GS model (necessary?), predict ebv (note the ebv's should be bad if markers are not exactly on QTL) 
	RGSC[[gen(0)]] <- newPop(founderPop)
	if(nFam > nInd(RGSC[[gen(0)]])) RGSC[[gen(0)]] <- selectCross(RGSC[[gen(0)]], nInd = nInd(RGSC[[gen(0)]]), use = "rand", simParam = simParam, nCrosses = nFam) 
	GSmodel[[gen(0)]] <- RRBLUP(RGSC[[gen(0)]], traits = 1, use = "pheno", snpChip = 1, simParam = simParam)
	if(selF2) RGSC[[gen(0)]] <- self(RGSC[[gen(0)]], nProgeny = nF2, simParam = simParam)
	RGSC[[gen(0)]] <- setEBV(RGSC[[gen(0)]], GSmodel[[gen(0)]], simParam = simParam)
	
	# run program for nYr years
	for(i in 1:(nYr + nTrial - 1)) { 
		lastRGSCgen <- names(RGSC)[length(RGSC)]
		lastGSmodel <- if(i <= nYr) gen(i-1) else gen(nYr)
		# nFami <- min(nInd(RGSC[[lastRGSCgen]]), nFam)
		if(i <= nYr){ # I need to fix this! the index still tries to set phenotypes for generations that dont exist
			if(verbose) cat("Year: ", i, "\n")
			# i = 1
			if(i > 1) {
				# predict latest RGSC with updated GS model 
				RGSC[[lastRGSCgen]] <- setEBV(RGSC[[lastRGSCgen]], GSmodel[[lastGSmodel]], simParam = simParam)
			}
			# select out of RGSC, on mean, expected quantile, etc...
			if(!is.null(selFunc)) {
				if(identical(selFunc, getExpDist)) {
					if(i > nTrial) {
						intensity <- (mean(gv(VDP[["variety"]][[gen(i - nTrial)]])) - mean(gv(VDP[["trial1"]][[gen(i - nTrial)]]))) / sqrt(varA(VDP[["trial1"]][[gen(i - nTrial)]])[[1]])
						# i should probably use ebv instead. need to set ebv for varieties... or use pheno. 
						# Also, should note that pheno of varieties will always be biased upward due to select on error?
						qInt <- pnorm(intensity)
					} else {
						qInt <- 1 - selectTrials[nTrial] / selectTrials[1] 
					} 
				} else {
					qInt <- useQuantile
				}
				selCrit <- selFunc(RGSC[[lastRGSCgen]], GSmodel[[lastGSmodel]], quant = qInt)
				parSel <- getSel(selCrit, nFam)
				if(is.matrix(parSel)){
					if(ncol(parSel) == 2){
						# if(famSize > 1) parSel <- parSel[rep(1:nFam, each = famSize), ] # dont think this is necessary, you dont need to make the cropss more than once, just DH/self. 
						selGStoP <- makeCross(RGSC[[lastRGSCgen]], crossPlan = parSel) # perhaps this should be set to the previous RGSC generation??
					} else {
						stop("Something is wrong with parent selection. Expecting 2 columns, p1 and p2 indicating parent pairs.")
					}
				} else {
					selGStoP <- RGSC[[lastRGSCgen]][parSel]
					if(any(!selGStoP@id %in% parSel)) stop("parent selection out of RGSC failed!")
				}
			} else {
				selGStoP <- selectInd(RGSC[[lastRGSCgen]], nInd = nFam, trait = 1, use = selectOutRGSC) 
			}

			# Determine number of individuals to select within familiy (default is all)
			nProgPerFam <- famSize / withinFamInt
			if((nProgPerFam) %% 1 != 0) {
				nProgPerFam <- round(nProgPerFam)
				nSelToTrial <- round(nProgPerFam * withinFamInt)
				cat("\nNOTE: Selection intensities within familiy have been rounded to the nearest integer resulting in", nSelToTrial, "progeny per family selected from", nProg, "progeny per family\n")
			}

			# make DH families or self
			VDP[[trials[1]]][[gen(i)]] <- if(ssd) self(selGStoP, nProgeny = nProgPerFam) else makeDH(selGStoP, nDH = nProgPerFam)

			#select within family if intensity < 1
			if(withinFamInt < 1) {
				VDP[[trials[1]]][[gen(i)]] <- setEBV(VDP[[trials[1]]][[gen(i)]], GSmodel[[lastGSmodel]], simParam = simParam)
				fams <- split(1:(nProgPerFam * nFam), rep(1:nFam, each = nProgPerFam))
				faml <- list()
				for(j in names(fams)){
					faml[[j]] <- selectInd(VDP[[trials[1]]][[gen(i)]][fams[[j]]], nInd = famSize, trait = 1, use = selectOutRGSC) 
				}
				VDP[[trials[1]]][[gen(i)]] <- mergePops(faml)		
			}
			# # print mean genotypic value of DH 
			# if(verbose) print(sapply(VDP[[trials[1]]], function(x) mean(gv(x))))
		}

		# get generation indices
		genI <- tail(1:i, min(5, i))
		genI <- genI[genI <= nYr]
		genBack <- abs(genI - i) + 1
		index = 1:length(genI)

		# phenotype, or predict if skipped
		for(g in index) {
			Vgi <- if(updateVg) varG(VDP[[trials[genBack[g]]]][[gen(genI[g])]])[[1]] else Vg
			if(!trials[genBack[g]] %in% skip) VDP[[trials[genBack[g]]]][[gen(genI[g])]] <- setPheno(VDP[[trials[genBack[g]]]][[gen(genI[g])]], varE = h2toVe(h2[genBack[g]], Vgi), reps = trialReps[genBack[g]] * trialLocs[genBack[g]])
		}

		if(i <= nYr){
			# run GS model to cycle through RGSC for year i
			for(j in GScylce){
				if(j != GScylce[1]) RGSC[[gen(j-1)]] <- setEBV(RGSC[[gen(j-1)]], GSmodel[[lastGSmodel]], simParam = simParam)
				predAcc[["RGSC"]][[gen(j-1)]] <- getAcc(RGSC[[gen(j-1)]])
				RGSC[[gen(j)]] <- selectCross(pop = RGSC[[gen(j-1)]], nInd = RGSC[[gen(j-1)]]@nInd * RGSCintensity, 
											   use = selectInRGSC,  trait = 1, simParam = simParam, nCrosses = nNuclear, nProgeny = 1) 
				if(selF2) RGSC[[gen(j)]] <- self(RGSC[[gen(j)]], nProgeny = nF2, simParam = simParam)
			}
			# update GScycle number
			GScylce <- GScylce + GScylcePerYr
			
			# so this removes the selections from the training population. 
			trnSet <- lapply(VDP[trials[!grepl("variety", trials)]], function(x) x[names(x) %in% gen(max(1, i-max(1, lgen)):i)])
			trnSet <- trnSet[sapply(trnSet, length) > 0]

			# concatenate training set and train GS model
	 		train <- mergePopsRec(trnSet) 
			cat("training set has ", train@nInd, "individuals...\n")	
			GSmodel[[gen(i)]] <- RRBLUP(train, traits = 1, use = "pheno", snpChip = 1, simParam=simParam)
		} 
		# check if final year
		if (i - nYr == 1) cat("\nFinal year reached, selecting on phenotypes / ebv trained with last year training set ...\n")
		
		for(g in index) {
			# set ebv if using to select / skip generations
			if(selectVDP == "ebv" | !is.null(skip)) {
				GSmodelVDP <- if(trials[genBack[g]] %in% skip) lastGSmodel else gen(min(i, nYr))
				VDP[[trials[genBack[g]]]][[gen(genI[g])]] <- setEBV(VDP[[trials[genBack[g]]]][[gen(genI[g])]], GSmodel[[GSmodelVDP]], simParam = simParam)
				predAcc[[trials[genBack[g]]]][[gen(genI[g])]] <- getAcc(VDP[[trials[genBack[g]]]][[gen(genI[g])]])
			}

			# select based on ebv or phenotype
			sel <- if(trials[genBack[g]] %in% skip) "ebv" else  selectVDP
			if(i - genI[g] < nTrial) VDP[[trials[genBack[g] + 1]]][[gen(genI[g])]] <- selectInd(VDP[[trials[genBack[g]]]][[gen(genI[g])]], nInd = selectTrials[genBack[g]], trait = 1, use = sel, returnPop = TRUE)
			if(ssd) VDP[[trials[genBack[g] + 1]]][[gen(genI[g])]] <- self(VDP[[trials[genBack[g] + 1]]][[gen(genI[g])]])
		}

		if(i <= nYr){
			# return lines from VDP into the RGSC 
			if(any(returnVDPtoRGSC > 0)){
				returnToRGSC <- genBack %in% which(returnVDPtoRGSC > 0)
				if(sum(returnToRGSC) > 0){	
					addToRGSC <- list()
					for(g in index[returnToRGSC]) {
						addToRGSC[[gen(g)]] <- selectInd(VDP[[trials[genBack[g]]]][[gen(genI[g])]], nInd = returnVDPtoRGSC[genBack[g]], trait = 1, use = returnVDPcrit, returnPop = TRUE) 
					}
					if(length(addToRGSC) > 0){
						addToRGSC <- Reduce(c, addToRGSC)
						RGSC[[gen(GScylce[1] - 1)]] <- c(RGSC[[gen(GScylce[1] - 1)]], addToRGSC)
					}
				}
			}
		}
	}

	rL <- list(SP = SP, paramL = paramL, RGSC = RGSC, VDP = VDP, GSmodel = GSmodel, predAcc = predAcc)
	returnFunc(rL)
}
