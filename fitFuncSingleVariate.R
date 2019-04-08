

# makeCrossPlan <- function(crossPlanFunc, pop, nInd){


defArgs <- list(seed = NULL,
			   founderRData = "founderPop/testAlphaSimR1000SegSite.RData",
			   simFunc = "fitFuncSingleVariate.R", 
			   nThreads = 10, 
			   simName = "simNoTitle",
			   selectInRGSC = "ebv", # ebv, rand
			   selectOutRGSC = "ebv", # ebv, var, exp?
			   selectInFam = "none", # none, 
			   selectVDP = "pheno", # ebv, pheno
			   returnVDPcrit = "pheno", # ebv?
			   kinship = "SNP",
			   updateVg = FALSE,
			   selF2 = FALSE,
			   nF2 = 10,
			   selFunc = getExpDist,
			   ssd = FALSE,
			   simpleFounder = FALSE,
			   nFounder = 10,
			   nNuclear = 20,
			   nFam = 10 ,
			   famSize = 10,
			   nChrom = 10,
			   nLoci = 100,
			   nM = floor(c(10*(9:1), 4) / 4),
			   nQTL = c(2*(9:1), 1),
			   Vg = 1,
			   Vgxe = 1,
			   founderh2 = 0.3,
			   h2 = c(0.1, 0.3, 0.3, 0.3, 0.3),
			   nYr = 7,
			   selectTrials = c(0.5, 0.5, 0.5, 0.5, 0.5),
			   trialReps = c(1, 1, 2, 3, 3),
			   trialLocs = c(1, 1, 2, 5, 5),
			   GScylcePerYr = 2 ,
			   returnVDPtoRGSC = c(0, 0.5, 0, 0, 0, 1), # default to rep(0, nTrial)?
			   lgen = 5,
			   nGen = 20,
			   RGSCintensity = 0.2,
			   reps = 10,
			   skip = NULL
)

# paramL = defArgs; simParam <- SP; select = "pheno"; returnFunc = identity; verbose = TRUE; skip = NULL; selQuantile = TRUE; checkParam = FALSE
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
	if(selF2 & !selQuantile) warning("Selection on F2 is being performed, but not on expected quantiles, you probably want to set selQuantile = TRUE")

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

	# if(!is.null(returnVDPtoRGSC)) if(!all(returnVDPtoRGSC %in% trials)) stop("somthing is wrong with returnVDPtoRGSC...") 

	# define cycles
	GScylce <- 1:GScylcePerYr

	# initialize lists to store populations
	RGSC <- list()
	GSmodel <- list()
	predAcc <- list()
	names(trials) <- trials
	VDP <- lapply(trials, function(x) list())
	# train <- NULL

	# initialize nuclear population, train GS model (necessary?), predict ebv (note the ebv's should be bad if markers are not exactly on QTL) 
	RGSC[[gen(0)]] <- newPop(founderPop)
	GSmodel[[gen(0)]] <- RRBLUP(RGSC[[gen(0)]], traits = 1, use = "pheno", snpChip = 1, simParam = simParam)
	if(selF2) RGSC[[gen(0)]] <- self(RGSC[[gen(0)]], nProgeny = nF2, simParam = simParam)
	RGSC[[gen(0)]] <- setEBV(RGSC[[gen(0)]], GSmodel[[gen(0)]], simParam = simParam)

	# run program for nYr years
	for(i in 1:(nYr + nTrial - 1)) { 
	# for(i in 1:6) { 
		lastRGSCgen <- names(RGSC)[length(RGSC)]
		if(i <= nYr){ # I need to fix this! the index still tries to set phenotypes for generations that dont exist
			if(verbose) cat("Year:", i, "\n")
			# i = 1
			if(i > 1) {
				# predict latest RGSC with updated GS model 
				# RGSC[[gen(GScylce[1]-1)]] <- setEBV(RGSC[[gen(GScylce[1]-1)]], GSmodel[[gen(i-1)]], simParam = simParam)
				RGSC[[lastRGSCgen]] <- setEBV(RGSC[[lastRGSCgen]], GSmodel[[gen(i-1)]], simParam = simParam)
			}
			# select out of RGSC, on mean, expected quantile, etc...
			if(!is.null(selFunc)) {
				if(identical(selFunc, getExpDist)) {
					if(i > nTrial) {
						intensity <- (mean(gv(VDP[["variety"]][[gen(i - 5)]])) - mean(gv(VDP[["trial1"]][[gen(i - 5)]]))) / sqrt(varA(VDP[["trial1"]][[gen(i - 5)]])[[1]])
						# i should probably use ebv instead. need to set ebv for varieties... or use pheno. 
						# Also, should note that pheno of varieties will always be biased upward due to select on error?
						# intensity <- (mean(ebv(VDP[["variety"]][[gen(i - 5)]])) - mean(ebv(VDP[["trial1"]][[gen(i - 5)]]))) / sqrt(varA(VDP[["trial1"]][[gen(i - 5)]])[[1]])
						qInt <- pnorm(intensity)
					} else {
						qInt <- 1 - selectTrials[nTrial] / selectTrials[1] 
					} 
				} else {
					qInt <- 0.9
				}
				selCrit <- selFunc(RGSC[[lastRGSCgen]], GSmodel[[gen(i-1)]], quant = qInt)
				parSel <- getSel(selCrit, nFam)
				if(is.matrix(parSel)){
					if(ncol(parSel) == 2){
						if(famSize > 1) parSel <- parSel[rep(1:nFam, each = famSize), ]
						selGStoP <- makeCross(RGSC[[lastRGSCgen]], crossPlan = parSel)
					} else {
						stop("Something is wrong with parent selection. Expecting 2 columns, p1 and p2 indicating parent pairs.")
					}
				} else {
					selGStoP <- RGSC[[lastRGSCgen]][parSel]
				}
				if(any(!selGStoP@id %in% parSel)) stop("parent selection out of RGSC failed!")
				# selGStoP <- selectInd(RGSC[[length(RGSC)]], nInd = nFam, trait = dummyFunc, use = , retrn = expQuant) 
			} else {
				selGStoP <- selectInd(RGSC[[lastRGSCgen]], nInd = nFam, trait = 1, use = selectOutRGSC) 
			}
			# make DH families or self
			VDP[[trials[1]]][[gen(i)]] <- if(ssd) self(selGStoP, nProgeny = famSize) else makeDH(selGStoP, nDH = famSize)
			# print mean genotypic value of DH 
			if(verbose) print(sapply(VDP[[trials[1]]], function(x) mean(gv(x))))
		}

		# get generation indices
		genI <- tail(1:i, min(5, i))
		genI <- genI[genI <= nYr]
		genBack <- abs(genI - i) + 1
		index = 1:length(genI)
		for(g in index) {
			gi <- genI[g]
			gb <- genBack[g]
			ti <- trials[gb]
			#phenotype if not skipped
			Vgi <- if(updateVg) varG(VDP[[ti]][[gen(gi)]])[[1]] else Vg
			if(!ti %in% skip) VDP[[ti]][[gen(gi)]] <- setPheno(VDP[[ti]][[gen(gi)]], varE = h2toVe(h2[gb], Vgi), reps = trialReps[gb] * trialLocs[gb])

			# set ebv (does this use phenotypes if not set above? need to check...), Yes if those phenotypes were in the trainning pop for RRBLUP fit. 
			if(select == "ebv" | !is.null(skip)) {
				VDP[[ti]][[gen(gi)]] <- setEBV(VDP[[ti]][[gen(gi)]], GSmodel[[gen(i-1)]], simParam = simParam)
				predAcc[[ti]][[gen(gi)]] <- getAcc(VDP[[ti]][[gen(gi)]])
			}

			# NEED TO MOVE SELECTION TILL AFTER GS MODEL UPDATED!!! why? they should get removed
			# select indviduals for next years trial based on ebv and/or phenotype
			sel <- if(ti %in% skip) "ebv" else  select
			if(i - gi < nTrial) VDP[[trials[gb + 1]]][[gen(gi)]] <- selectInd(VDP[[ti]][[gen(gi)]], nInd = selectTrials[gb], trait = 1, use = sel, returnPop = TRUE)
			if(ssd) VDP[[trials[gb + 1]]][[gen(gi)]] <- self(VDP[[trials[gb + 1]]][[gen(gi)]])
		}

		if(i <= nYr){
			# run GS model to cycle through RGSC for year i
			for(j in GScylce){
				if(j != GScylce[1]) RGSC[[gen(j-1)]] <- setEBV(RGSC[[gen(j-1)]], GSmodel[[gen(i-1)]], simParam = simParam)
				predAcc[["RGSC"]][[gen(j-1)]] <- getAcc(RGSC[[gen(j-1)]])
				RGSC[[gen(j)]] <- selectCross(pop = RGSC[[gen(j-1)]], nInd = RGSC[[gen(j-1)]]@nInd * RGSCintensity, 
											   use = selectInRGSC,  trait = 1, simParam = simParam, nCrosses = nNuclear, nProgeny = 1) 
				if(selF2) RGSC[[gen(j)]] <- self(RGSC[[gen(j)]], nProgeny = nF2, simParam = simParam)
			}
			# update GScycle number
			GScylce <- GScylce + GScylcePerYr

			# return lines from VDP into the RGSC 
			if(any(returnVDPtoRGSC > 0)){
				returnToRGSC <- genBack %in% which(returnVDPtoRGSC > 0)
				if(sum(returnToRGSC) > 0){	
					addToRGSC <- list()
					for(g in index[returnToRGSC]) {
						gi <- genI[g]
						gb <- genBack[g]
						ti <- trials[gb]
						# returns selection out of selection
						addToRGSC[[gen(g)]] <- selectInd(VDP[[ti]][[gen(gi)]], nInd = returnVDPtoRGSC[gb], trait = 1, use = returnVDPcrit, returnPop = TRUE) 
					}
					if(length(addToRGSC) > 0){
						addToRGSC <- Reduce(c, addToRGSC)
						RGSC[[gen(GScylce[1] - 1)]] <- c(RGSC[[gen(GScylce[1] - 1)]], addToRGSC)
					}
				}
			}
			
			# so this removes the selections from the training population. 
			trnSet <- lapply(VDP[trials[!grepl("variety", trials)]], function(x) x[names(x) %in% gen(max(1, i-max(1, lgen)):i)])
	 		trnSet[-1] <- lapply(trnSet[-1], function(x) x[-length(x)])
			trnSet <- trnSet[sapply(trnSet, length) > 0]

			# concatenate training set and train GS model
	 		# train <- Reduce(c, lapply(trnSet, function(x) Reduce(c, x)))
	 		train <- mergePops(lapply(trnSet, mergePops)) # apparently there is a function to do this...
	 		# train <- rlapply(trnSet, level = 2, f = mergePops, combine = mergePops) # This should work too....

			if(is.list(train)) train <- Reduce(c, train)
			cat("training set has ", train@nInd, "individuals...\n")	
			GSmodel[[gen(i)]] <- RRBLUP(train, traits = 1, use = "pheno", snpChip = 1, simParam=simParam)
		} else {
			if (i - nYr == 1) cat("Final year reached, selecting on phenotypes / ebv trained with last year training set ...\n")
			GSmodel[[gen(i)]] <- GSmodel[[gen(i-1)]]
		}
	}

	rL <- list(SP = SP, paramL = paramL, RGSC = RGSC, VDP = VDP, GSmodel = GSmodel, predAcc = predAcc)
	returnFunc(rL)
}

# simDHdist <- function(pop, returnQuantile = 0.9){

# }

