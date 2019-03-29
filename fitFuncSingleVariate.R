# simParam <- SP; select = "ebv"; returnFunc = identity; verbose = TRUE; skip = NULL; selDHquantile = TRUE
sim <- function(founderPop, simParam = SP, select = "ebv", returnFunc = identity, verbose = TRUE, skip = NULL){
	
	simDHdist <- function(pop, DHsampleSize = 200, GSfit = GSmodel[[1]], returnQuantile = 0.9){
		qDH <- function(i){
			DH <- makeDH(pop[i], nDH = DHsampleSize)
			DH <- setEBV(DH, GSfit, simParam = simParam)
			DHebv <- ebv(DH)
			quantile(DHebv, returnQuantile)[[1]]
		}
		sapply(pop0@id, qDH)
	}

	if(!all(selectTrials > 0) | (any(selectTrials < 1) & any(selectTrials > 1))) stop("'selectTrials' must have elements between 0 and 1 or positiive integers")
	
	nDH <- nDHfam * DHfamSize
	if(all(selectTrials < 1)) selectTrials <- nDH * cumprod(selectTrials)

	nTrial <- length(selectTrials)
	trials <- c(paste0("trial", 1:nTrial), "variety")
	returnVDPtoRGSC <- trials[returnVDPtoRGSC]
	if(!is.null(skip)) skip <- trials[skip]
	# trials <- c("headrow", "prelim", "advance", "elite1", "elite2", "variety")

	if(!is.null(returnVDPtoRGSC)) if(!all(returnVDPtoRGSC %in% trials[-length(trials)])) stop("'returnVDPtoRGSC' argument can only take values of 'headrow', 'prelim', 'advance', 'elite1', or 'elite2', to return selected lines out of those trials into the RGSC, e.g. to return identified varieties, use 'elite2'") 
	pop0 <- newPop(founderPop)
# i commented this out, not sure it needs to be set?
	# pheno(pop0)
	GScylce <- 1:GScylcePerYr

	RGSC <- list()
	GSmodel <- list()
	names(trials) <- trials
	VDP <- lapply(trials, function(x) list())

	# initialize nuclear population
	RGSC[[gen(0)]] <- pop0
	GSmodel[[gen(0)]] <- RRBLUP(RGSC[[gen(0)]], traits = 1, use = "pheno", snpChip = 1, simParam = simParam)
	RGSC[[gen(0)]] <- setEBV(RGSC[[gen(0)]], GSmodel[[gen(0)]], simParam = simParam)
	# getAcc(RGSC[[gen(0)]])
	# for(i in 1:4){
	for(i in 1:(nYr + nTrial)){
		if(i <= nYr){
			if(verbose) cat("Year:", i, "\n")
			# i = 1
			# predict latest RGSC with updated GS model 
			if(i > 1) {
				# RGSC[[gen(GScylce[1]-1)]] <- setEBV(RGSC[[gen(GScylce[1]-1)]], GSmodel[[gen(i-1)]], simParam = simParam)
				RGSC[[gen(GScylce[1]-1)]] <- setEBV(RGSC[[gen(GScylce[1]-1)]], GSmodel[[length(GSmodel)]], simParam = simParam)
			 }

			# select on mean or expected quantile?
			seltr <- if(selDHquantile) simDHdist(RGSC[[length(RGSC)]], GSmodel[[length(RGSC)]]) else 1
			# make selections for DH parents
			selGStoP <- selectInd(RGSC[[length(RGSC)]], nInd = nDHfam, trait = 1, use = "ebv") # does this select from specific families? Almost certainly.
			# selGStoP@id
			

			# make DH families
			VDP[[trials[1]]][[gen(i)]] <- makeDH(selGStoP, nDH = DHfamSize)
			# print mean genotypic value of DH 
			if(verbose) print(sapply(VDP[[trials[1]]], function(x) mean(gv(x))))
		}

		# get generation indicies
		genI <- tail(1:i, min(5, i))
		genBack <- abs(genI - i) + 1
		index = 1:length(genI)
		for(g in index) {
			gi <- genI[g]
			gb <- genBack[g]
			ti <- trials[gb]
			#phenotype if not skipped
			if(!ti %in% skip) VDP[[ti]][[gen(gi)]] <- setPheno(VDP[[ti]][[gen(gi)]], varE = h2toVe(h2[gb], Vg), reps = trialReps[gb] * trialLocs[gb]) 	
			# set ebv (does this use phenotypes if not set above? need to check...)
			if(select == "ebv" | !is.null(skip)) VDP[[ti]][[gen(gi)]] <- setEBV(VDP[[ti]][[gen(gi)]], GSmodel[[gen(i-1)]], simParam = simParam)
			sel <- if(ti %in% skip) "ebv" else  select
			#select indviduals for next years trial based on ebv and/or phenotype
			if(i - gi < nTrial) VDP[[trials[gb + 1]]][[gen(gi)]] <- selectInd(VDP[[ti]][[gen(gi)]], nInd = selectTrials[gb], trait = 1, use = sel, returnPop = TRUE)
		}

		if(i <= nYr){
			# run GS model to cycle through RGSC for year i
			for(j in GScylce){
				if(j != GScylce[1]) RGSC[[gen(j - 1)]] <- setEBV(RGSC[[gen(j-1)]], GSmodel[[gen(i - 1)]], simParam = simParam)
				RGSC[[gen(j)]] <- selectCross(pop = RGSC[[gen(j-1)]], nInd = RGSC[[gen(j-1)]]@nInd * RGSCintensity, 
											   use = RGSCuse,  trait = 1, simParam = simParam, nCrosses = nNuclear, nProgeny = 1) 
			}
			GScylce <- GScylce + GScylcePerYr

			# return lines from VDP into the RGSC 
			if(!is.null(returnVDPtoRGSC)){
				returnToRGSC <- trials[genBack] %in% returnVDPtoRGSC
				if(sum(returnToRGSC) > 0){	
					addToRGSC <- list()
					for(g in index[returnToRGSC]) {
						gi <- genI[g]
						ti <- trials[genBack[g]]
						addToRGSC[[gen(g)]] <- VDP[[ti]][[gen(gi)]]
					}
					if(length(addToRGSC) > 0){
						addToRGSC <- Reduce(c, addToRGSC)
						RGSC[[gen(GScylce[1] - 1)]] <- c(RGSC[[gen(GScylce[1] - 1)]], addToRGSC)
					}
				}
			}

			# add new phenotypes to training set and retrain GS model
			trnSet <- lapply(VDP[trials[!grepl("variety", trials)]], function(x) x[names(x) %in% gen((i-max(1, lgen)):i)])
	 		hasPop <- sapply(trnSet, length) > 0
			train <- Reduce(c, lapply(trnSet[hasPop], function(x) Reduce(c, x)))
			cat("training set has ", train@nInd, "individuals...\n")	
			GSmodel[[gen(i)]] <- RRBLUP(train, traits = 1, use = "pheno", snpChip = 1, simParam=simParam)
		} else {
			cat("final year reached, selecting on phenotypes / ebv trainde with last year training set ...\n")
			GSmodel[[gen(i)]] <- GSmodel[[gen(i-1)]]
		}
	}
	rL <- returnFunc(list(RGSC = RGSC, VDP = VDP, GSmodel = GSmodel))
	return(rL)
}

# simDHdist <- function(pop, returnQuantile = 0.9){

# }

