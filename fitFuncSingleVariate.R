# simParam <- SP; select = "ebv"; returnFunc = identity; verbose = TRUE; skip = NULL; selQuantile = TRUE; SSD = FALSE
sim <- function(founderPop, simParam = SP, select = "ebv", returnFunc = identity, verbose = TRUE, skip = NULL){
	if(selF2 & GScylcePerYr > 1) warning("Selection on F2 is being performed, and there is more than 1 GS cycle per year. You may want to reduce 'GScylcePerYr' to 1")
	if(selF2 & !selQuantile) warning("Selection on F2 is being performed, but not on expected quantiles, you probably want to set selQuantile = TRUE")
	
	dummyFunc <- function(x, retrn) retrn
	# function to return expected quantiles from sampling DH individuals
	simDHdist <- function(pop, GSfit = GSmodel[[1]], DHsampleSize = 200, returnQuantile = 0.9){
		qDH <- function(i){
			DH <- makeDH(pop[i], nDH = DHsampleSize)
			DH <- setEBV(DH, GSfit, simParam = simParam)
			DHebv <- ebv(DH)
			quantile(DHebv, returnQuantile)[[1]]
		}
		sapply(pop0@id, qDH)
	}

	# select inside RGSC?
	RGSCuse <- if(RGSCselect) "ebv" else  "rand"

	# check selectTrials
	if(!all(selectTrials > 0) | (any(selectTrials < 1) & any(selectTrials > 1))) stop("'selectTrials' must have elements between 0 and 1 or positiive integers")
	
	# define number of individuals per cycle, and number to select at each stage
	nI <- nFam * famSize
	if(all(selectTrials < 1)) selectTrials <- nI * cumprod(selectTrials)
	
	# count and rename trials
	nTrial <- length(selectTrials)
	trials <- c(paste0("trial", 1:nTrial), "variety")
	returnVDPtoRGSC <- trials[returnVDPtoRGSC]
	if(!is.null(skip)) skip <- trials[skip]

	if(!is.null(returnVDPtoRGSC)) if(!all(returnVDPtoRGSC %in% trials[-length(trials)])) stop("somthing is wrong ith returnVDPtoRGSC...") 

	# define cycles
	GScylce <- 1:GScylcePerYr

	# initialize lists to store populations
	RGSC <- list()
	GSmodel <- list()
	names(trials) <- trials
	VDP <- lapply(trials, function(x) list())

	# initialize nuclear population
	RGSC[[gen(0)]] <- newPop(founderPop)
	GSmodel[[gen(0)]] <- RRBLUP(RGSC[[gen(0)]], traits = 1, use = "pheno", snpChip = 1, simParam = simParam)
	RGSC[[gen(0)]] <- setEBV(RGSC[[gen(0)]], GSmodel[[gen(0)]], simParam = simParam)
	if(selF2) RGSC[[gen(0)]] <- self(RGSC[[gen(0)]], nProgeny = nF2, simParam = simParam)
	# getAcc(RGSC[[gen(0)]])

	# run prgram for nYr years
	for(i in 1:(nYr + nTrial)){
		if(i <= nYr){
			if(verbose) cat("Year:", i, "\n")
			# i = 1
			if(i > 1) {
				# predict latest RGSC with updated GS model 
				RGSC[[gen(GScylce[1]-1)]] <- setEBV(RGSC[[gen(GScylce[1]-1)]], GSmodel[[length(GSmodel)]], simParam = simParam)
			}
			# select on mean or expected quantile
			if(selQuantile) {
				expQuant <- simDHdist(RGSC[[length(RGSC)]], GSmodel[[length(GSmodel)]])
				selGStoP <- selectInd(RGSC[[length(RGSC)]], nInd = nFam, trait = dummyFunc, use = "ebv", retrn = expQuant) 
			} else {
				selGStoP <- selectInd(RGSC[[length(RGSC)]], nInd = nFam, trait = 1, use = "ebv") 
			}
			# make DH families
			VDP[[trials[1]]][[gen(i)]] <- if(SSD) self(selGStoP, nProgeny = famSize) else makeDH(selGStoP, nDH = famSize)
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
			if(SSD) VDP[[trials[gb + 1]]][[gen(gi)]] <- self(VDP[[trials[gb + 1]]][[gen(gi)]])
		}

		if(i <= nYr){

			# run GS model to cycle through RGSC for year i
			for(j in GScylce){
				if(j != GScylce[1]) RGSC[[gen(j - 1)]] <- setEBV(RGSC[[gen(j-1)]], GSmodel[[gen(i - 1)]], simParam = simParam)
				RGSC[[gen(j)]] <- selectCross(pop = RGSC[[gen(j-1)]], nInd = RGSC[[gen(j-1)]]@nInd * RGSCintensity, 
											   use = RGSCuse,  trait = 1, simParam = simParam, nCrosses = nNuclear, nProgeny = 1) 
				if(selF2) RGSC[[gen(j)]] <- self(RGSC[[gen(j)]], nProgeny = nF2, simParam = simParam)
			}
			# update GScycle number
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
			cat("final year reached, selecting on phenotypes / ebv trained with last year training set ...\n")
			GSmodel[[gen(i)]] <- GSmodel[[gen(i-1)]]
		}
	}
	rL <- returnFunc(list(RGSC = RGSC, VDP = VDP, GSmodel = GSmodel))
	return(rL)
}

# simDHdist <- function(pop, returnQuantile = 0.9){

# }

