
# save x quantile value!!!!!!!!!!!! check
# burn in RGSC random crosses CHECK
# use mu + x sigma inside RGSC too! --> mu or mu + x sigma or x sigma # check?
# add truely traditional VDP -- how to do? CHECK

# x drops as VDP gets small -- CDF of normal distribution. 
# x as function of heritability and VDP size/intensity
# VDP size relative to the heritability for fixed x
# selection efficiency of VDP

# add sparse selection



# k = 2; paramL = defArgs; simParam <- SP; verbose = TRUE; checkParam = FALSE; GSfunc = RRBLUP; nGenOut = NULL; nGenInbr = NULL
sim <- function(k = 1, founderPop, paramL, simParam = SP, returnFunc = identity, verbose = TRUE, checkParam = FALSE, GSfunc = NULL, switchGSfunc = 4, ...){

	# parameter checks and warnings.
	if (checkParam){
		paramNames <- c("maxIter", "lgen", "useTrue", "traditional", "founderBurnIn", "selectRGSC", "nProgenyPerCrossIn", "nProgenyPerCrossOut", 
						"useIn", "selectOut", "useVDP", "returnVDPcrit", "selFuncOut", "selFuncIn", "withinFamInt", "setXint", "skip", 
						"nFounder", "nNuclear", "nFam", "famSize", "ssd", "selF2", "nF2", "Vg", "updateVg", "h2", "nYr", "selectTrials", "trialReps", 
						"trialLocs", "cyclePerYr", "returnVDPtoRGSC", "nChrom", "nLoci", "nM", "nQTL")
		if (!all(paramNames %in% names(paramL))) stop("not all parameters in 'paramL'! Please include all these parameters in parameter list:\n", paste0(paramNames, "\n"))
	}
	
	for (p in names(paramL)) assign(p, paramL[[p]])
	
	selFuncStop <- c(" is a list, but of the wrong structure. Please provide either a single function, a list of functions for for each year, each cycle, or a nested list of length 'nYr' with each element of length 'cyclePerYr'")
	
	# funcs <- list(selFuncIn = selFuncIn, selFuncOut = selFuncOut, inbreedFunc = inbreedFunc)
	# for(f in names(funcs)){
	# 	if(is.list(funcs[[f]])) {
	# 		if (length(funcs[[f]]) %in% c(cyclePerYr, nYr, cyclePerYr * nYr)) {
	# 			if(all(sapply(funcs[[f]], class) == "list")){
	# 				if (!all(unique(sapply(funcs[[f]], length)) == cyclePerYr)) stop(paste0(f, funcListStop))
	# 			}
	# 			if(length(funcs[[f]]) == cyclePerYr & f == "selFuncIn") funcs[[f]] <- rep(list(funcs[[f]]), nYr)
	# 			if(length(funcs[[f]]) == nYr & f == "selFuncIn") funcs[[f]] <- lapply(1:nYr, function(x) rep(list(funcs[[f]][[x]]), cyclePerYr))
	# 		} else {
	# 			stop(paste0(f, funcListStop))
	# 		}
	# 	} else {
	# 		if(!is.null(funcs[[f]])) {
	# 			if(class(funcs[[f]]) == "function") funcs[[f]] <- rep(list(rep(list(funcs[[f]]), cyclePerYr)), nYr) else stop(paste0(f, funcListStop))
	# 		}
	# 	}
	# }
	# for (f in names(funcs)) assign(f, funcs[[f]])

	if(is.null(pullCycle)) pullCycle <- cyclePerYr

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

	# if(any(sapply(unlist(selFuncIn), function(x) identical(x, solqp)))) require(LowRankQP)

	if(!(length(selectRGSC) == 1 | length(selectRGSC) == nYr)) stop("'selectRGSC' must be of length 1 or nYr!")
	if(all(selectRGSC <= 1)) {
		selectRGSC <- nNuclear * selectRGSC
		if(selectRGSC %% 1 != 0) selectRGSC <- round(selectRGSC)
		deltaRGSC <- if(length(selectRGSC) == nYr) TRUE else FALSE
	}

	if(nFounder < nInd(founderPop)) {
		subFounder <-  if (is.null(founderSamples)) sample(1:nInd(founderPop), nFounder) else founderSamples[[k]]
		founderPop <- founderPop[subFounder]
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
		msg(0, "NOTE: Selection intensities have been rounded to the nearest integer:\n", selectTrials, "\nThese correspond to selection intensities of:\n", actInt)
	}
	if(length(returnVDPtoRGSC) != length(selectTrials) + 1) stop("returnVDPtoRGSC must be of length(selectTrials) + 1 (for variety)!")
	if (all(returnVDPtoRGSC <= 1)) returnVDPtoRGSC <- returnVDPtoRGSC * c(nI, selectTrials) 

	# count and rename trials
	nTrial <- length(selectTrials)
	trials <- c(paste0("trial", 1:nTrial), "variety")
	if (!is.null(skip)) skip <- trials[skip]

	if(traditional) {
		useIn <- "pheno"
		# useOut <- "pheno" # can still use markers in traditional program...
		cyclePerYr <- 1
	}
	# if use true values
	if(useTruth > 0) {
		msg(0, "NOTICE: using QTL as markers.")
		useTrue <- TRUE
	}
	if(useTruth > 1){
		# get true marker effects?
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
	if(verbose) msg(1, "Founder population has genetic variance of:", popVar(getTrueBV(RGSC[[gen(0)]], simParam = simParam)))
	
	# burn in using random crosses
	printBurnin <- TRUE
	while (nFam > nInd(RGSC[[gen(0)]]) | founderBurnIn > 0) {
		if(nFam > nInd(RGSC[[gen(0)]]) & verbose) cat("nFounder < nFam. Random mating to make nFam parents...")
		if(printBurnin & founderBurnIn > 0 & verbose) cat("Running", founderBurnIn, "burn-in cycles of random mating...")
		RGSC[[gen(0)]] <- selectCross(RGSC[[gen(0)]], nInd = nInd(RGSC[[gen(0)]]), use = "rand", simParam = simParam, nCrosses = max(nFam, nInd(RGSC[[gen(0)]]))) 
		founderBurnIn <- founderBurnIn - 1
		if(printBurnin) printBurnin <- FALSE
	}
	# train GS model and set EBV (after selfing if ssd)
	GSmodel[[gen(0)]] <- RRBLUP(RGSC[[gen(0)]], traits = 1, use = gen0use, snpChip = 1, simParam = simParam, useQTL = useTrue)
	# GSmodel[[gen(0)]] <- do.call(GSfunc, getArgs(GSfunc, pop = RGSC[[gen(0)]], traits = 1, use = gen0use, snpChip = 1, simParam = simParam, useQTL = useTrue, maxIter = 1000L, ...))
	if (selF2) RGSC[[gen(0)]] <- self(RGSC[[gen(0)]], nProgeny = nF2, simParam = simParam)
	RGSC[[gen(0)]] <- setEBV(RGSC[[gen(0)]], GSmodel[[gen(0)]], simParam = simParam)
	
	xInt <- if(is.null(setXint)) 1 - selectTrials[nTrial] / selectTrials[1] else setXint 
	if(verbose) msg(0, "Estimated selection intensity:", round(qnorm(xInt, sd = sqrt(Vg)), 3))

	# pop <- RGSC[[gen(0)]]; GSfit <- GSmodel[[gen(0)]]
	# pop <- RGSC[[lastRGSCgen]]; GSfit <- GSmodel[[lastGSmodel]]

	pullRGSCgen <- gen(0)
	pullGSmodel <- gen(0)
	
	# run program for nYr years
	for (i in 1:(nYr + nTrial - 1)) { 
	# for (i in 1:6) { 
		# i = 1
		lastRGSCgen <- names(RGSC)[length(RGSC)]
		lastGSmodel <- if (i <= nYr) gen(i-1) else gen(nYr)
		selectRGSCi <- if(deltaRGSC) selectRGSC[i] else selectRGSC[1]
		# nProgenyPerCrossIni <- if(identical(expDistPairs, selFuncIn)) nNuclear / selectRGSCi  * nProgenyPerCrossIn else nProgenyPerCrossIn

		if (i <= nYr){ 
			if (verbose) msg(0, "Year: ", i)

			# predict latest RGSC with updated GS model 
			if (i > 1) {
				RGSC[[lastRGSCgen]] <- setEBV(RGSC[[lastRGSCgen]], GSmodel[[lastGSmodel]], simParam = simParam)
				if (pullCycle < cyclePerYr) {
					pullGSmodel <- gen(i - 2) # note, this is defaulting to lastGSmodel, which is ok for this design. 
					RGSC[[pullRGSCgen]] <- setEBV(RGSC[[pullRGSCgen]], GSmodel[[pullGSmodel]], simParam = simParam)
				}
			}

			# estimate VDP selection intensity from previous generations
			if (i > nTrial) {
				intensity[[gen(i)]] <- estIntensity(VDP, i, nT = nTrial, start = "trial1", end = "variety", estFunc = estIntFunc, Gvar = Gvar)
				if(is.null(setXint)) xInt <- pnorm(intensity[[gen(i)]]) 
				if(verbose) msg(1, "Realized selection intensity:", round(intensity[[gen(i)]], 3))
				if(verbose) msg(1, "best variety:", round(max(gv(VDP[["variety"]][[gen(i - nTrial)]])), 3))
			}
			
			# use selected lines from last trial to make new crosses. Not sure this timeline is realistic...
			if(i == 1 | is.null(selFuncOut)){
				selToP <- selectInd(RGSC[[lastRGSCgen]], nInd = nFam, trait = 1, use = useOut)
			} else {
				RGSC[[pullRGSCgen]]@ebv
			# select out of RGSC, on mean, expected quantile, etc...
				# selToP <- do.call(selFuncOut, getArgs(selFuncOut, nSel = nFam, pop = RGSC[[lastRGSCgen]], GSfit = GSmodel[[lastGSmodel]], 
								  # trait = 1, use = selectOut, quant = xInt, nProgeny = nProgenyPerCrossOut, Gvar = Gvar))#, ...))
				# cyclePerYr - pullCycle 
				# cyclePerYr - pullCycle + delay * cyclePerYr
				if(is.null(nGenOut)) nGenOut <- cyclePerYr - pullCycle # if longer you need to count!
				selToP <- do.call(selFuncOut[[i]], getArgs(selFuncOut[[i]], pop = RGSC[[pullRGSCgen]], GSfit = GSmodel[c(pullGSmodel, lastGSmodel)], nSel = nFam, nGenOut = nGenOut, nGenThisYr = cyclePerYr - pullCycle, 
												  # trait = 1, use = useOut, quant = xInt, nProgeny = nProgenyPerCrossOut, Gvar = Gvar, verbose = verbose, ...))
												  # trait = 1, use = useOut, quant = xInt, nProgeny = nProgenyPerCrossOut, Gvar = Gvar, fthreshOut = 0.1))#, ...))
												  trait = 1, use = useOut, quant = xInt, nProgeny = nProgenyPerCrossOut, Gvar = Gvar, fthreshOut = 0.2))#, ...))
												  # trait = 1, use = useOut, quant = xInt, nProgeny = nProgenyPerCrossOut, Gvar = Gvar, fthreshOut = list(0.1, 0.2)))#, ...))
												  # trait = 1, use = useOut, quant = xInt, nProgeny = nProgenyPerCrossOut, Gvar = Gvar, fthreshOut = list(0.1, NULL), truncqpOut = list(NULL, nFam)))#, ...))
												  # trait = 1, use = useOut, quant = xInt, nProgeny = nProgenyPerCrossOut, Gvar = Gvar, lambdaOut = 0))
												  # trait = 1, use = useOut, quant = xInt, nProgeny = nProgenyPerCrossOut, Gvar = Gvar, truncqpOut = nFam))#, ...))
												  # trait = 1, use = useOut, quant = xInt, nProgeny = nProgenyPerCrossOut, Gvar = Gvar, fthreshOut = 0.1))#, ...))
												  # pullGeno = pullGenoFunc, w = weight))
			}
			if(verbose) msg(1, "VDP input Vg:", round(varA(selToP), 6))
			if(verbose) msg(1, "VDP input pop mean:", round(mean(gv(selToP)), 6))
			# if(nInd(selToP) != nFam) stop("selToP is wrong...")
			famSizei <- round(nFam / nInd(selToP) * famSize) 
			
			if(traditional & i > 1) {
				# lastVDPSel <- tail(names(VDP)[sapply(VDP, length) > 0], 1)
				if(verbose) msg(1, nInd(selToP), "lines selected out of VDP", lastVDPSel, "for making crosses")
			} else {
				# THIS DOES NOT MKE SENSE TO ME!!!!
				# if(verbose) cat("        ", nInd(RGSC[[lastRGSCgen]]), "individuals selected out of", nNuclear, " RGSC population \n")
				if(verbose) msg(1, nInd(selToP), "crosses selected out of", nInd(RGSC[[pullRGSCgen]]), " RGSC population for VDP")
			}
			
			# (mean(ebv(selToP)) - mean(ebv(RGSC[[lastRGSCgen]]))) / sd(ebv(RGSC[[lastRGSCgen]]))
			# Determine number of individuals to select within family (default is all)
			# nProgPerFam <- famSizei / withinFamInt
			# if ((nProgPerFam) %% 1 != 0) {
			# 	nProgPerFam <- round(nProgPerFam)
			# 	nSelToTrial <- round(nProgPerFam * withinFamInt)
			# 	cat("\nNOTE: Selection intensities within familiy have been rounded to the nearest integer resulting in", nSelToTrial, "progeny per family selected from", nProg, "progeny per family\n")
			# }


			#INTERT function here for material to VDP

			# make DH families or self
			if(is.null(inbreedFunc)) {
				nProgPerFam <- famSizei / withinFamInt
				if ((nProgPerFam) %% 1 != 0) {
					nProgPerFam <- round(nProgPerFam)
					nSelToTrial <- round(nProgPerFam * withinFamInt)
					msg(0, "NOTICE: Selection intensities within familiy have been rounded to the nearest integer resulting in", nSelToTrial, "progeny per family selected from", nProg, "progeny per family")
				}
				VDP[[trials[1]]][[gen(i)]] <- if (ssd) self(selToP, nProgeny = nProgPerFam) else makeDH(selToP, nDH = nProgPerFam)
				#select within family if intensity < 1
				if (withinFamInt < 1) {
					VDP[[trials[1]]][[gen(i)]] <- setEBV(VDP[[trials[1]]][[gen(i)]], GSmodel[[lastGSmodel]], simParam = simParam)
					fams <- split(1:(nProgPerFam * nFam), rep(1:nFam, each = nProgPerFam))
					faml <- list()
					for (j in names(fams)){
						faml[[j]] <- selectInd(VDP[[trials[1]]][[gen(i)]][fams[[j]]], nInd = famSizei, trait = 1, use = useOut) 
					}
					VDP[[trials[1]]][[gen(i)]] <- mergePops(faml)		
				}
				# # print mean genotypic value of DH 
				# if (verbose) print(sapply(VDP[[trials[1]]], function(x) mean(gv(x))))
				# if(nInd(VDP[[trials[1]]][[gen(i)]]) != nFam * famSizei) stop("selToP is wrong...")
			} else {
				if(is.null(nGenInbr)) nGenInbr <- cyclePerYr  
				VDP[[trials[1]]][[gen(i)]] <- do.call(inbreedFunc[[i]], getArgs(inbreedFunc[[i]], pop = selToP, GSfit = GSmodel[[lastGSmodel]], # use last GS model here because selectOut will burn through the rest of the year. Need to make flexible if not. 
					# trait = 1, use = useInbreed, int = withinFamInt, ssd = ssd, nProgeny = famSizei, nGenInbr = nGenInbr, simParam = simParam, ...))
					trait = 1, use = useInbreed, int = withinFamInt, ssd = ssd, nProgeny = famSizei, nGenInbr = nGenInbr, simParam = simParam))#, ...))
					# trait = 1,  use = useIn, nCrosses = nNuclear, nProgeny = nProgenyPerCrossIn, quant = xInt, verbose = verbose, Gvar = Gvar, fthresh = 0.5))
			}

			msg(1, "best new line:", round(max(gv(VDP[[trials[1]]][[gen(i)]])), 6))

			# make DH families or self
			# VDP[[trials[1]]][[gen(i)]] <- if (ssd) self(selToP, nProgeny = nProgPerFam) else makeDH(selToP, nDH = nProgPerFam)

			# #select within family if intensity < 1
			# if (withinFamInt < 1) {
			# 	VDP[[trials[1]]][[gen(i)]] <- setEBV(VDP[[trials[1]]][[gen(i)]], GSmodel[[lastGSmodel]], simParam = simParam)
			# 	fams <- split(1:(nProgPerFam * nFam), rep(1:nFam, each = nProgPerFam))
			# 	faml <- list()
			# 	for (j in names(fams)){
			# 		faml[[j]] <- selectInd(VDP[[trials[1]]][[gen(i)]][fams[[j]]], nInd = famSizei, trait = 1, use = useOut) 
			# 	}
			# 	VDP[[trials[1]]][[gen(i)]] <- mergePops(faml)		
			# }
			# # # print mean genotypic value of DH 
			# # if (verbose) print(sapply(VDP[[trials[1]]], function(x) mean(gv(x))))
			# if(nInd(VDP[[trials[1]]][[gen(i)]]) != nFam * famSizei) stop("selToP is wrong...")
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
		}

		if (i <= nYr){
			# if(traditional) 
			for (j in cycle){
				jp <- which(cycle == j)
				if(traditional) {
					# Ok, this is a problem. I need to add 'returnToRGSC' here. I htink I need to restructure how the RGSC works for a traditional program. It just needs to hold cross candidates. 
					lastVDPSel <- tail(names(VDP)[sapply(VDP, length) > 0], 1)
					selPop <- VDP[[lastVDPSel]][[length(VDP[[lastVDPSel]])]] 

				} else {
					RGSC[[gen(j-1)]] <- setEBV(RGSC[[gen(j-1)]], GSmodel[[lastGSmodel]], simParam = simParam)
					# if (j != cycle[1]) RGSC[[gen(j-1)]] <- setEBV(RGSC[[gen(j-1)]], GSmodel[[lastGSmodel]], simParam = simParam)
					# if (j != cycle[1] & pullCycle != cyclePerYr) RGSC[[gen(j-1)]] <- setEBV(RGSC[[gen(j-1)]], GSmodel[[lastGSmodel]], simParam = simParam)
###################################
					# SOMTHING IS BROKEN HERE! 	NOT PREDICTED??????!!!!!!
					# msg(3, cor(ebv(RGSC[[gen(j-1)]]), bv(RGSC[[gen(j-1)]])))
					# predAcc[["RGSC"]][[gen(j-1)]] <- cor(ebv(RGSC[[gen(j-1)]]), bv(RGSC[[gen(j-1)]]))
					selPop <- RGSC[[gen(j-1)]]
					predAcc[["RGSC"]][[gen(j-1)]] <- getAcc(selPop)
					# print(cor(ebv(selPop), bv(selPop)))
					# slotNames(selPop)
					# bv(selPop)
					# predAcc[["RGSC"]][[gen(j-1)]] <- GPaccur(selPop)
					# print(getAcc(pop = RGSC[[gen(j-1)]], simParam = simParam))
###################################
				}
			# run GS model to cycle through RGSC for year i
				if(is.null(selFuncIn)){
					RGSC[[gen(j)]] <- selectCross(pop = selPop, nInd = min(selectRGSCi, nInd(selPop)), use = useIn,  trait = 1, simParam = simParam, nCrosses = nNuclear, nProgeny = nProgenyPerCrossIn)
				} else {
					RGSC[[gen(j)]] <- do.call(selFuncIn[[i]][[jp]], getArgs(selFuncIn[[i]][[jp]], nSel = selectRGSCi, pop = selPop, GSfit = GSmodel[[lastGSmodel]],
					# trait = 1,  use = useIn, nCrosses = nNuclear, nProgeny = nProgenyPerCrossIn, quant = xInt, verbose = verbose, Gvar = Gvar, ...))
					trait = 1,  use = useIn, nCrosses = nNuclear, nProgeny = nProgenyPerCrossIn, quant = xInt, verbose = verbose, Gvar = Gvar, fthresh = 0.01))
					if(nInd(RGSC[[gen(j)]]) != nNuclear) msg(2, "Only", nInd(RGSC[[gen(j)]]), "crosses made in RGSC...")
				}
				# would be good to be able to select within f2 family if f2 > 1
				if (selF2) RGSC[[gen(j)]] <- self(RGSC[[gen(j)]], nProgeny = nF2, simParam = simParam)
				if(jp == pullCycle) {
					pullRGSCgen <- gen(j)
				}
				if(verbose) msg(1, "RGSC Vg:", round(varA(RGSC[[gen(j)]]), 6))
				if(verbose) msg(1, "RGSC Pop Mean:", round(mean(RGSC[[gen(j)]]@gv), 6))
			}
			# update GScycle number
			cycle <- cycle + cyclePerYr
			
			# so this removes the selections from the training population. 
			trnSet <- lapply(VDP[trials[!grepl("variety", trials)]], function(x) x[names(x) %in% gen(max(1, i-max(1, lgen)):i)])
			trnSet <- trnSet[sapply(trnSet, length) > 0]

			# concatenate training set and train GS model
	 		train <- mergePopsRec(trnSet) 
			msg(1, "Training set has", train@nInd, "individuals...")	
			# GSmodel[[gen(i)]] <- GSfunc(train, traits = 1, use = useGS, snpChip = 1, simParam=simParam, useQTL = useTrue)
			if(genParam(train)$varG == 0) {
				GSmodel[[gen(i)]] <- GSmodel[[gen(i - 1)]]
			} else {
				if(!is.null(GSfunc)){
					GSmodel[[gen(i)]] <- do.call(GSfunc, getArgs(GSfunc, pop = train, traits = 1, use = useGS, snpChip = 1, simParam=simParam, useQTL = useTrue, maxIter = maxIter))#, ...))
				} else {
					if(nInd(train) > simParam$snpChips[[1]]@nLoci * switchGSfunc) {
						GSmodel[[gen(i)]] <- do.call(RRBLUP2, getArgs(RRBLUP2, pop = train, traits = 1, use = useGS, snpChip = 1, simParam=simParam, useQTL = useTrue, maxIter = maxIter))#, ...))
					} else {
						GSmodel[[gen(i)]] <- do.call(RRBLUP, getArgs(RRBLUP, pop = train, traits = 1, use = useGS, snpChip = 1, simParam=simParam, useQTL = useTrue, maxIter = maxIter))#, ...))
					}
				}
			}
		}
		# check if final year
		if (i - nYr == 1) msg(0, "Final year reached, selecting on phenotypes / ebv trained with last year training set ...")
		
		# Below needs to be checked for consistancy!!!!!!!!!!!!!!!!!!!!!!!!
		for (g in index) {
			# set ebv if using to select / skip generations
			if (useVDP == "ebv" | !is.null(skip)) {
				GSmodelVDP <- if (trials[genBack[g]] %in% skip) lastGSmodel else gen(min(i, nYr))
				VDP[[trials[genBack[g]]]][[gen(genI[g])]] <- setEBV(VDP[[trials[genBack[g]]]][[gen(genI[g])]], GSmodel[[GSmodelVDP]], simParam = simParam)
				predAcc[[trials[genBack[g]]]][[gen(genI[g])]] <- GPaccur(pop = VDP[[trials[genBack[g]]]][[gen(genI[g])]])
			}

			# select based on ebv or phenotype
			sel <- if (trials[genBack[g]] %in% skip) "ebv" else  useVDP
			if (i - genI[g] < nTrial) VDP[[trials[genBack[g] + 1]]][[gen(genI[g])]] <- selectInd(VDP[[trials[genBack[g]]]][[gen(genI[g])]], nInd = min(selectTrials[genBack[g]], nInd(VDP[[trials[genBack[g]]]][[gen(genI[g])]])), trait = 1, use = sel, returnPop = TRUE)
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
						# addToRGSC <- Reduce(c, addToRGSC) # double check this works! not sure how long that Reduce was used here...
						addToRGSC <- mergePopsRec(addToRGSC)
						RGSC[[gen(cycle[1] - 1)]] <- c(RGSC[[gen(cycle[1] - 1)]], addToRGSC)
					}
				}
			}
		}
		if(verbose) cat("\n")
	}

	rL <- list(SP = SP, paramL = paramL, RGSC = RGSC, VDP = VDP, GSmodel = GSmodel, predAcc = predAcc, selIntensity = intensity)
	do.call(returnFunc, getArgs(returnFunc, resultL = rL, ...))
}
