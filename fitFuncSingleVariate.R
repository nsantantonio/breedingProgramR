
# save x quantile value!!!!!!!!!!!! check
# burn in RGSC random crosses CHECK
# use mu + x sigma inside RGSC too! --> mu or mu + x sigma or x sigma # check?
# add truely traditional VDP -- how to do? CHECK

# x drops as VDP gets small -- CDF of normal distribution. 
# x as function of heritability and VDP size/intensity
# VDP size relative to the heritability for fixed x
# selection efficiency of VDP

# add sparse selection



sim <- function(k = 1, founderPop, paramL, simParam = SP, returnFunc = identity, verbose = TRUE, checkParam = FALSE, GSfunc = NULL, switchGSfunc = 4, ...){
# k = 1; paramL = defArgs; simParam <- SP; acrossInt = 0.5; withinInt = 0.2; verbose = TRUE; checkParam = FALSE; GSfunc = RRBLUP; nGenOut = NULL; nGenInbr = NULL; returnFunc = getPopStats
	# parameter checks and warnings.
	if (checkParam){
		paramNames <- c("maxIter", "lgen", "useTrue", "traditional", "founderBurnIn", "selectRGSC", "nProgenyPerCrossIn", "nProgenyPerCrossOut", 
						"useIn", "selectOut", "useVDP", "returnVDPcrit", "selFuncOut", "selFuncIn", "withinFamInt", "setXint", "skip", 
						"nFounder", "nNuclear", "nFam", "famSize", "ssd", "selF2", "nF2", "Vg", "updateVg", "h2", "nYr", "selectTrials", "trialReps", 
						"trialLocs", "cyclePerYr", "returnVDPtoRGSC", "nChrom", "nLoci", "nM", "nQTL")
		if (!all(paramNames %in% names(paramL))) stop("not all parameters in 'paramL'! Please include all these parameters in parameter list:\n", paste0(paramNames, "\n"))
	}
	
	for (p in names(paramL)) assign(p, paramL[[p]])
	if(traditional > 0) msg(1, "Intensity across families:", intAcross, "Intensity within family:", intWithin)
	
	selFuncStop <- c(" is a list, but of the wrong structure. Please provide either a single function, a list of functions for for each year, each cycle, or a nested list of length 'nYr' with each element of length 'cyclePerYr'")
	if(phenoRGSC < 0 | phenoRGSC > 1) stop("'phenoRGSC' must be between 0 and 1 representing the proportion of the first VDP trial dedicated to phenotyping the RGSC!")

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

	if (withinFamInt > 1) withinFamInt <- famSize / ((nI + withinFamInt) / nFam)

	# count and rename trials
	nTrial <- length(selectTrials)
	trials <- c(paste0("trial", 1:nTrial), "variety")
	if (!is.null(skip)) skip <- trials[skip]

	if(traditional > 0) {
		# useIn <- "pheno" # can still use markers in traditional program...
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
	if(is.null(pullCycle)) pullCycle <- cyclePerYr
	cycle <- 1:cyclePerYr

	# ignore GS models if only phenotypes used. 
	noGS <- if(useIn == "pheno" & useOut == "pheno" & withinFamInt == 1) TRUE else FALSE 
	if(noGS) msg(1, "noGS is true!")
	# initialize lists to store populations
	RGSC <- list()
	if(phenoRGSC > 0){
		RGSCtoVDP <- list()
		RGSCtoVDPmodel <- list()
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
	RGSC[[gen(0)]] <- newPop(founderPop)
	if(verbose) msg(1, "Founder population has genetic variance of:", popVar(getTrueBV(RGSC[[gen(0)]], simParam = simParam)))
	
	# burn in using random crosses
	printBurnin <- TRUE
	while (nFam > nInd(RGSC[[gen(0)]]) | founderBurnIn > 0) {
		if(nFam > nInd(RGSC[[gen(0)]]) & verbose) msg(1, "nFounder < nFam. Random mating to make nFam parents...")
		if(printBurnin & founderBurnIn > 0 & verbose) msg(1, "Running", founderBurnIn, "burn-in cycle(s) of random mating...")
		RGSC[[gen(0)]] <- selectCross(RGSC[[gen(0)]], nInd = nInd(RGSC[[gen(0)]]), use = "rand", simParam = simParam, nCrosses = max(nFam, nInd(RGSC[[gen(0)]]))) 
		founderBurnIn <- founderBurnIn - 1
		if(printBurnin) printBurnin <- FALSE
	}
	# phenotype founders
	RGSC[[gen(0)]] <- setPheno(RGSC[[gen(0)]], varE = h2toVe(founderh2, Vg), reps = founderReps)

	if(!noGS){	
		# train GS model and set EBV (after selfing if ssd)
		if (selF2) RGSC[[gen(0)]] <- self(RGSC[[gen(0)]], nProgeny = nF2, simParam = simParam)
		GSmodel[[gen(0)]] <- RRBLUP(RGSC[[gen(0)]], traits = 1, use = gen0use, snpChip = 1, simParam = simParam, useQtl = useTrue)
		# GSmodel[[gen(0)]] <- do.call(GSfunc, getArgs(GSfunc, pop = RGSC[[gen(0)]], traits = 1, use = gen0use, snpChip = 1, simParam = simParam, useQtl = useTrue, maxIter = 1000L, ...))
		RGSC[[gen(0)]] <- setEBV(RGSC[[gen(0)]], GSmodel[[gen(0)]], simParam = simParam)
		predAcc[["RGSC"]][[gen(0)]] <- getAcc(pop = RGSC[[gen(0)]], simParam = simParam)
		
		# initSNPAlFreq <- attributes(scale(pullSnpGeno(RGSC[[gen(0)]])))[["scaled:center"]] / 2
		# initQTLAlFreq <- attributes(scale(pullQtlGeno(RGSC[[gen(0)]])))[["scaled:center"]] / 2
	}
	
	xInt <- if(is.null(setXint)) 1 - selectTrials[nTrial] / selectTrials[1] else setXint 
	if(verbose) msg(0, "Estimated selection intensity:", round(qnorm(xInt, sd = sqrt(Vg)), 3))

	# txtdensity(initSNPAlFreq)
	# pop <- RGSC[[gen(0)]]; GSfit <- GSmodel[[gen(0)]]
	# pop <- RGSC[[lastRGSCgen]]; GSfit <- GSmodel[[lastGSmodel]]

	pullRGSCgen <- gen(0)
	pullGSmodel <- gen(0)

	if(!exists("nElite")) nElite <- nFam 
	if(traditional > 0) elite[[gen(0)]] <- selectInd(RGSC[[gen(0)]], nInd = min(nInd(RGSC[[gen(0)]]), nElite), use = useIn)
	# sapply(RGSC, function(x){mean(gv(x))})
	# rlapply(VDP, function(x){mean(gv(x))}, level = 2, combine = c)
	# run program for nYr years
	for (i in 1:(nYr + nTrial - 1)) { 
	# for (i in 1:5) { 
		# i = 1
		lastRGSCgen <- names(RGSC)[length(RGSC)]
		lastGSmodel <- if (i <= nYr) gen(i-1) else gen(nYr)
		selectRGSCi <- if(deltaRGSC) selectRGSC[i] else selectRGSC[1]
		# selectTrialsi <- selectTrials
		# nProgenyPerCrossIni <- if(identical(expDistPairs, selFuncIn)) nNuclear / selectRGSCi  * nProgenyPerCrossIn else nProgenyPerCrossIn

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

			# predict latest RGSC with updated GS model 
			if (!noGS & i > 1) {
				RGSC[[lastRGSCgen]] <- setEBV(RGSC[[lastRGSCgen]], GSmodel[[lastGSmodel]], simParam = simParam)
				# if (pullCycle < cyclePerYr) {
					# pullGSmodel <- gen(i - 1) # note, this should actually be defaulting to lastGSmodel
				# 	RGSC[[pullRGSCgen]] <- if(separateTrain) setEBV(RGSC[[pullRGSCgen]], RGSCtoVDPmodel[[pullGSmodel]], simParam = simParam) else setEBV(RGSC[[pullRGSCgen]], GSmodel[[pullGSmodel]], simParam = simParam)
				# }
			}

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
				# msg(1, "Correlation of genetic value and phennotype of selected individuals", cor(gv(selToP), pheno(selToP)))
				# mean(gv(selToP))
				# mean(pheno(selToP))

				selToP <- if(nFam < nInd(RGSC[[lastRGSCgen]])) selectInd(RGSC[[lastRGSCgen]], nInd = nFam, trait = 1, use = useOut) else RGSC[[lastRGSCgen]]
			} else {
				pullGSmodel <- gen(i - 1) # note, this should actually be defaulting to lastGSmodel
				# select out of RGSC
				# cyclePerYr - pullCycle + delay * cyclePerYr
				if(is.null(nGenOut)) nGenOut <- cyclePerYr - pullCycle # if longer you need to count!
				# GSmodelOut <- if(separateTrain) RGSCtoVDPmodel[c(pullGSmodel, lastGSmodel)] else GSmodel[c(pullGSmodel, lastGSmodel)]
				GSmodelOut <- if(separateTrain) RGSCtoVDPmodel[[pullGSmodel]] else GSmodel[[lastGSmodel]]
				selToP <- do.call(selFuncOut[[i]], getArgs(selFuncOut[[i]], pop = RGSC[[pullRGSCgen]], GSfit = GSmodelOut, nSel = nFam, nGenOut = nGenOut, nGenThisYr = cyclePerYr - pullCycle, 
												  trait = 1, use = useOut, quant = xInt, nProgeny = nProgenyPerCrossOut, Gvar = Gvar, simParam = simParam, ...))
												  # trait = 1, use = useOut, quant = xInt, nProgeny = nProgenyPerCrossOut, Gvar = Gvar, simParam = simParam, fthreshOut = 0.2))
				# double check this is the correct GS model!!!! Maybe needs to be pullGSmodel???
			}
			# check GP accuracy for material to VDP
			if(!noGS) {
				selToP <- setEBV(selToP, GSmodel[[lastGSmodel]], simParam = simParam) #
				predAcc[["RGSCout"]][[gen(i)]] <- getAcc(pop = selToP, simParam = simParam)
			}
			# if(useOut == "pheno"){

			# }
			
			if(verbose) msg(1, "VDP input Vg:", round(varA(selToP), 6))
			if(verbose) msg(1, "VDP input pop mean:", round(mean(gv(selToP)), 6))
			famSizei <- round(nFam / nInd(selToP) * famSize) 

			# Save some proportion of RGSC to put into VDP to update training models. 
			if(phenoRGSC > 0) {
				RGSCtoVDPgen <- lastRGSCgen
				# see how many lines from RGSC you can sample, then add back to famSizei
				newFamSizei <- round(famSizei * (1 - phenoRGSC))
				availPlots <- (nInd(selToP) * famSizei) - (nInd(selToP) * newFamSizei)
				RGSCfamSize <- availPlots / nInd(RGSC[[RGSCtoVDPgen]])
				if(RGSCfamSize < 1) {
					nRGSC <- RGSCfamSize * nInd(RGSC[[RGSCtoVDPgen]]) 
					RGSCfamSize <- 1
					whichRGSCtoVDP <- sample(nInd(RGSC[[RGSCtoVDPgen]]), nRGSC)
				} else {
					nRGSC <- nInd(RGSC[[RGSCtoVDPgen]]) 
					RGSCfamSize <- floor(RGSCfamSize)
					whichRGSCtoVDP <- 1:nInd(RGSC[[RGSCtoVDPgen]])
				}
				availPlots <- availPlots - RGSCfamSize * nRGSC
				famSizei <- newFamSizei + floor(availPlots / nInd(selToP))
				
				if(ssd){
					RGSCtoVDP[[gen(i)]] <- self(RGSC[[RGSCtoVDPgen]][whichRGSCtoVDP], nProgeny = RGSCfamSize)
					for(i in 1:(cyclePerYr-1)) RGSCtoVDP <- self(RGSC[[RGSCtoVDPgen]], nProgeny = 1)
				} else {
					RGSCtoVDP[[gen(i)]] <- makeDH(RGSC[[RGSCtoVDPgen]][whichRGSCtoVDP], nDH = RGSCfamSize)
				}
				if(verbose) {
					msg(1, nInd(RGSCtoVDP[[gen(i)]]), "inbred lines from", nRGSC, "RGSC included in VDP for making crosses")
					msg(1, "Family size reduced from", famSize, "to", famSizei)
				}
				# selectTrialsi[1] <-  selectTrialsi[1] # dont update, allows for no selection in first
			}

			if(traditional > 0 & i > 1) {
				# if(verbose) msg(1, nInd(selToP), "lines selected out of VDP", VDPsel, "for making crosses")
				if(verbose) msg(1, nInd(selToP), "crosses made out of VDP", VDPsel)
			} else if(i > 1) {
				if(verbose) msg(1, nInd(selToP), "crosses selected out of", nInd(RGSC[[pullRGSCgen]]), "RGSC population for VDP")
			} else {
				if(verbose) msg(1, nInd(selToP), "crosses selected out of", nInd(RGSC[[gen(0)]]), "founder population for VDP")
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
			if(phenoRGSC > 0 & i > 1 & g == 1 & i <= nYr) RGSCtoVDP[[gen(i)]] <- setPheno(RGSCtoVDP[[gen(i)]], varE = h2toVe(h2[genBack[g]], Vgi), reps = trialReps[genBack[g]] * trialLocs[genBack[g]])
		}

		if (i <= nYr){
			if(SP$isTrackPed) ped[[gen(i)]] <- cbind(VDP[[trials[1]]][[gen(i)]]@mother, VDP[[trials[1]]][[gen(i)]]@father, VDP[[trials[1]]][[gen(i)]]@id)
			# if(traditional) 
			for (j in cycle){
				jp <- which(cycle == j)
				if(traditional > 0) {
					# Ok, this is a problem. I need to add 'returnToRGSC' here. I htink I need to restructure how the RGSC works for a traditional program. It just needs to hold cross candidates. 
					whichTrials <- which(sapply(VDP, length) > 0)
					VDPsel <- if(is.logical(traditional) | max(whichTrials) < traditional) tail(names(VDP)[whichTrials], 1) else paste0("trial", traditional)
					oriGen <- i - as.numeric(gsub("[A-z]", "", VDPsel)) + 1 # this is hacky, but it works...
					
					# VDPsel <- paste0("trial", min(i, traditional))
					# oriGen <- i - min(i, traditional) + 1
					selPop <- VDP[[VDPsel]][[oriGen]] 
					families <- split(VDP[[trials[1]]][[gen(oriGen)]]@id, rep(1:nFam, each = famSizei)) # this assumes famSizei is a single integer, need to allow for diff number ind per fam!!!!!!!!!!!!!!!
					if(verbose) msg(1, "Individuals selected out of", VDPsel, "from generation", oriGen, "for crossing...")
				
				} else {
					RGSC[[gen(j-1)]] <- setEBV(RGSC[[gen(j-1)]], GSmodel[[lastGSmodel]], simParam = simParam)
					# if (j != cycle[1]) RGSC[[gen(j-1)]] <- setEBV(RGSC[[gen(j-1)]], GSmodel[[lastGSmodel]], simParam = simParam)
					selPop <- RGSC[[gen(j-1)]]
					predAcc[["RGSC"]][[gen(j-1)]] <- getAcc(selPop)
				}
				# run GS model to cycle through RGSC for year i
 				if(traditional > 0) {
					RGSC[[gen(j)]] <- do.call(tradSelCross2, getArgs(tradSelCross2, pop = selPop, elite = elite[[gen(i)]], families = families, nFam = nFam, famSize = famSizei, use = useIn, trait = 1, simParam = simParam, 
						nCrosses = nFam, nProgeny = nProgenyPerCrossIn, verbose = verbose, ...))
						# nCrosses = nFam, nProgeny = nProgenyPerCrossIn, verbose = verbose))
				} else if(is.null(selFuncIn)){
					RGSC[[gen(j)]] <- selectCross(pop = selPop, nInd = min(selectRGSCi, nInd(selPop)), use = useIn,  trait = 1, simParam = simParam, nCrosses = nNuclear, nProgeny = nProgenyPerCrossIn)
				} else {
					RGSC[[gen(j)]] <- do.call(selFuncIn[[i]][[jp]], getArgs(selFuncIn[[i]][[jp]], nSel = selectRGSCi, pop = selPop, GSfit = GSmodel[[lastGSmodel]],
					trait = 1,  use = useIn, nCrosses = nNuclear, nProgeny = nProgenyPerCrossIn, quant = xInt, verbose = verbose, Gvar = Gvar, simParam = simParam, ...))
					# trait = 1,  use = useIn, nCrosses = nNuclear, nProgeny = nProgenyPerCrossIn, quant = xInt, verbose = verbose, Gvar = Gvar, simParam = simParam, fthresh = 0.01))
					if(nInd(RGSC[[gen(j)]]) != nNuclear) msg(2, "Only", nInd(RGSC[[gen(j)]]), "crosses made in RGSC...")
				}
				# would be good to be able to select within f2 family if f2 > 1
				if (selF2) RGSC[[gen(j)]] <- self(RGSC[[gen(j)]], nProgeny = nF2, simParam = simParam)
				# if(jp == pullCycle) {
				# 	pullRGSCgen <- gen(j) # should I push material out earlier?
				# }
				if(verbose) msg(1, "RGSC Vg:", round(varA(RGSC[[gen(j)]]), 6))
				if(verbose) msg(1, "RGSC Pop Mean:", round(mean(RGSC[[gen(j)]]@gv), 6))
			}
			pullRGSCgen <- if(pullCycle > 0) gen(cycle[[pullCycle]]) else lastRGSCgen
			# update GScycle number
			cycle <- cycle + cyclePerYr
			
			if(!noGS){
				# so this removes the selections from the training population. 
				trnSet <- lapply(VDP[trials[!grepl("variety", trials)]], function(x) x[names(x) %in% gen(max(1, i-max(1, lgen)):i)])
				trnSet <- trnSet[sapply(trnSet, length) > 0]
				names(trnSet)
				# include founders in prediction model
				if(founderKeep >= i) trnSet[["founder"]] <- RGSC[[gen(0)]]

				# Inlcude RGSC individuals into training
				if(phenoRGSC > 0) {
					if(separateTrain) {
						RGSCtoVDPtrian <- RGSCtoVDP
						# RGSCtoVDPtrian <- RGSCtoVDP[names(RGSCtoVDP) %in% gen(max(1, i-max(1, lgen)):i)]  
						if(founderKeep >= i) RGSCtoVDPtrian["founder"] <- RGSC[gen(0)]
						RGSCtoVDPtrian <- mergePopsRec(RGSCtoVDPtrian) 
						msg(1, "Training set for RGSCout branch has", RGSCtoVDPtrian@nInd, "individuals...")	
						if(!is.null(GSfunc)){
							RGSCtoVDPmodel[[gen(i)]] <- do.call(GSfunc, getArgs(GSfunc, pop = RGSCtoVDPtrian, traits = 1, use = useGS, snpChip = 1, simParam=simParam, useQtl = useTrue, maxIter = maxIter))#, ...))
						} else {
							if(nInd(RGSCtoVDPtrian) > simParam$snpChips[[1]]@nLoci * switchGSfunc) {
								RGSCtoVDPmodel[[gen(i)]] <- do.call(RRBLUP2, getArgs(RRBLUP2, pop = RGSCtoVDPtrian, traits = 1, use = useGS, snpChip = 1, simParam=simParam, useQtl = useTrue, maxIter = maxIter))#, ...))
							} else {
								RGSCtoVDPmodel[[gen(i)]] <- do.call(RRBLUP, getArgs(RRBLUP, pop = RGSCtoVDPtrian, traits = 1, use = useGS, snpChip = 1, simParam=simParam, useQtl = useTrue, maxIter = maxIter))#, ...))
							}
						}
					} else {
						trnSet[["RGSCtoVDP"]] <- RGSCtoVDP[names(RGSCtoVDP) %in% gen(max(1, i-max(1, lgen)):i)]
					}
				}
				
				# concatenate training set and train GS model
		 		train <- mergePopsRec(trnSet) 
				msg(1, "Training set has", train@nInd, "individuals...")	
				# GSmodel[[gen(i)]] <- GSfunc(train, traits = 1, use = useGS, snpChip = 1, simParam=simParam, useQtl = useTrue)
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
		
		# Below needs to be checked for consistancy! 
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

	# sapply(RGSC, function(x){mean(gv(x))})
	# rlapply(VDP, function(x){mean(gv(x))}, level = 2, combine = c)
	# txtplot(1:{i - nTrial + 1}, sapply(VDP[["variety"]], function(x){mean(gv(x))}))
	# txtplot(1:i, sapply(VDP[["trial1"]], function(x){mean(gv(x))}))

	# txtplot(1:nYr, sapply(VDP[["trial1"]], function(x){mean(gv(x))}))
	# txtplot(1:nYr, sapply(VDP[["variety"]], function(x){mean(gv(x))}))

	# M <- do.call(rbind, lapply(VDP[["variety"]], pullSnpGeno))
	# rownames(M) <- paste0(rep(names(VDP[["variety"]]), times = sapply(VDP[["variety"]], nInd)), "_", rownames(M))
	# K <- genCov(M)
	# pdf("figures/multiGenProblem.pdf")
	# image(K[, ncol(K):1])
	# dev.off()
	
	rL <- list(SP = SP, paramL = paramL, RGSC = RGSC, VDP = VDP, GSmodel = GSmodel, predAcc = predAcc, selIntensity = intensity)
	do.call(returnFunc, getArgs(returnFunc, resultL = rL, ...))
}
