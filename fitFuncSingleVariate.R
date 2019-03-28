# simParam <- SP; select = "ebv"; returnFunc = identity; verbose = TRUE; skip = NULL
sim <- function(founderPop, simParam = SP, select = "ebv", returnFunc = identity, verbose = TRUE, skip = NULL){
	trials <- paste0("trial", 1:nTrial)
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

	for(i in 1:(nYr + length(trials))){
		if(verbose) cat("Year:", i, "\n")
		# i = 1
		# predict latest RGSC with updated GS model 
		if(i > 1) {
			RGSC[[gen(GScylce[1]-1)]] <- setEBV(RGSC[[gen(GScylce[1]-1)]], GSmodel[[gen(i-1)]], simParam = simParam)
		 }
		# make selections for DH parents
		selGStoP <- selectInd(RGSC[[length(RGSC)]], nInd = nDHfam, trait = 1, use = "ebv") # does this select from specific families? Almost certainly.
		# selGStoP@id
		
		# get generation indicies
		genBack <- if(i > 1) tail(5:1, min(5, i-1)) else 0
		genI <- if(i > 1) tail(1:(i-1), min(5, i-1)) else 0
		index <- if(i > 1) 1:length(genI) else 1

		# make DH families
		VDP[[trials[1]]][[gen(i)]] <- makeDH(selGStoP, nDH = DHfamSize)
		
		# phenotype DH in first trial if not skipped
		if(!trials[1] %in% skip) VDP[[trials[1]]][[gen(i)]] <- setPheno(VDP[[trials[1]]][[gen(i)]], varE = h2toVe(h2[[1]], Vg), reps = 1)

		# print mean genotypic value of DH 
		if(verbose) print(sapply(VDP[[trials[1]]], function(x) mean(gv(x))))

		# 
		if(select == "ebv" | !is.null(skip)) if(verbose) cat("updating ebvs for generations:", genI, "\n")
		
		if(i > 1) {
			ii <- i
			while(ii >= max(1, i - nTrial)){
				ii <- ii - 1
				# count <- count + 1
				# set EBVs using the GS model from previous year
				for(g in index) {
					gi <- genI[g]
					ti <- trials[genBack[g]]
					tip1 <- trials[genBack[g] + 1]
					h2i <- h2[genBack[g]]
					li <- trialReps[genBack[g]]
					ri <- trialLocs[genBack[g]]
					# if(select == "ebv" | !is.null(skip)) VDP[[ti]][[gen(gi)]] <- setEBV(VDP[[ti]][[gen(gi)]], GSmodel[[gen(ii)]], simParam = simParam)
					if(select == "ebv" | !is.null(skip)) VDP[[ti]][[gen(gi)]] <- setEBV(VDP[[ti]][[gen(gi)]], GSmodel[[gen(ii)]], simParam = simParam)
					sel <- if(ti %in% skip) "ebv" else  select
					VDP[[ti]][[gen(ii)]] <- selectInd(VDP[[ti]][[gen(ii)]], nInd = selOutOfHR, trait = 1, use = sel, returnPop = TRUE)
					if(!tip1 %in% skip) VDP[[tip1]][[gen(ii)]] <- setPheno(VDP[[tip1]][[gen(ii)]], varE = h2toVe(h2i, Vg), reps = li * ri) 	
				}
			}
		} 

		# # count <- 0
		# if(i == 1) {

		# 	if(select == "ebv" | !is.null(skip)) VDP[[ti]][[gen(gi)]] <- setEBV(VDP[[ti]][[gen(gi)]], GSmodel[[gen(ii)]], simParam = simParam)
		# 	sel <- if(ti %in% skip) "ebv" else  select
		# 	VDP[[ti]][[gen(ii)]] <- selectInd(VDP[[ti]][[gen(ii)]], nInd = selOutOfHR, trait = 1, use = sel, returnPop = TRUE)
		# 	if(!tip1 %in% skip) VDP[[tip1]][[gen(ii)]] <- setPheno(VDP[[tip1]][[gen(ii)]], varE = h2toVe(h2i, Vg), reps = li * ri) 	

		# } else {
		# 	ii <- i
		# 	while(ii >= max(1, i - nTrial)){
		# 		ii <- ii - 1
		# 		# count <- count + 1
		# 		# set EBVs using the GS model from previous year
		# 		for(g in index) {
		# 			gi <- genI[g]
		# 			ti <- trials[genBack[g]]
		# 			tip1 <- trials[genBack[g] + 1]
		# 			h2i <- h2[genBack[g]]
		# 			li <- trialReps[genBack[g]]
		# 			ri <- trialLocs[genBack[g]]
		# 			if(select == "ebv" | !is.null(skip)) VDP[[ti]][[gen(gi)]] <- setEBV(VDP[[ti]][[gen(gi)]], GSmodel[[gen(ii)]], simParam = simParam)
		# 			sel <- if(ti %in% skip) "ebv" else  select
		# 			VDP[[ti]][[gen(ii)]] <- selectInd(VDP[[ti]][[gen(ii)]], nInd = selOutOfHR, trait = 1, use = sel, returnPop = TRUE)
		# 			if(!tip1 %in% skip) VDP[[tip1]][[gen(ii)]] <- setPheno(VDP[[tip1]][[gen(ii)]], varE = h2toVe(h2i, Vg), reps = li * ri) 	
		# 		}
		# 	}
		# } 

		# # make DH families
		# VDP[["headrow"]][[gen(i)]] <- makeDH(selGStoP, nDH = DHfamSize)

		# # phenotype population for 1 year
		# # initial trials (headrows), traits are assumed to have half te heritability of normal plots. 
		# if(!"headrow" %in% skip) VDP[["headrow"]][[gen(i)]] <- setPheno(VDP[["headrow"]][[gen(i)]], varE = h2toVe(h2hr, Vg), reps = 1)
		# if(verbose) print(sapply(VDP[["headrow"]], function(x) mean(gv(x))))

		# # So far this assumes that we only consider phenotypes from the latest trial... 
		# # I dont know how to keep multiple generations without having a bunch of correlated traits and using a (evenly weighted) selection index.
		# # The problem is then there will be missing values, and this program cant handle that (at least I dont think so)



		if(i > 1) {
			genBack <- tail(5:1, min(5, i-1))
			genI <- tail(1:(i-1), min(5, i-1))
			index <- 1:length(genI)
			# predict performance if seleciton on ebv or skip stages 
			if(select == "ebv" | !is.null(skip)){
				if(verbose) cat("updating ebvs for generations:", genI, "\n")
				for(g in index) {
					gi <- genI[g]
					ti <- trials[genBack[g]]
					VDP[[ti]][[gen(gi)]] <- setEBV(VDP[[ti]][[gen(gi)]], GSmodel[[gen(i-1)]], simParam = simParam)
				}
			}
			sel <- if("headrow" %in% skip) "ebv" else  select
			VDP[["prelim"]][[gen(i-1)]] <- selectInd(VDP[["headrow"]][[gen(i-1)]], nInd = selOutOfHR, trait = 1, use = sel, returnPop = TRUE)
			if(!"prelim" %in% skip) VDP[["prelim"]][[gen(i-1)]] <- setPheno(VDP[["prelim"]][[gen(i-1)]], varE = h2toVe(h2, Vg), reps = nLocPrelim * nRepsPerLocPrelim) # need to adjust replicaitons!
		}
		if(i > 2) {
			sel <- if("prelim" %in% skip) "ebv" else  select
			VDP[["advance"]][[gen(i-2)]] <- selectInd(VDP[["prelim"]][[gen(i-2)]], nInd = selOutOf1plot, trait = 1, use = sel, returnPop = TRUE)
			if(!"advance" %in% skip) VDP[["advance"]][[gen(i-2)]] <- setPheno(VDP[["advance"]][[gen(i-2)]], varE = h2toVe(h2, Vg), reps = nLocAdv * nRepsPerLocAdv)
		}
		if(i > 3) {
			sel <- if("advance" %in% skip) "ebv" else  select
			VDP[["elite1"]][[gen(i-3)]] <- selectInd(VDP[["advance"]][[gen(i-3)]], nInd = selOutOf3plot, trait = 1, use = sel, returnPop = TRUE)
			if(!"elite1" %in% skip) VDP[["elite1"]][[gen(i-3)]] <- setPheno(VDP[["elite1"]][[gen(i-3)]], varE = h2toVe(h2, Vg), reps = nLocMET * nRepsPerLocMET)
		}
		if(i > 4) {
			sel <- if("elite1" %in% skip) "ebv" else  select
			VDP[["elite2"]][[gen(i-4)]] <- VDP[["elite1"]][[gen(i-3)]]
			if(!"elite2" %in% skip) VDP[["elite2"]][[gen(i-4)]] <- setPheno(VDP[["elite1"]][[gen(i-4)]], varE = h2toVe(h2, Vg), reps = nLocMET * nRepsPerLocMET)
		}
		if(i > 5) {
			sel <- if("elite1" %in% skip | "elite2" %in% skip) "ebv" else  select
			VDPelite <- c(VDP[["elite1"]][[gen(i-4)]], VDP[["elite2"]][[gen(i-4)]])
			VDP[["variety"]][[gen(i-5)]] <- selectInd(VDPelite, nInd = selVariety, trait = 1, use = sel, returnPop = TRUE) 
		}

		if(i <= nYr){
			# run GS model to cycle through RGSC for year i
			for(j in GScylce){
				# if(j > 1) RGSC[[gen(j - 1)]] <- setEBV(RGSC[[gen(j-1)]], GSmodel[[gen(i - 1)]], simParam = simParam)
				if(j != GScylce[1]) RGSC[[gen(j - 1)]] <- setEBV(RGSC[[gen(j-1)]], GSmodel[[gen(i - 1)]], simParam = simParam)
				RGSC[[gen(j)]] <- selectCross(pop = RGSC[[gen(j-1)]], nInd = RGSC[[gen(j-1)]]@nInd * RGSCintensity, 
											   use = RGSCuse,  trait = 1, simParam = simParam, nCrosses = nNuclear, nProgeny = 1) 
			}
			GScylce <- GScylce + GScylcePerYr

			# return lines from VDP into the RGSC 
			if(!is.null(returnVDPtoRGSC) & i > 1){
				returnToRGSC <- trials[genBack] %in% returnVDPtoRGSC
				addToRGSC <- list()
				for(g in index[returnToRGSC]) {
					gi <- genI[g]
					ti <- trials[genBack[g] + 1]
					addToRGSC[[g]] <- VDP[[ti]][[gen(gi)]]
				}
				addToRGSC <- Reduce(c, addToRGSC)
				RGSC[[gen(GScylce[1] - 1)]] <- c(RGSC[[gen(GScylce[1] - 1)]], addToRGSC)
			}

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
