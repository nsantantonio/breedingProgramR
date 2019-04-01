getArgs <- function(defaultArgs) {
	args <- commandArgs(TRUE)
	isAssn <- grepl("=", args)
	userArgs <- args[isAssn]
	argSplit <- strsplit(userArgs, "=")
	argList <- lapply(argSplit, "[[", 2)
	names(argList) <- lapply(argSplit, "[[", 1)
	print(argList)
	defaultArgs[names(argList)] <- argList
	defaultArgs
}

# get parent directory
parDir <- getwd()

if(system("hostname", intern = TRUE) == "cbsurobbins.biohpc.cornell.edu") system("export OMP_NUM_THREADS=1")

# set default arguments
defArgs <- list(paramFile = "testparam.R", founderRData = "founderPop/testAlphaSimR1000SegSite.RData", 
				simFunc = "fitFuncSingleVariate.R", nThreads = 10, simName = "simNoTitle")
defArgs <- getArgs(defArgs)
attach(defArgs)

# load libraries
library(AlphaSimR)
library(doMC)
library(txtplot) 
source("alphaTools.R")
source(paramFile)
source(simFunc)


# simName <- "RGSCnoselF2redo"
# RGSCselect <- FALSE
# selF2 <- TRUE
# nF2 <- 100
# selQuantile = TRUE
# ssd = FALSE



# make directory to store all results
simDir <- paste0(parDir, "/", simName) 
system(paste0("simdir=", simDir, "
if [ ! -d $simdir ]; then
	mkdir $simdir
else 
	echo 'Directory already exists! Any results therein will be overwritten without further notice...'
fi"))


# set seed
if(!is.null(seed)) set.seed(seed)


# create simple founder pop, either simple or from loaded MaCS sim
if(simpleFounder) {
	founderPop <- quickHaplo(nFounder, nChrom, nLoci, genLen = 1, ploidy = 2L, inbred = FALSE)
} else {
	load(founderRData)
	founderPop <- founderPop[1:nFounder]
}

# Setting Simulation Parameters
SP = SimParam$new(founderPop)
# SP$nThreads <- 40

# add trait 
SP$addTraitA(nQtlPerChr=nQTL, mean = 0, var = Vg)
# SP$addTraitAG(nQtlPerChr=nQTL, mean = 0, var = Vg, varGxE = Vgxe, varE = h2toVe(h2), corGxE = )

# add error variance as function of heritability
SP$setVarE(h2 = founderh2)

# add genotype platform
SP$addSnpChip(nM)

loci <- pullLoci(SP)
# Reduce(intersect, loci)
# simParam <- SP

# nYr <- 6
# run1 <- sim(founderPop, simParam = SP, select = "ebv")



registerDoMC(nThreads)

simrun <- foreach(r = 1:reps) %dopar% sim(founderPop, simParam = SP, select = "pheno")
save(simrun, SP, file = paste0(simName, ".RData"))

# simParam <- SP; select = "ebv"; returnFunc = identity; verbose = TRUE; skip = NULL


# h2toVe <- function(h2, Vg = 1) Vg * (1-h2) / h2
# gen <- function(i) paste0("gen", i)
# logSelInd <- function(pop, sel) pop@id %in% sel  
# rSel <- function(sel) Reduce("&", sel) 
# maxBv <- function(simParam, traits = 1) sapply(simParam$traits[traits], function(x) sum(abs(x@addEff)))
# maxBv(SP)

# fitBP <- function(founderPop, simParam, select = "ebv", returnFunc = identity, verbose = TRUE, skip = NULL){
# 	if(!is.null(returnVDPtoRGSC)) if(!(returnVDPtoRGSC %in% trials[-length(trials)])) stop("'returnVDPtoRGSC' argument can only take values of 'headrow', 'prelim', 'advance', 'elite1', or 'elite2', to return selected lines out of those trials into the RGSC, e.g. to return identified varieties, use 'elite2'") 
# 	pop0 <- newPop(founderPop)
# # i commented this out, not sure it needs to be set?
# 	# pheno(pop0)
# 	GScylce <- 1:GScylcePerYr

# 	RGSC <- list()
# 	GSmodel <- list()
# 	trials <- c("headrow", "prelim", "advance", "elite1", "elite2", "variety")
# 	names(trials) <- trials
# 	VDP <- lapply(trials, function(x) list())

# 	# initialize nuclear population
# 	RGSC[[gen(0)]] <- pop0
# 	GSmodel[[gen(0)]] <- RRBLUP(RGSC[[gen(0)]], traits = 1, use = "pheno", snpChip = 1, simParam = simParam)
# 	RGSC[[gen(0)]] <- setEBV(RGSC[[gen(0)]], GSmodel[[gen(0)]], simParam = simParam)
# 	# getAcc(RGSC[[gen(0)]])

# 	for(i in 1:(nYr + 5)){
# 	# for(i in 1:(nYr)){
# 		if(verbose) cat("Year:", i, "\n")
# 		# i = 1
# 		# predict latest RGSC with updated GS model 
# 		if(i > 1) {
# 			RGSC[[gen(GScylce[1]-1)]] <- setEBV(RGSC[[gen(GScylce[1]-1)]], GSmodel[[gen(i-1)]], simParam = simParam)
# 		 }
# 		# make selections for DH parents
# 		selGStoP <- selectInd(RGSC[[length(RGSC)]], nInd = nDHfam, trait = 1, use = "ebv") # does this select from specific families? Almost certainly.
# 		# selGStoP@id
		
# 		# make DH families
# 		VDP[["headrow"]][[gen(i)]] <- makeDH(selGStoP, nDH = DHfamSize)


# 		# phenotype population for 1 year
# 		# initial trials (headrows), traits are assumed to have half te heritability of normal plots. 
# 		if(!"headrow" %in% skip) VDP[["headrow"]][[gen(i)]] <- setPheno(VDP[["headrow"]][[gen(i)]], varE = h2toVe(h2hr, Vg), reps = 1)
# 		if(verbose) print(sapply(VDP[["headrow"]], function(x) mean(gv(x))))

# 		# So far this assumes that we only consider phenotypes from the latest trial... 
# 		# I dont know how to keep multiple generations without having a bunch of correlated traits and using a (evenly weighted) selection index.
# 		# The problem is then there will be missing values, and this program cant handle that (at least I dont think so)

# 		if(i > 1) {
# 			genBack <- tail(5:1, min(5, i-1))
# 			genI <- tail(1:(i-1), min(5, i-1))
# 			index <- 1:length(genI)
# 			# predict performance if seleciton on ebv or skip stages 
# 			if(select == "ebv" | !is.null(skip)){
# 				if(verbose) cat("updating ebvs for generations:", genI, "\n")
# 				for(g in index) {
# 					gi <- genI[g]
# 					ti <- trials[genBack[g]]
# 					VDP[[ti]][[gen(gi)]] <- setEBV(VDP[[ti]][[gen(gi)]], GSmodel[[gen(i-1)]], simParam = simParam)
# 				}
# 			}
# 			sel <- if("headrow" %in% skip) "ebv" else  select
# 			VDP[["prelim"]][[gen(i-1)]] <- selectInd(VDP[["headrow"]][[gen(i-1)]], nInd = selOutOfHR, trait = 1, use = sel, returnPop = TRUE)
# 			if(!"prelim" %in% skip) VDP[["prelim"]][[gen(i-1)]] <- setPheno(VDP[["prelim"]][[gen(i-1)]], varE = h2toVe(h2, Vg), reps = nLocPrelim * nRepsPerLocPrelim) # need to adjust replicaitons!
# 		}
# 		if(i > 2) {
# 			sel <- if("prelim" %in% skip) "ebv" else  select
# 			VDP[["advance"]][[gen(i-2)]] <- selectInd(VDP[["prelim"]][[gen(i-2)]], nInd = selOutOf1plot, trait = 1, use = sel, returnPop = TRUE)
# 			if(!"advance" %in% skip) VDP[["advance"]][[gen(i-2)]] <- setPheno(VDP[["advance"]][[gen(i-2)]], varE = h2toVe(h2, Vg), reps = nLocAdv * nRepsPerLocAdv)
# 		}
# 		if(i > 3) {
# 			sel <- if("advance" %in% skip) "ebv" else  select
# 			VDP[["elite1"]][[gen(i-3)]] <- selectInd(VDP[["advance"]][[gen(i-3)]], nInd = selOutOf3plot, trait = 1, use = sel, returnPop = TRUE)
# 			if(!"elite1" %in% skip) VDP[["elite1"]][[gen(i-3)]] <- setPheno(VDP[["elite1"]][[gen(i-3)]], varE = h2toVe(h2, Vg), reps = nLocMET * nRepsPerLocMET)
# 		}
# 		if(i > 4) {
# 			sel <- if("elite1" %in% skip) "ebv" else  select
# 			VDP[["elite2"]][[gen(i-4)]] <- VDP[["elite1"]][[gen(i-3)]]
# 			if(!"elite2" %in% skip) VDP[["elite2"]][[gen(i-4)]] <- setPheno(VDP[["elite1"]][[gen(i-4)]], varE = h2toVe(h2, Vg), reps = nLocMET * nRepsPerLocMET)
# 		}
# 		if(i > 5) {
# 			sel <- if("elite1" %in% skip | "elite2" %in% skip) "ebv" else  select
# 			VDPelite <- c(VDP[["elite1"]][[gen(i-4)]], VDP[["elite2"]][[gen(i-4)]])
# 			VDP[["variety"]][[gen(i-5)]] <- selectInd(VDPelite, nInd = selVariety, trait = 1, use = sel, returnPop = TRUE) 
# 		}

# 		if(i <= nYr){
# 			# run GS model to cycle through RGSC for year i
# 			for(j in GScylce){
# 				# if(j > 1) RGSC[[gen(j - 1)]] <- setEBV(RGSC[[gen(j-1)]], GSmodel[[gen(i - 1)]], simParam = simParam)
# 				if(j != GScylce[1]) RGSC[[gen(j - 1)]] <- setEBV(RGSC[[gen(j-1)]], GSmodel[[gen(i - 1)]], simParam = simParam)
# 				RGSC[[gen(j)]] <- selectCross(pop = RGSC[[gen(j-1)]], nInd = RGSC[[gen(j-1)]]@nInd * RGSCintensity, 
# 											   use = RGSCuse,  trait = 1, simParam = simParam, nCrosses = nNuclear, nProgeny = 1) 
# 			}
# 			GScylce <- GScylce + GScylcePerYr

# 			# return lines from VDP into the RGSC 
# 			if(!is.null(returnVDPtoRGSC) & i > 1){
# 				returnToRGSC <- trials[genBack] %in% returnVDPtoRGSC
# 				addToRGSC <- list()
# 				for(g in index[returnToRGSC]) {
# 					gi <- genI[g]
# 					ti <- trials[genBack[g] + 1]
# 					addToRGSC[[g]] <- VDP[[ti]][[gen(gi)]]
# 				}
# 				addToRGSC <- Reduce(c, addToRGSC)
# 				RGSC[[gen(GScylce[1] - 1)]] <- c(RGSC[[gen(GScylce[1] - 1)]], addToRGSC)
# 			}

# 			trnSet <- lapply(VDP[trials[!grepl("variety", trials)]], function(x) x[names(x) %in% gen((i-max(1, lgen)):i)])
# 	 		hasPop <- sapply(trnSet, length) > 0
# 			train <- Reduce(c, lapply(trnSet[hasPop], function(x) Reduce(c, x)))
# 			cat("training set has ", train@nInd, "individuals...\n")	
# 			GSmodel[[gen(i)]] <- RRBLUP(train, traits = 1, use = "pheno", snpChip = 1, simParam=simParam)
# 		} else {
# 			cat("final year reached, selecting on phenotypes / ebv trainde with last year training set ...\n")
# 			GSmodel[[gen(i)]] <- GSmodel[[gen(i-1)]]
# 		}
# 	}
# 	rL <- returnFunc(list(RGSC = RGSC, VDP = VDP, GSmodel = GSmodel))
# 	return(rL)
# }















# setMKLthreads(1)      
# simrun <- foreach(r = 1:reps) %dopar% fitBP(founderPop, simParam, select = "ebv")


# rlapply()

# RGSCgv <- sapply(RGSC, function(x) mean(gv(x)))
# HRgv <- sapply(VDP[["headrow"]], function(x) mean(gv(x)))
# VARgv <- sapply(VDP[["variety"]], function(x) mean(gv(x)))

# library(doMC)
# registerDoMC(2)
# foreach(i=1:10) %dopar% crossprod(matrix(rnorm(10000 * 10000), 10000, 10000))[[1]]


# foreach()

# function(select)


# 		if(i > 1) {
# 			# gen1 <- VDP[[gen(i-1)]]
# 			# initPheno <- pheno(gen1)
# 			# gen1 <- setPheno(gen1, varE = h2toVe(h2, Vg), reps = nRepsPerLoc)
# 			# secPheno <- pheno(gen1)
# 			# cor(secPheno, initPheno)

# 			VDP[[gen(i-1)]] <- setPheno(VDP[[gen(i-1)]], varE = h2toVe(h2, Vg), reps = nRepsPerLoc)

# 			# still need to work on this part. 
# 			GSmodel[[gen(0)]] <- RRBLUP(RGSC[[gen(0)]], traits = 1, use = "pheno", snpChip = 1, simParam=SP)

# 			VDP[[gen(j-1)]] <- setEBV(RGSC[[gen(j-1)]], GSmodel[[gen(i - 1)]], simParam = SP)

# 			selectInd(VDP[[gen(j-1)]]) 

# 		}


# 		train <- do.call(c, VDP)
# 		if(class(train) == "list") train <- train[[1]]
# 		GSmodel[[gen(i)]] = RRBLUP(train, traits = 1, use = "pheno", snpChip = 1, simParam=SP)



# 		# FOR SOME REASON h2 IS HIGHER THAN EXPECTED?
# 		# r <- NULL
# 		# while(length(r) < 100){
# 		# 	# VDP[[gen(i)]] <- setPheno(VDP[[gen(i)]], varE = h2toVe(h2 / 2, Vg), reps = 1)
# 		# 	VDP[[gen(i)]] <- setPheno(VDP[[gen(i)]], varE = h2toVe(h2, Vg), reps = 1)
# 		# 	r <- c(r, cor(VDP[[gen(i)]]@pheno, VDP[[gen(i)]]@gv)^2)
# 		# }
# 		# mean(r)

# 		# VDP[[gen(i)]]@pheno

# 		# 1 location trials


# 		# 3 location trials
# 		if(i > 2){
# 			VDP[[gen(i-1)]] <- setPheno(VDP[[gen(i-1)]], h2 = h2, reps = nRepsPerLoc)

# 			setPhen0(x)
# 		}
			
# 		# MET1
# 		if(i > 3)

# 		# MET2
# 		if(i > 4)


# 		# 1 location trials
# 		if(i > 1){
# 			train <- do.call(c, VDP)

# 			VDP[[gen(i-1)]] <- setPheno(VDP[[gen(i-1)]], h2 = h2, reps = nRepsPerLoc)
# 			# still need to work on this part. 
# 			GSmodel[[gen(0)]] <- RRBLUP(RGSC[[gen(0)]], traits = 1, use = "pheno", snpChip = 1, simParam=SP)

# 			VDP[[gen(j-1)]] <- setEBV(RGSC[[gen(j-1)]], GSmodel[[gen(i - 1)]], simParam = SP)

# 			selectInd(VDP[[gen(j-1)]]) 

# 		}

# 		# 3 location trials
# 		if(i > 2){
# 			VDP[[gen(i-1)]] <- setPheno(VDP[[gen(i-1)]], h2 = h2, reps = nRepsPerLoc)

# 			setPhen0(x)
# 		}
			
# 		# MET1
# 		if(i > 3)

# 		# MET2
# 		if(i > 4)

# 		# update GS model
# 		# should include original phenotypes from founder pop, or only use for first round?
# 		train <- do.call(c, VDP)
# 		if(class(train) == "list") train <- train[[1]]
# 		train@nInd
# 		GSmodel[[gen(i)]] = RRBLUP(train, traits = 1, use = "pheno", snpChip = 1, simParam=SP)

# 		setEBV RGSC[[length(RGSC)]]


# 		bvsIndex[[r]][[i]] <- colMeans(gv(poplIndex[[i]]))
# 		bvs1trait[[r]][[i]] <- colMeans(gv(popl1trait[[i]]))

# 		allPheno <- do.call(c, VDP)

# 		if(i %% 3 == 0){
# 			pop1 = selectCross(pop=pop0, nInd=pop0@nInd * intensity, use = "ebv",  trait=1, simParam=SP, nCrosses=10, nProgeny = 1) 

# 			if(skipStage %in% 1){
# 					RRBLUP
# 				} else {
# 					if(sparseStage %in% 1)
# 					setPheno
# 				}


# 			if(skipStage %in% 2){
# 					RRBLUP
# 				} else {
# 					if(sparseStage %in% 2)
# 					setPheno
# 				}


# 			if(skipStage %in% 3){
# 					RRBLUP
# 				} else {
# 					if(sparseStage %in% 3)
# 					setPheno
# 				}

# 			# update training pop	!
# 		}


# 	}
# 	bvs1trait[[r]] <- data.frame(rep = r, gen = generation, do.call(rbind, bvs1trait[[r]]))
# 	bvsIndex[[r]] <- data.frame(rep = r, gen = generation, do.call(rbind, bvsIndex[[r]]))
# }
# bvs1traitDf <- do.call(rbind, bvs1trait)
# bvsIndexDf <-  do.call(rbind, bvsIndex)

# meanGV1trait <- data.frame(trait1 = tapply(bvs1traitDf$X1, bvs1traitDf$gen, mean),
# 			   trait1sd = tapply(bvs1traitDf$X1, bvs1traitDf$gen, sd),
# 			   trait2 = tapply(bvs1traitDf$X2, bvs1traitDf$gen, mean),
# 			   trait2sd = tapply(bvs1traitDf$X2, bvs1traitDf$gen, sd),
# 			   trait3 = tapply(bvs1traitDf$X3, bvs1traitDf$gen, mean),
# 			   trait3sd = tapply(bvs1traitDf$X3, bvs1traitDf$gen, sd))
# meanGVindex <- data.frame(trait1 = tapply(bvsIndexDf$X1, bvsIndexDf$gen, mean),
# 			   trait1sd = tapply(bvsIndexDf$X1, bvsIndexDf$gen, sd),
# 			   trait2 = tapply(bvsIndexDf$X2, bvsIndexDf$gen, mean),
# 			   trait2sd = tapply(bvsIndexDf$X2, bvsIndexDf$gen, sd),
# 			   trait3 = tapply(bvsIndexDf$X3, bvsIndexDf$gen, mean),
# 			   trait3sd = tapply(bvsIndexDf$X3, bvsIndexDf$gen, sd))




# ylims <- range(cbind(meanGV1trait[paste0("trait", 1:3)],meanGVindex[paste0("trait", 1:3)] )) * 1.1

# pdf("selIndexVs1Trait.pdf")
# plot(NA, xlim = c(0, 20), ylim = ylims, ylab = "Mean Genetic Value", xlab = "Generation")
# for(i in 1:nTrait) {
# 	lines(generation, meanGV1trait[[paste0("trait", i)]], col = cols[i], lty = 1, lwd = 1.75)
# 	lines(generation, meanGVindex[[paste0("trait", i)]], col = cols[i], lty = 2, lwd = 1.75)
# }
# legend("topleft", legend = c("trait 1", "trait 2", "trait 3", "index") , col = cols[c(1:3, 1)], lty = c(1, 1, 1, 2))
# dev.off()





# pop1 = selectCross(pop=pop0, nInd=pop0@nInd * intensity, use = "ebv",  trait=selIndex, simParam=SP, b=b, nCrosses=1000, nProgeny = 1) # number of families and projeny also needs to be considered here. Especially within some reasonable amount of crosses to be made
# # set EBV based on pop0 
# pop1 = setEBV(pop1, pop0pred, simParam=SP)
# ebv(pop1)


# # pop = selectInd(pop, nInd=100, trait=selIndex, simParam=SP, b=b)
# pop = selectInd(pop, nInd=pop0@nInd * intensity, trait=selIndex, simParam=SP, b=b)

# pop1 




# M <- pullSnpGeno(pop0, simParam=SP)
# G <- pullQtlGeno(pop0, simParam=SP)
# dim(M)
# dim(G)


# founderPop


# G <- pullIbdHaplo(pop1, simParam=SP)


# # nt <- 2

# # if(nt == 2){
# # 	Rho <- diag(2)
# # 	Rho[upper.tri(Rho)] <- c(-0.3)
# # 	Rho[lower.tri(Rho)] <- t(Rho)[lower.tri(Rho)]
# # 	econWt = c(0.5, 0.5)
# # 	SP$addTraitA(nQtlPerChr=c(2*(9:1), 1), mean = rep(0, 2), var = rep(1, 2), corA = Rho)
# # 	SP$setVarE(h2=c(0.3, 0.5))
# # } else {
# # }



# # Modeling the Breeding Program
# pop <- newPop(founderPop)

# pop@genMap


# genMean <- meanG(pop)
# generation <- 1:20
# for(i in generation){
# 	SP$addSnpChip(10)


# 	setPheno()

# 	b = smithHazel(econWt, varG(pop), varP(pop))

# 	# pop = selectInd(pop, nInd=100, trait=selIndex, simParam=SP, b=b)
# 	pop = selectInd(pop, nInd=100, trait=selIndex, simParam=SP, b=b)
#   # pop = selectCross(pop=pop, nFemale=500, nMale=25, use="gv", nCrosses=1000)
#   pop = selectCross(pop=pop, nInd=100, trait=selIndex, simParam=SP, b=b, nCrosses=1000)
#   genMean = rbind(genMean, meanG(pop))
# }

# ylims <- range(genMean) * 1.1
# pdf("gainTestAlphaSR.pdf")
# plot(NA, ylim = ylims, xlim = c(0, nGen), xlab = "Generation", ylab = "Value", main = "Genetic Gain")

# nGen <- nrow(genMean)
# nTrait <- ncol(genMean)
# legend("topleft", legend = c("Quant", "Min", "Int", "Bin"), lwd = 2, col = cols[1:nTrait], bty = "n")
# for(i in 1:nTrait){
# 	lines(1:nGen - 1, genMean[, i], col = cols[[i]], lwd = 2)
# }
# dev.off()

# # Examining the Results
# txtplot(0:20, genMean)#, xlab="Generation", ylab="Mean Genetic Value")












# # SP$setGender("yes_sys") # not way to make both?



# # Modeling the Breeding Program
# pop <- newPop(founderPop)
# genMean <- meanG(pop)
# generation <- 1:20
# for(i in generation){
# 	b = smithHazel(econWt, varG(pop), varP(pop))

# 	# pop = selectInd(pop, nInd=100, trait=selIndex, simParam=SP, b=b)
# 	pop = selectInd(pop, nInd=100, trait=selIndex, simParam=SP, b=b)
#   # pop = selectCross(pop=pop, nFemale=500, nMale=25, use="gv", nCrosses=1000)
#   pop = selectCross(pop=pop, nInd=100, trait=selIndex, simParam=SP, b=b, nCrosses=1000)
#   genMean = rbind(genMean, meanG(pop))
# }

# ylims <- range(genMean) * 1.1
# pdf("gainTestAlphaSR.pdf")
# plot(NA, ylim = ylims, xlim = c(0, nGen), xlab = "Generation", ylab = "Value", main = "Genetic Gain")

# nGen <- nrow(genMean)
# nTrait <- ncol(genMean)
# legend("topleft", legend = c("Quant", "Min", "Int", "Bin"), lwd = 2, col = cols[1:nTrait], bty = "n")
# for(i in 1:nTrait){
# 	lines(1:nGen - 1, genMean[, i], col = cols[[i]], lwd = 2)
# }
# dev.off()

# # Examining the Results
# txtplot(0:20, genMean)#, xlab="Generation", ylab="Mean Genetic Value")



# G = 1.5*diag(2)-0.5 #Genetic correlation matrix
# SP$addTraitA(10, mean=c(0,0), var=c(1,1), corA=G)
# SP$setVarE(h2=c(0.5,0.5))


# for (i in 1:nGen) {
# 	if (i %% 3 == 0) {

# 	}  
# } 