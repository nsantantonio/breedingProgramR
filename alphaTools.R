# getArgs <- function(defaultArgs) {
# 	args <- commandArgs(TRUE)
# 	isAssn <- grepl("=", args)
# 	userArgs <- args[isAssn]
# 	argSplit <- strsplit(userArgs, "=")
# 	argList <- lapply(argSplit, "[[", 2)
# 	names(argList) <- lapply(argSplit, "[[", 1)
# 	print(argList)
# 	defaultArgs[names(argList)] <- argList
# 	defaultArgs
# }


h2toVe <- function(h2, Vg = 1) Vg * (1-h2) / h2
gen <- function(i) paste0("gen", i)
logSelInd <- function(pop, sel) pop@id %in% sel  
rSel <- function(sel) Reduce("&", sel) 
maxBv <- function(simParam, traits = 1) sapply(simParam$traits[traits], function(x) sum(abs(x@addEff)))
sdUnCor <- function(x) sqrt(mean(x^2) - mean(x)^2)
getRRh2 <- function(rrFit) solve(pop0pred@Vu + pop0pred@Ve) %*% pop0pred@Vu
getAcc <- function(pop) cor(gv(pop), ebv(pop))

# get locus index
pullLoci <- function(simParam, snpChip = 1, asList = FALSE) 	{
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
	if(returnMatrix) geno <- do.call(cbind, geno)
	geno
}


rlapply <- function(l, f = identity, level = 1, combine = list, counter = 1, ...){
	args <- list(...)
	if(counter < level){
		do.call(lapply, c(list(X = l, FUN = rlapply, f = f, level = level, combine = combine, counter = counter + 1), args))
	} else {
		result <- do.call(lapply, c(list(X = l, FUN = f), args))
		if(identical(combine, list)) return(result) else return(do.call(combine, result))
	}
}



getSel <- function(selCrit, n) names(sort(selCrit, decreasing = TRUE))[1:n]
dummyFunc <- function(x, retrn) retrn

# function to return expected quantiles from sampling DH individuals
simDHdist <- function(pop, GSfit, DHsampleSize = 200, returnQuantile = 0.9) {
	qDH <- function(i){
		DH <- makeDH(pop[i], nDH = DHsampleSize)
		DH <- setEBV(DH, GSfit, simParam = simParam)
		DHebv <- ebv(DH)
		quantile(DHebv, returnQuantile)[[1]]
	}
	sapply(pop@id, qDH)
}


getPopMeanVar <- function(parVal, parCov, Vg){
	if(ncol(parVal) > 1) stop("cannot use more than 1 trait...") 
	pbar <- combn(parVal, 2, mean)
	pCovar <- parCov[lower.tri(parCov)] * Vg
	pVarSum <- combn(diag(parCov) * Vg, 2, sum) 
	crossvar <- pVarSum - 2 * pCovar
	list(pbar = pbar, crossvar = crossvar, pcovar = pCovar)
}

vanRaden1 <- function(M){
	Z <- scale(M, scale = FALSE)
	p <- attributes(Z)[["scaled:center"]] / 2
	ZZt <- tcrossprod(Z)
	ZZt / (2 * crossprod(p, 1-p)[[1]])
}

getExpDist <- function(pop, GSfit, eSelInt, pullGeno = pullSnpGeno, Gvar = varA, DHsampleSize = 200, returnQuantile = 0.9) {
	parVal <- ebv(pop)
	rownames(parVal) <- pop@id
	K <- vanRaden1(pullGeno(pop))
	Vg <- Gvar(pop)
	if(prod(dim(Vg)) > 1) stop("can only handle a single trait!") else Vg <- Vg[[1]]
	pE <- getPopMeanVar(parVal, K, Vg)
	do.call(rbind, combn(pop@id, 2, simplify = FALSE))
	pE$pbar + qnorm(1 - eSelInt, sd = sqrt(pE$crossvar)) # does this make sense?
	# cor(pE$pbar, pE$pbar + qnorm(1 - eSelInt, sd = sqrt(pE$crossvar)))
}


crossPlanFunc <- function(pop, nFam, famSize = 1){ # note this is just a random sampler, to illustrate how one might build a function to pick pairs. 
	allCrosses <- combn(pop@id, 2)
	resample <- if(nFam > ncol(allCrosses)) TRUE else FALSE 
	crosses <- allCrosses[, sample(1:ncol(allCrosses), nFam, replace = resample)]
	if(famSize > 1) crosses[, rep(1:nFam, each = famSize)]
	t(crosses)
}