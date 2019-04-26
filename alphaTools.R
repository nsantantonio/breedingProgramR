getComArgs <- function(defaultArgs = NULL) {
  defaults <- !is.null(defaultArgs)
  args <- commandArgs(TRUE)
  isAssn <- grepl("=", args)
  userArgs <- args[isAssn]
  needEval <- grepl("\\(|\\)|\\:|\\'", userArgs) 
  argSplit <- strsplit(userArgs, "=")
  argList <- lapply(argSplit, "[[", 2)
  names(argList) <- lapply(argSplit, "[[", 1)
  argList[needEval] <- lapply(argList[needEval], function(x) eval(parse(text = x)))
  argList[!needEval] <- lapply(argList[!needEval], function(x) strsplit(x, ",")[[1]])
  argList[!needEval] <- type.convert(argList[!needEval], as.is = TRUE)
  print(argList)
  if (defaults){
    defaultArgs[names(argList)] <- argList
    return(defaultArgs)
  } else {
    return(argList)
  }
}


getArgs <- function(f, ...) list(...)[names(list(...)) %in% names(formals(f))]

# popList <- list(1, 2, list(3, 4, list(5, 6, list(7))), list(8, list(9, 10)))

# cR <- function(popList){
#   unlist(lapply(popList, function(x) if (is.list(x)) cR(x) else x), recursive = FALSE)
# }
# cR(popList)

mergePopsRec <- function(popList) { mergePops(lapply(popList, function(x) if (is.list(x)) mergePopsRec(x) else x)) }

getAF <- function(pop, pullGeno = pullSnpGeno) { colMeans(pullGeno(pop)) / 2 }
h2toVe <- function(h2, Vg = 1) { Vg * (1-h2) / h2 }
gen <- function(i) { paste0("gen", i) }
logSelInd <- function(pop, sel) { pop@id %in% sel   }
rSel <- function(sel) { Reduce("&", sel)  }
maxBv <- function(simParam, traits = 1) { sapply(simParam$traits[traits], function(x) sum(abs(x@addEff))) }
sdUnCor <- function(x) { sqrt(mean(x^2) - mean(x)^2) }
getRRh2 <- function(rrFit) { solve(pop0pred@Vu + pop0pred@Ve) %*% pop0pred@Vu }
getAcc <- function(pop) { cor(gv(pop), ebv(pop)) }
dummyFunc <- function(x, retrn) { retrn }

# get locus index
pullLoci <- function(simParam, snpChip = 1, asList = FALSE) {
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
	if (returnMatrix) geno <- do.call(cbind, geno)
	geno
}

rlapply <- function(l, f = identity, level = 1, combine = list, counter = 1, ...){
	args <- list(...)
	if (counter < level){
		do.call(lapply, c(list(X = l, FUN = rlapply, f = f, level = level, combine = combine, counter = counter + 1), args))
	} else {
		result <- do.call(lapply, c(list(X = l, FUN = f), args))
		if (identical(combine, list)) return(result) else return(do.call(combine, result))
	}
}

estIntensity <- function(VDP, i, nT = nTrial, start = "trial1", end = "variety", estFunc = pheno, Gvar = estVg) {
	S <- mean(pheno(VDP[[end]][[gen(i - nT)]])) - mean(pheno(VDP[[start]][[gen(i - nT)]]))
	i <- S / sqrt(Gvar(VDP[[start]][[gen(i - nT)]]))
	if(is.matrix(i) & prod(dim(i)) == 1) i <- i[[1]] else cat("intensity has dimensions:", dim(i), "\n")
	i
}


# getSel <- function(selCrit, n, high = TRUE, variable = "selCrit") {
# 	len <- if(is.data.frame(selCrit)) nrow(selCrit) else length(selCrit)
# 	if(len < n) n <- nrow(selCrit)
# 	if (is.data.frame(selCrit)){
# 		selCrit <- selCrit[order(selCrit[[variable]], decreasing = high), ]
# 		sel <- as.matrix(selCrit[1:n, c("p1", "p2")])
# 	} else {
# 		sel <- names(sort(selCrit, decreasing = high))[1:n]
# 	}
# 	sel
# }

# expand.grid(p1 = c("1", "2", "3", "4"), p2 = c("1", "2", "3", "4"))[sample(1:16), ]
# p <- t(combn(as.character(1:6), 2))
# colnames(p) <- c("p1", "p2")
# dF <- data.frame(p[sample(1:15), ], selCrit = 1:15)

# dF[dFSel(dF, maxP = 2), ]


# tmpDf <- dF
# dF <- tmpDf
dFSel <- function(dF, limit = 1, val = "selCrit", parentCols = c("p1", "p2"), returnPar = TRUE) {
	dFord <- order(dF[[val]])
	if(!(all(dFord == 1:nrow(dF)) | all(dFord == nrow(dF):1))) stop("dF must be sorted in order to select!")
	
	parMat <- as.matrix(dF[, parentCols])
	parents <- sort(unique(c(parMat)))
	parCount <- rep(0, length(parents))
	names(parCount) <- parents
	rows <- NULL

	for (i in 1:nrow(parMat)) {
		parenti <- parMat[i, ]
		parCount[parenti] <- parCount[parenti] + 1
		if(all(parCount[parenti] <= limit)) rows <- c(rows, i)
	}
	if(returnPar) parMat[rows, ] else rows
}

# dFSel(dF, maxP = 2)


getSel <- function(selCrit, n, high = TRUE, variable = "selCrit", parentCols = c("p1", "p2"), maxP = NULL) {
	# add inbreeding limits here!
	len <- if(is.data.frame(selCrit)) nrow(selCrit) else length(selCrit)
	if(len < n) n <- nrow(selCrit)
	if (is.data.frame(selCrit)){
		selCrit <- selCrit[order(selCrit[[variable]], decreasing = high), ]
		if(!is.null(maxP)){
			sel <- dFSel(selCrit, limit = maxP, val = variable, parentCols = parentCols)
			sel <- sel[1:min(nrow(sel), n), ] 
		} else {
			sel <- as.matrix(selCrit[1:n, parentCols])
		}
	} else {
		sel <- names(sort(selCrit, decreasing = high))[1:n]
	}
	sel
}

getTrueQTLeff <- function(simParam, trait = 1) { simParam$traits[[trait]]@addEff }

getTrueBV <- function(pop, simParam, trait = 1) {
	M <- pullQtlGeno(pop)	
	u <- getTrueQTLeff(simParam)
	M %*% u
} 

getSelfVar <- function(M, u, fdiff = NULL) {
	H <- M == 1
	Hu <- H %*% u^2 
	if(!is.null(fdiff)) Hu <- Hu * (1 - 2^(-fdiff))
	# if(ncol(Hu) == 1) Hu <- c(Hu)
	Hu
}

weightedQuantile <- function(mu, sigmasq, quant = 0.9, w = 0.5) {
	cat("            weight parameter 'w' has value:", w, "\n")
	w * mu + (1-w) * qnorm(quant, sd = sqrt(sigmasq))
}

estVg <- function(pop, GSfit = NULL, GSfunc = RRBLUP) {
	if(is.null(GSfit)) GSfit <- GSfunc(pop)
	af <- getAF(pop)
	Vg <- GSfit@Vu
	Vg <- if(is.matrix(Vg) & prod(dim(Vg)) == 1) Vg[[1]]
	Vg * sum(2 * af * (1-af))
}

geth2 <- function(pop, GSfit) {
	Vg <- estVg(pop, GSfit)
	Vg / sum(Vg, GSfit@Ve)
}
# function to return expected quantiles from sampling DH individuals
# select individuals
simDHdist <- function(nSel, pop, GSfit, retQuant = FALSE, quant = 0.9, nDH = 200, w = 0.5, nProgeny = 1, returnPop = TRUE, bigmem = FALSE, ...) {
	if(bigmem) {
		DH <- makeDH(pop, nDH = nDH)
		DH <- setEBV(DH, GSfit)
		DHebv <- ebv(DH)
		simDist <- split(DHebv, rep(1:nInd(pop), each =  nDH))
	} else {
		DHdist <- function(i){
			DH <- makeDH(pop[i], nDH = nDH)
			DH <- setEBV(DH, GSfit)
			ebv(DH)
		}
		simDist <- lapply(pop@id, DHdist)
	}
	expVar <- if(retQuant) sapply(simDist, quantile, probs = quant) else weightedQuantile(sapply(simDist, mean), sapply(simDist, var), quant, w = w)
	names(expVar) <- pop@id
	if(returnPop){	
		selection <- getSel(expVar, n = nSel)
		if(nProgeny > 1) selection <- rep(selection, each = nProgeny)
		return(pop[selection])
	} else {
		return(expVar)
	}
}

# select pairs and cross
simDHdistPairs <- function(nSel, pop, GSfit, nCrosses, retQuant = FALSE, quant = 0.9, nDH = 200, nSimCrosses = 1, nProgeny = 1, verbose = FALSE, ...) {
	nCombos <- choose(nInd(pop), 2) 
	nEx <- if(nCombos < nCrosses) ceiling(nCrosses / nCombos) else 1 
	parents <- do.call(rbind, combn(pop@id, 2, simplify = FALSE))
	colnames(parents) <- c("p1", "p2")
	crosses <- rep(1:nrow(parents), each = nSimCrosses)
	popX <- makeCross(pop, parents[crosses,])
	# popX <- setEBV(popX, GSfit)
	if(verbose) cat("simulating distribution of", nDH, "DH for", nSimCrosses, "crosses for each of", nCombos, "parental pairs\n")
	simQuant <- simDHdist(nSel = nInd(popX), pop = popX, GSfit = GSfit, retQuant = retQuant, quant = quant, nDH = nDH, returnPop = FALSE)
	# if(nSimCrosses > 1) simSd <- tapply(simQuant, crosses, sd)
	if(nSimCrosses > 1) simQuant <- tapply(simQuant, crosses, mean)

	selCrit <- data.frame(parents, selCrit = simQuant)
	selection <- getSel(selCrit, n = nCrosses)
	if(nEx > 1) selection <- selection[rep(1:nrow(selection), times = nEx)[1:nCrosses], ]

	if(nProgeny > 1) selection <- selection[rep(1:nrow(selection), each = nProgeny), ] 
	# makeCross(pop, crossPlan = selection) 
	makeCross(pop, crossPlan = selection) 
}

# select individuals
expDist <- function(nSel, pop, GSfit, use, quant, returnQuant = TRUE, pullGeno = pullSnpGeno, updateEBV = FALSE, w = 0.5, nProgeny = 1, ...){
	expVar <- do.call(getSelfVar, getArgs(getSelfVar, M = pullGeno(pop), u = GSfit@markerEff, ...))
	if(returnQuant) {
		if(updateEBV) pop <- setEBV(pop, GSfit)
		parVal <- ebv(pop)
		expVar <- weightedQuantile(mu = parVal, sigmasq = expVar, quant = quant, w = w)
		# expVar <- w * parVal + (1-w) * qnorm(quant, sd = sqrt(expVar))
	} # need top check why this seems to work so poorly...
	if(ncol(expVar) == 1) expVar <- expVar[, 1]
	selection <- getSel(expVar, n = nSel)
	if(nProgeny > 1) selection <- rep(selection, each = nProgeny)
	pop[selection]
}

# pop = RGSC[[lastRGSCgen]]; GSfit = GSmodel[[lastGSmodel]]; nSel = 20; nCrosses = nNuclear; use = ebv; pullGeno = pullSnpGeno; weightLoci = FALSE; maxCrossPerParent = 1; 
maxVar <- function(pop, GSfit, nSel, nCrosses, use, weightLoci = FALSE, pullGeno = pullSnpGeno, maxCrossPerParent = 1, verbose = FALSE, nProgeny = 1, ...){
	n <- nInd(pop)
	if (n < nSel) nSel <-  n
	nCombos <- choose(nSel, 2)
	nEx <- if(nCombos < nCrosses) ceiling(nCrosses / nCombos) else 1 
	maxP <- if(maxCrossPerParent == 0 | nCombos <  nCrosses) nCrosses else maxCrossPerParent
	
	if(nSel < n) pop <- truncSel(pop, nSel = nSel, use = use)

	M <- pullGeno(pop)
	K <- if(weightLoci) genCov(M, u = c(GSfit@markerEff)) else genCov(M)
	covL <- data.frame(which(lower.tri(K), arr.ind = TRUE), selCrit = K[lower.tri(K)])
	selCrit <- data.frame(covL, p1 = colnames(K)[covL$col], p2 = colnames(K)[covL$row]) 
	
	lenSel <- 0
	while(lenSel < nCrosses / nEx){
		if(lenSel > 0) {
			cat("Not enough possible crosses with maxP =", maxP, "! Increasing maxP to", maxP + 1, "and retrying...\n")
			maxP <- maxP + 1
		}
		selection <- getSel(selCrit, n = nCrosses, high = FALSE, maxP = maxP)
		lenSel <- nrow(selection)
	}
	if(verbose) print(table(selection))
	
	if(nEx > 1) selection <- selection[rep(1:nrow(selection), times = nEx)[1:nCrosses], ]
	if(nProgeny > 1) selection <- selection[rep(1:nrow(selection), each = nProgeny), ] 
	makeCross(pop, crossPlan = selection) 
}


# select pairs
# nSel = selectRGSCi; pop = RGSC[[lastRGSCgen]]; GSfit = GSmodel[[lastGSmodel]]; quant = xInt; nCrosses = nNuclear
# returnQuant = TRUE; weightLoci = FALSE; pullGeno = pullSnpGeno; Gvar = varA; w = 0.5; nProgeny = 1
expDistPairs <- function(pop, GSfit, nSel, quant, nCrosses, use, returnQuant = TRUE, weightLoci = FALSE, pullGeno = pullSnpGeno, maxCrossPerParent = 0, Gvar = estVg, w = 0.5, nProgeny = 1, ...) {
	n <- nInd(pop)
	if (n < nSel) nSel <-  n
	nCombos <- choose(nSel, 2)
	nEx <- if(nCombos < nCrosses) ceiling(nCrosses / nCombos) else 1 
	maxP <- if(maxCrossPerParent == 0 | nCombos <  nCrosses) nCrosses else maxCrossPerParent
	
	if(nSel < n) pop <- truncSel(pop, nSel = nSel, use = use)
	parVal <- ebv(pop)
	rownames(parVal) <- pop@id

	M <- pullGeno(pop)
	K <- if(weightLoci) genCov(M, u = c(GSfit@markerEff)) else genCov(M)
	# txtdensity(K[lower.tri(K)])

	Vg <- Gvar(pop = pop, GSfit = GSfit)
	if (prod(dim(Vg)) > 1) stop("can only handle a single trait!") else Vg <- Vg[[1]]
	parents <- do.call(rbind, combn(pop@id, 2, simplify = FALSE))
	colnames(parents) <- c("p1", "p2")
	pE <- getPopMeanVar(parVal, K, Vg)
	Eq <- weightedQuantile(mu = pE$pbar, sigmasq = pE$crossvar, quant = quant, w = w)
	# Eq <- w * pE$pbar + (1-w) * qnorm(quant, sd = sqrt(pE$crossvar)) # does this make sense?
	selCrit <- data.frame(parents, selCrit = Eq)

	lenSel <- 0
	while(lenSel < nCrosses / nEx){
		if(lenSel > 0) {
			cat("Not enough possible crosses with maxP =", maxP, "! Increasing maxP to", maxP + 1, "and retrying...\n")
			maxP <- maxP + 1
		}
		selection <- getSel(selCrit, n = nCrosses, high = FALSE, maxP = maxP)
		lenSel <- nrow(selection)
	}

	selection <- getSel(selCrit, n = nCrosses)
	if(nEx > 1) selection <- selection[rep(1:nrow(selection), times = nEx)[1:nCrosses], ]
	if(nProgeny > 1) selection <- selection[rep(1:nrow(selection), each = nProgeny), ] 
	makeCross(pop, crossPlan = selection) 
}

# select and cross
# expDistSel <- function(nSel, pop, GSfit, quant, nProgeny = 1, distFunc = expDist, pullGeno = pullSnpGeno, Gvar = varA, w = 0.5, ...) {
# 	expVar <- do.call(distFunc, getArgs(distFunc, pop = pop, GSfit = GSfit, quant = quant, w = w, ...))
# 	# expVar <- do.call(distFunc, getArgs(distFunc, pop = pop, GSfit = GSfit, quant = quant))
# 	selection <- getSel(expVar, n = nSel)
# 	if(is.data.frame(selection)){
# 		if(nProgeny > 1) selection <- selection[rep(1:nrow(selection), each = nProgeny), ] 
# 		selection <- makeCross(pop, crossPlan = selection) 
# 	} else {
# 		if(nProgeny > 1) selection <- rep(selection, each = nProgeny)
# 		selection <- pop[selection]
# 	}
# 	selection
# }

randomCross <- function(pop, nFam, nProgeny = 1){ # note this is just a random sampler, to illustrate how one might build a function to pick pairs. 
	allCrosses <- combn(pop@id, 2)
	resample <- if (nFam > ncol(allCrosses)) TRUE else FALSE 
	crosses <- allCrosses[, sample(1:ncol(allCrosses), nFam, replace = resample)]
	if (famSize > 1) crosses[, rep(1:nFam, each = nProgeny)]
	t(crosses)
}

# truncSel <- function(pop, nSel, traits = 1, ...) selectInd(pop, nInd = nSel, trait = traits, ...)

# something needs to be fixed here, ebv vs "ebv"
selectInd2 <- function(pop, nSel, use, trait = 1, selFunc = identity){
	if(is.character(use)) use <- match.fun(use)
	sel <- use(pop)
	names(sel) <- pop@id
	selection <- getSel(selFunc(sel), n = nSel)
	pop[selection]
}

truncSel <- function(pop, nSel, use, traits = 1, ...) do.call(selectInd2, getArgs(selectInd2, pop = pop, nSel = nSel, use = use, trait = traits, ...))

truncCross <- function(pop, nSel, nCrosses, use, nProgeny = 1, crossFunc = randomCross, traits = 1, ...) {
	if(nSel < nInd(pop)) selPop <- do.call(selectInd2, getArgs(selectInd2, pop = pop, nSel = nSel, use = use, trait = traits, ...)) else selPop <- pop
	# if(nInd(pop) * nProgeny < nCrosses)
	selection <- do.call(crossFunc, getArgs(crossFunc, pop = selPop, nFam = nCrosses, nProgeny = nProgeny, ...))
	makeCross(pop, selection)
}

checkFit <- function(pop){
	# geth2(GSfit)
	GSfit <- GSfunc(pop, traits = 1, use = "pheno", snpChip = 1, simParam = simParam)
	pop <- setEBV(pop, GSfit, simParam = simParam)
	M <- pullSnpGeno(pop)

	cat("alphaSimR RRBLUP iterations:", GSfit@iter, "\n")
	
	cat("alphaSimR RR-BLUP Vu:",  GSfit@Vu, ", Ve:", GSfit@Ve, "\n")

	require(EMMREML)
	rr <- emmreml(y = pheno(pop), X = matrix(1, nInd(pop), 1), Z = pullSnpGeno(pop), K = diag(simParam$snpChips[[1]]@nLoci))
	af <- getAF(pop)
	cat("emreml RR-BLUP Vu:",  rr$Vu, ", Ve:", rr$Ve, "\n")
	cat("Correlation of marker effect estimates:", cor(rr$u, GSfit@markerEff), "\n")
	cat("emreml RR-BLUP Vg:", rr$Vu * sum(2 * af * (1-af)), "\n")
	M %*% rr$u - ebv(pop) 

	# max(abs(rr$u - GSfit@markerEff))
	K <- vanRaden1()
	mean(diag(K))
	gblup <- emmreml(y = pheno(pop), X = matrix(1, nInd(pop), 1), Z = diag(nInd(pop)), K = K)
	cat("emreml GBLUP Vg:",  gblup$Vu, ", Ve:", gblup$Ve, "\n")

	cat("Prediction accuracy alphaSimR:",  getAcc(pop), "\n")
	cat("Prediction accuracy emmreml:",  cor(gblup$u, bv(pop)), "\n")
	cat("Correlation of emmreml and alphaSimR genetic effect estimates:",  cor(gblup$u, ebv(pop)), "\n")
	


	# all(ebv(pop) == M %*% GSfit@markerEff)
	# all(ebv(pop) == M %*% GSfit@markerEff + GSfit@fixEff[[1]] )
	# also note that they do not include the population mean!
	# Note that the blups from alphaSim are shifted... dont center genotypes? yes.
	# M %*% rr$u - gblup$u		
	# mean(gblup$u) 
	# mean(ebv(pop))
	# {gblup$u - ebv(pop)}[1]
	# gblup$betahat - GSfit@fixEff
	# max(abs(gblup$u - ebv(pop)))

}

getPopMeanVar <- function(parVal, parCov, Vg){
	if (ncol(parVal) > 1) stop("cannot use more than 1 trait...") 
	pbar <- combn(parVal, 2, mean)
	pCovar <- parCov[lower.tri(parCov)] * Vg
	pVarSum <- combn(diag(parCov) * Vg, 2, sum) 
	crossvar <- pVarSum - 2 * pCovar
	crossvar[crossvar < 0] <- 0
	list(pbar = pbar, crossvar = crossvar, pcovar = pCovar)
}

vanRaden1 <- function(M){
	Z <- scale(M, scale = FALSE)
	p <- attributes(Z)[["scaled:center"]] / 2
	ZZt <- tcrossprod(Z)
	ZZt / (2 * crossprod(p, 1-p)[[1]])
}

wrightsF <- function(M, returnNA = TRUE){
	n <- nrow(M)
	p <- colMeans(M) / 2
	d <- colSums(M == 1) / n # assumes M in 0, 1, 2
	v <- 2 * p *(1-p)
	wF <- 1 - d/v
	wF[is.na(wF)] <- if(returnNA) NA else 0
}



# # I must be crazy!
# u <- c(GSfit@markerEff)
# M <- pullGeno(pop)

# M <- M[, u !=0]
# dim(M)

# M 
# M <- t(unique(t(M))) 

# set.seed(123)

genCov <- function(M, u = NULL, absU = TRUE, sumVar = TRUE, scaleD = TRUE, inclm = TRUE){
	if(is.matrix(u)) u <- c(u)
	Z <- scale(M, scale = FALSE)
	m <- if(inclm) ncol(Z) else 1
	p <- attributes(Z)[["scaled:center"]] / 2
	v <- 2 * p * (1 - p)
	if(is.null(u) & sumVar){
		ZDZt <- tcrossprod(Z) / sum(v)
	} else {
		if(all(u == 1) & length(u) == 1) u <- rep(1, ncol(Z))
		d <- u
		if(absU) d <- abs(d)
		if(!sumVar) {
			seg <- v != 0
			d <- d[seg] / (v*m)[seg]
			Z <- Z[, seg]
		}
		if(scaleD) d <- d / mean(d)
		ZDZt <- tcrossprod(Z %*% diag(d), Z)
		if(sumVar) ZDZt <- ZDZt / sum(v)
	}
	ZDZt
}

# mean(diag(genCov(M)))
# mean(diag(genCov(M, u)))




getPopStats <- function(resultL, meanVariety = TRUE){
    VDPparam <- rlapply(resultL[["VDP"]], f = genParam, level = 2)
    RGSCparam <- lapply(resultL[["RGSC"]], genParam)
    # pop <- list(RGSC = RGSCparam, VDP = VDPparam)

    nYr <- length(resultL[["VDP"]][[1]])
    GScylcePerYr <- (length(resultL[["RGSC"]]) - 1) / nYr 
    yr <- 1:nYr
    Ryr <- yr * GScylcePerYr
    Rcyc <- c(0, 1:(GScylcePerYr * nYr))

    VgRGSC <- sapply(RGSCparam, "[[", "varG")
    # gvRGSC <- sapply(pop[["RGSC"]], function(x) mean(x$gv_a)) # this is misleading
    gvRGSC <- sapply(RGSCparam, function(x) mean(x$gv_a) + x$gv_mu) # this is correct
    gvVariety <- lapply(VDPparam[["variety"]], function(x) x$gv_a + x$gv_mu)
    SDgRGSC <- sqrt(VgRGSC)
   
    if (meanVariety){
      Yvariety <- sapply(gvVariety, mean)
      Xvariety <- Ryr[1:length(Yvariety)]
    } else {
      nVariety <- sapply(gvVariety, nrow)
      Yvariety <- unlist(gvVariety)
      Xvariety <- rep(Ryr[1:length(nVariety)], times = nVariety)
    }
    acc <- resultL$predAcc[["RGSC"]]
    theorMax <- maxBv(resultL$SP)
    return(list(SP = resultL$SP, paramL = resultL$paramL, Rcyc = Rcyc, Vg = VgRGSC, gv = gvRGSC, sd = SDgRGSC, vx = Xvariety, vy = Yvariety, RGSCyr = Ryr, acc = acc, theorMax = theorMax))
}


getYrange <- function(simR) range(c(simR$gv + simR$sd, simR$gv - simR$sd, simR$vy))

invertList <-  function(ll) {
    nms <- unique(unlist(lapply(ll, function(X) names(X))))
    ll <- lapply(ll, function(X) setNames(X[nms], nms))
    ll <- apply(do.call(rbind, ll), 2, as.list)
    lapply(ll, function(X) X[!sapply(X, is.null)])
}

plotPop <- function(simL, Rgen = RGSCgen, vLine = "none", popcol = "#000000", alpha = "0D", alphaMean = "0D", pch = 1){
	polycol <- paste0(popcol, alphaMean)
	popcol <- paste0(popcol, alpha)
	# attach(simL)
	xpoly <- c(Rgen, rev(Rgen), Rgen[1])
	ypoly <- c(simL$gv + simL$sd, rev(simL$gv - simL$sd), simL$gv[1] + simL$sd[1])

	polygon(x = xpoly, y = ypoly, col = polycol, border = NA)
	lines(x = Rgen, y = simL$gv, type = "l", col = popcol, lwd = 2)
	points(simL$vx, simL$vy, col = popcol, pch = pch)
    if (vLine == "linear") {
		abline(with(simL, lm(vy ~ vx)), col = popcol, lwd = 2)
    } else if (vLine %in% c("poly", "curve")){	    	
    	fit <- if(vLine == "poly") with(simL, lm(vy ~ poly(vx, 2))) else with(simL, loess(vy ~ vx))
		smx <- seq(min(simL$vx), max(simL$vx), by = 0.1)
		lines(smx, predict(fit, newdata = data.frame(vx = smx)), type = "l", col = popcol, lwd = 1)
  #   } else if (vLine == "curve"){	
		# smoothFit <- with(simL, loess(vy ~ vx)) # this is super annoying... why doesnt work with data = ?????
		# smx <- seq(min(simL$vx), max(simL$vx), by = 0.1)
		# lines(smx, predict(smoothFit, newdata = data.frame(vx = smx)), type = "l", col = popcol, lwd = 1)
    } else if (vLine == "connect"){
    	lines(simL$vx, simL$vy, col = popcol, lwd = 2)
    }
}


getIntensity <- function(x) {
    S <- x$vy - x$gv[x$RGSCyr]
	i <- S / x$Vg[x$RGSCyr]
	list(S = S, i = i)
}

# plotIntensity <- function(simL, x, popcol = "#000000"){

#     xInt <- getIntensity(simL)
#     lines(simL$Rcyc, simL$Vg, type = "l", lty = 1, col = popcol)
#     lines(simL$RGSCyr, xInt$S, type = "l", lty = 2, col = popcol)
#     legend("topright", legend = c("Vg", "S"), lty = c(1, 2))
# }

# popLabs = NULL; varLine = "connect"; meanVariety = TRUE; legendPos = "topleft"; plotReps = FALSE; plotVg = TRUE; plotSelInt = TRUE
simPlot <- function(popList, cols = "#000000", popLabs = NULL, varLine = "none", meanVariety = FALSE, legendPos = "topleft", plotReps = FALSE, plotVg = TRUE, plotSelInt = TRUE){
    avgInGen <- function(x) lapply(x, function(xx) tapply(xx, rep(1:(length(xx)/ nVar), each = nVar), mean))

    if (length(cols) != length(popList)) stop("cols must be same length as popList!")
    lineCol <- paste0(cols, "FF")
    ptCol <- paste0(cols, "FF")
    polyCol <- paste0(cols, "4D")

    if (is.null(popLabs)) popLabs <- names(popList)

  	nVar <- tail(with(popList[[1]][[1]][["paramL"]], nFam * famSize * cumprod(selectTrials)), 1)
    cyclePerYr <- popList[[1]][[1]][["paramL"]][["cyclePerYr"]]
    simStats <- lapply(popList, function(x) lapply(x, function(xx) xx[!names(xx) %in% c("SP", "paramL")]))
	# this is fucking hacky.................... UGH!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    if (meanVariety) {
    	for (i in 1:length(simStats)) {
    		for (j in 1:length(simStats[[i]])) simStats[[i]][[j]][c("vy", "vx")] <- avgInGen(simStats[[i]][[j]][c("vy", "vx")])
    	}
    }
    simStatsInv <- lapply(simStats, invertList)
    simReps <- rlapply(simStatsInv, level = 3, combine = rbind) 
    simAvg <- rlapply(simReps, f = colMeans, level = 2, na.rm = TRUE)

    RGSCyr <- simAvg[[1]][["RGSCyr"]]
    RGSCgen <- simAvg[[1]][["Rcyc"]]
    yr <- RGSCyr / cyclePerYr
    xlims <- range(c(0, RGSCgen))
    ylims <- range(sapply(simReps, getYrange)) * 1.1

    # plot means of RGSC and VDP output
    plot(NA, xlim = xlims, ylim = ylims, xaxt = "n", xlab = "generation", ylab = "standardized genetic value")
    axis(1, at = c(0, RGSCyr), labels = c(0, yr))

    for (i in 1:length(popList)){
    	if(plotReps) invisible(lapply(simStats[[i]], plotPop, Rgen = RGSCgen, popcol = cols[i]))
    	plotPop(simAvg[[i]], popcol = cols[i], alpha = "FF", alphaMean = "4D", Rgen = RGSCgen, vLine = varLine, pch = 16)
    }

    if (length(popList) > 1) {
        legend(legendPos, legend = popLabs, col=cols, lty = 1, lwd = 2, pch = 16)
      } else {
        legend(legendPos, legend = c("RGSC mean", expression(paste('RGSC ', sigma[g])), "Variety mean"), 
         bty = "n",
         col = c(lineCol, "black", ptCol),
         lty = c(1, 0, 0), lwd = c(2, 0, 0),
         pch = c(NA, 22, 16),
         pt.bg = c(NA, polyCol, ptCol),
        )
	}
	if(plotVg){
		ylims2 <- c(0, max(sapply(simAvg, "[[", "Vg")))
		plot(NA, xlim = xlims, ylim = ylims2, xaxt = "n", xlab = "generation", ylab = "Vg", main = "Genetic variance across generations")
	    axis(1, at = c(0, RGSCyr), labels = c(0, yr))
	    for (i in 1:length(simAvg)) {
	    	lines(simAvg[[i]][["Rcyc"]], simAvg[[i]][["Vg"]], type = "l", lwd = 2, lty = 1, col = cols[[i]])
	    }	
	    legend("topright", legend = popLabs, col=cols, lty = 1, lwd = 2, pch = 16)
	}


# check for selection intensity off by one
	if(plotSelInt){
		selInt <- lapply(simAvg, getIntensity)
		S <- lapply(selInt, "[[", "S")
		int <- lapply(selInt, "[[", "i")
		ylims3 <- range(unlist(S))
		# ylims4 <- range(unlist(int))
		ylims4 <- c(-2, max(unlist(int)))

		plot(NA, xlim = xlims, ylim = ylims3, xaxt = "n", xlab = "generation", ylab = "S", main = "Selection differential across generations")
	    axis(1, at = c(0, RGSCyr), labels = c(0, yr))
	    
	    for (i in 1:length(S)) {
	    	lines(RGSCyr, S[[i]], type = "l", lwd = 2, lty = 1, col = cols[[i]])
	    }	
	    legend("topright", legend = popLabs, col=cols, lty = 1, lwd = 2, pch = 16)

	   	plot(NA, xlim = xlims, ylim = ylims4, xaxt = "n", xlab = "generation", ylab = "S/Vg", main = "Selection intensity across generations")
	    axis(1, at = c(0, RGSCyr), labels = c(0, yr))
	    
	    for (i in 1:length(int)) {
	    	lines(RGSCyr, int[[i]], type = "l", lwd = 2, lty = 1, col = cols[[i]])
	    }
	    legend("topleft", legend = popLabs, col=cols, lty = 1, lwd = 2, pch = 16)
	}
}





asremlRfit <- function(pop, modelL){
	if (!all(names(modelL) %in% c("fixed", "random", "rcov"))) stop("please provide a list of formulas with elements 'fixed', 'random', and 'rcov' to 'modelL' argument")
	if (!all(sapply(modelL, class) == "formula")) modelL <- sapply(modelL, as.formula)

	if(grepl("giv", modelL[["random"]])){
		pullSnpGeno()
	}
} 


modGSfit <- function(GSfit, ebv){

}

