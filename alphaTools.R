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

# selectGSfunc <- function(pop, simParam, snpChip = 1, xtimes = 2, ...) {
# 	n <- nInd(pop)
# 	m <- simParam$snpChips[[1]]@nLoci
# 	gsF <- if(n > m * xtimes) RRBLUP2 else RRBLUP
# 	print(list(...))
# 	do.call(gsF, getArgs(gsF, pop = pop, simParam = simParam, ...))
# }


dFSel <- function(dF, limit = 1, val = "selCrit", parentCols = c("p1", "p2"), returnPar = TRUE) {
	dFord <- list(fwd = order(dF[[val]]), rev = order(-dF[[val]]))
	dup <- duplicated(dF[[val]])
	dFord <- lapply(dFord, function(x) x[!dup])
	if(!(list({1:nrow(dF)}[!dup]) %in% dFord | list({nrow(dF):1}[!rev(dup)]) %in% dFord)) {
		stop("dF must be sorted in order to select!")
	}
	
	parMat <- as.matrix(dF[, parentCols])
	parents <- sort(unique(c(parMat)))
	parCount <- rep(0, length(parents))
	names(parCount) <- parents
	rows <- NULL

	for (i in 1:nrow(parMat)) {
		parenti <- parMat[i, ]
		if(all(parCount[parenti] < limit)) {
			parCount[parenti] <- parCount[parenti] + 1
			rows <- c(rows, i)
		}
	}
	if(returnPar) parMat[rows, , drop = FALSE] else rows
}
# dFSel(dF, maxP = 2)


getSel <- function(selCrit, n, high = TRUE, variable = "selCrit", parentCols = c("p1", "p2"), maxP = NULL) {
	len <- if(is.data.frame(selCrit)) nrow(selCrit) else length(selCrit)
	if(len < n) n <- nrow(selCrit)
	if (is.data.frame(selCrit)){
		selCrit <- selCrit[order(selCrit[[variable]], decreasing = high), ]
		if(!is.null(maxP)){
			sel <- dFSel(selCrit, limit = maxP, val = variable, parentCols = parentCols)
			sel <- sel[1:min(nrow(sel), n), , drop = FALSE] 
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
	cat("            weight parameter has value:", w, "\n")
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
simDHdist <- function(nSel, pop, GSfit, retQuant = FALSE, quant = 0.9, nDH = 200, weight = 0.5, nProgeny = 1, returnPop = TRUE, bigmem = FALSE, ...) {
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
	expVar <- if(retQuant) sapply(simDist, quantile, probs = quant) else weightedQuantile(sapply(simDist, mean), sapply(simDist, var), quant, w = weight)
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
# pop <- RGSC[[lastRGSCgen]]; GSfit <- GSmodel[[lastGSmodel]]; use = ebv; nSel = selectRGSCi; nCrosses = nNuclear; maxCrossPerParent = 1; nDH = 200; quant = xInt; retQuant = FALSE;  w = 0.5; weight = 0.5; nSimCrosses = 1
simDHdistPairs <- function(nSel, pop, GSfit, nCrosses, use, retQuant = FALSE, quant = 0.9, nDH = 200, weight = 0.5, maxCrossPerParent = 0, nSimCrosses = 1, nProgeny = 1, verbose = FALSE, ...) {
	n <- nInd(pop)
	if (n < nSel) nSel <-  n
	nCombos <- choose(nSel, 2) 
	nEx <- if(nCombos < nCrosses) ceiling(nCrosses / nCombos) else 1 
	maxP <- if(maxCrossPerParent == 0 | nCombos <  nCrosses) nCrosses else maxCrossPerParent

	if(nSel < n) pop <- truncSel(pop, nSel = nSel, use = use)

	parents <- do.call(rbind, combn(pop@id, 2, simplify = FALSE))
	colnames(parents) <- c("p1", "p2")
	crosses <- rep(1:nrow(parents), each = nSimCrosses)
	popX <- makeCross(pop, parents[crosses, , drop = FALSE])

	if(verbose) cat("simulating distribution of", nDH, "DH for", nSimCrosses, "crosses for each of", nCombos, "parental pairs\n")
	simVar <- do.call(simDHdist, getArgs(simDHdist, nSel = nInd(popX), pop = popX, GSfit = GSfit, retQuant = retQuant, quant = quant, nDH = nDH, weight = weight, returnPop = FALSE, ...))
	if(nSimCrosses > 1) simVar <- tapply(simVar, crosses, mean)

	selCrit <- data.frame(parents, selCrit = simVar)
	lenSel <- 0
	while(lenSel < nCrosses / nEx){
		if(lenSel > 0) {
			cat("Not enough possible crosses with maxP =", maxP, "! Increasing maxP to", maxP + 1, "and retrying...\n")
			maxP <- maxP + 1
		}
		selection <- getSel(selCrit, n = nCrosses, high = TRUE, maxP = maxP)
		lenSel <- nrow(selection)
	}
	if(verbose) print(table(selection))

	if(nEx > 1) selection <- selection[rep(1:nrow(selection), times = nEx)[1:nCrosses], ]
	if(nProgeny > 1) selection <- selection[rep(1:nrow(selection), each = nProgeny), ] 
	makeCross(pop, crossPlan = selection) 
}

# select individuals
expDist <- function(nSel, pop, GSfit, use, quant, returnQuant = TRUE, pullGeno = pullSnpGeno, updateEBV = FALSE, weight = 0.5, nProgeny = 1, ...){
	expVar <- do.call(getSelfVar, getArgs(getSelfVar, M = pullGeno(pop), u = GSfit@markerEff, ...))
	if(returnQuant) {
		if(updateEBV) pop <- setEBV(pop, GSfit)
		parVal <- ebv(pop)
		expVar <- weightedQuantile(mu = parVal, sigmasq = expVar, quant = quant, w = weight)
	} 
	if(ncol(expVar) == 1) expVar <- expVar[, 1]
	selection <- getSel(expVar, n = nSel, high = TRUE)
	if(nProgeny > 1) selection <- rep(selection, each = nProgeny)
	pop[selection]
}

# pop = RGSC[[lastRGSCgen]]; GSfit = GSmodel[[lastGSmodel]]; nSel = selectRGSCi; nCrosses = nNuclear; use = ebv; pullGeno = pullSnpGeno; weightLoci = FALSE; maxCrossPerParent = 1; 
maxVar <- function(pop, GSfit, nSel, nCrosses, use, weightLoci = FALSE, pullGeno = pullSnpGeno, maxCrossPerParent = 0, verbose = FALSE, nProgeny = 1, ...){
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
# use = ebv; nSel = selectRGSCi; pop = RGSC[[lastRGSCgen]]; GSfit = GSmodel[[lastGSmodel]]; quant = xInt; nCrosses = nNuclear
# returnQuant = TRUE; weightLoci = FALSE; pullGeno = pullSnpGeno; Gvar = estVg; w = 0.5; nProgeny = 1
expDistPairs <- function(pop, GSfit, nSel, quant, nCrosses, use, returnQuant = TRUE, weightLoci = FALSE, pullGeno = pullSnpGeno, maxCrossPerParent = 0, Gvar = estVg, weight = 0.5, nProgeny = 1, verbose = FALSE, ...) {
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

	Vg <- do.call(Gvar, getArgs(Gvar, pop = pop, GSfit = GSfit, ...))
	if (prod(dim(Vg)) > 1) stop("can only handle a single trait!") else Vg <- Vg[[1]]
	parents <- do.call(rbind, combn(pop@id, 2, simplify = FALSE))
	colnames(parents) <- c("p1", "p2")
	pE <- getPopMeanVar(parVal, K, Vg)
	Eq <- weightedQuantile(mu = pE$pbar, sigmasq = pE$crossvar, quant = quant, w = weight)
	# Eq <- w * pE$pbar + (1-w) * qnorm(quant, sd = sqrt(pE$crossvar)) # does this make sense?
	selCrit <- data.frame(parents, selCrit = Eq)

	lenSel <- 0
	while(lenSel < nCrosses / nEx){
		if(lenSel > 0) {
			cat("Not enough possible crosses with maxP =", maxP, "! Increasing maxP to", maxP + 1, "and retrying...\n")
			maxP <- maxP + 1
		}
		# selection <- getSel(selCrit, n = nCrosses, high = FALSE, maxP = maxP) # WAIT!!!!!! should this be low???????????
		selection <- getSel(selCrit, n = nCrosses, high = TRUE, maxP = maxP) # rerun with high = TRUE??
		lenSel <- nrow(selection)
	}
	if(verbose) print(table(selection))
	
	selection <- getSel(selCrit, n = nCrosses)
	if(nEx > 1) selection <- selection[rep(1:nrow(selection), times = nEx)[1:nCrosses], ]
	if(nProgeny > 1) selection <- selection[rep(1:nrow(selection), each = nProgeny), ] 
	makeCross(pop, crossPlan = selection)
}
##############

#############

pherFuncID <- function(x, p) {x}
pherFuncMax <- function(x, p) {(x / max(x))^p}
quant <- function(x, sigma, w) {w * mean(x) + (1 - w) * sigma * sd(x)}
popQuant <- function(sel, ebvs, sigma = 2, w = 0.5) {quant(x = ebvs[sel], sigma = sigma, w = 0.5)}
initFunc <- function(x) {rep(1, length(x))} # can change later 
gsQuant <- function(sel, ebvs, sigma = 2, w = 0.5, nCrosses = 100, nProgenyPerCross = 1) {
	if(!all(ebv(pop[sel]) == ebvs[sel])) stop("ebvs of pop dont match those of ant!")
	simCrosses <- randomCross(pop[sel], nFam = nCrosses, nProgeny = nProgenyPerCross)
	simPop <- makeCross(pop[sel], simCrosses)
	simPop <- setEBV(simPop, GSfit)
	quant(ebv(simPop), sigma = sigma, w = w)
}

acOpt <- function(x, n, targetFunc, pherFunc, xAt0 = FALSE, evapRate = 0.05, nAnts = 500, maxIter = 300, countThresh = 1000, showResult = FALSE, pherPower = 1, dumbAnt = FALSE, plateau = 100, returnStats = FALSE){
	evap <- function(oldP, newP, evapRate) {(1 - evapRate) * oldP + newP }

	if(is.matrix(x)) {if(ncol(x) == 1) x <- c(x) else stop("I cant handle multiple traits yet!")}
	
	if(xAt0){
		minx <- min(x)
		x <- x - minx
	}
	N <- length(x)
	pheromone <- initFunc(x)

	bestAnt <- NULL
	noP = rep(0, N)
	bestestPath = 0
	lastBestPath <- bestestPath
	pathCounter = 0 
	iter = 0
	if(showResult) {
		plotBest <- NULL 
		plotMean <- NULL 
	}

	while(pathCounter <= countThresh & iter < maxIter){
		iter = iter + 1
		ant <- list()
		path <- list()
		antPheromone <- list()
		for(i in 1:nAnts){
			ant[[i]] <- sample(1:N, n, prob = pheromone / sum(pheromone))
			path[[i]] <- targetFunc(sel = ant[[i]], ebvs = x)
			# path[[i]] <- quant(x[ant[[i]]], sigma = 2, w = 0.5)
			Pi <- noP
			Pi[ant[[i]]] <- path[[i]]
			antPheromone[[i]] <- Pi
		}
		bestPath <- which.max(path)
		if (path[[bestPath]] == lastBestPath) {
			pathCounter <- pathCounter + 1
		} else {
			lastBestPath <- path[[bestPath]]
			pathCounter <- 0
		}
		if(path[[bestPath]] >= bestestPath) {
			bestestPath <- path[[bestPath]] 
			bestAnt <- ant[[bestPath]]
		}
		newPheromone <- pherFunc(Reduce("+", antPheromone), p = pherPower)	
		pheromone <- evap(pheromone, newPheromone, evapRate = evapRate)
	}
	return(bestAnt)
}

# pop = RGSC[[lastRGSCgen]]; GSfit = GSmodel[[lastGSmodel]]; acTrunc = 1; evapRate = 0.05; nAnts = 500; pherPower = 1.5; nSel = selectRGSCi; nCrosses = nNuclear; use = ebv; pullGeno = pullSnpGeno; weightLoci = FALSE; maxCrossPerParent = 1; nProgeny = 1
ACquant <- function(pop, GSfit, nSel, nCrosses, use, acTrunc = 1, evapRate = 0.05, nAnts = 500, pherPower = 1.5, verbose = FALSE, nProgeny = 1, ...){
	n <- nInd(pop)
	if (n < nSel) nSel <-  n
	nCombos <- choose(nSel, 2)
	nEx <- if(nCombos < nCrosses) ceiling(nCrosses / nCombos) else 1 
	# maxP <- if(maxCrossPerParent == 0 | nCombos <  nCrosses) nCrosses else maxCrossPerParent
	popVar <- varA(pop)
	popMean <- mean(gv(pop))

	if (acTrunc < 1) pop <- truncSel(pop, nSel = n * acTrunc, use = use)
	if (n == nSel) {
		selectedParents <- 1:n
	} else {
		selectedParents <- acOpt(ebv(pop), n = nSel, xAt0 = TRUE, targetFunc = popQuant, pherFunc = pherFuncMax, evapRate = evapRate, nAnts = nAnts, pherPower = pherPower)
	}
	selection <- randomCross(pop[selectedParents], nFam = nCrosses, nProgeny = nProgeny)

	if(nEx > 1) selection <- selection[rep(1:nrow(selection), times = nEx)[1:nCrosses], ]
	if(nProgeny > 1) selection <- selection[rep(1:nrow(selection), each = nProgeny), ] 
	newpop <- makeCross(pop, crossPlan = selection) 
	if(verbose){
		cat("Selected population variance diff:", {varA(newpop) - popVar}, "\n")
		cat("Selected population mean diff:", {mean(gv(newpop)) - popMean}, "\n")
	}
	newpop
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
	if (nProgeny > 1) crosses[, rep(1:nFam, each = nProgeny)]
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


getTotalIntensity <- function(x) {
    S <- x$vy - x$gv[x$RGSCyr - 1]
	i <- S / x$Vg[x$RGSCyr - 1]
	list(S = S, i = i)
}


# note:
# all(abs((x$bv + x$mu) -  (x$gv_a + x$gv_mu)) < 1e-14)

# resultL <- simrun[[1]]
getPopStats <- function(resultL, meanVariety = TRUE, verbose = FALSE){
    VDPparam <- rlapply(resultL[["VDP"]], f = genParam, level = 2)
    RGSCparam <- lapply(resultL[["RGSC"]], genParam)
    # pop <- list(RGSC = RGSCparam, VDP = VDPparam)

    nYr <- length(resultL[["VDP"]][[1]])
    GScylcePerYr <- (length(resultL[["RGSC"]]) - 1) / nYr 
    yr <- 1:nYr
    Ryr <- yr * GScylcePerYr
    Rcyc <- c(0, 1:(GScylcePerYr * nYr))

    VgRGSC <- sapply(RGSCparam, "[[", "varG")
    # gvRGSC <- sapply(pop[["RGSC"]], function(x) mean(x$gv_a)) # this is misleading. Why?
    gvRGSC <- sapply(RGSCparam, function(x) mean(x$gv_a) + x$gv_mu) # this is correct
    sRGSC <- gvRGSC[-1] - gvRGSC[-length(gvRGSC)]
    iRGSC <- sRGSC / sqrt(VgRGSC[-length(VgRGSC)])

    VgVDP <- rlapply(VDPparam, "[[", i = "varG", level = 2, combine = c)
    gvVDP <- rlapply(VDPparam, function(x) mean(x$gv_a) + x$gv_mu, level = 2, combine = c)
    sVDP <- gvVDP$variety - gvVDP$trial1
    iVDP <- sVDP / sqrt(VgVDP$trial1)
  
  	RyrIndex <- Ryr - GScylcePerYr + 1
  	sTotal <- gvVDP$variety - gvRGSC[RyrIndex]
	iTotal <- sTotal / sqrt(VgRGSC[RyrIndex])

    gvVariety <- lapply(VDPparam[["variety"]], function(x) x$gv_a + x$gv_mu)
    SDgRGSC <- sqrt(VgRGSC)

    varMean <- gvVDP[["variety"]]
	nVariety <- sapply(gvVariety, nrow)
	Yvariety <- unlist(gvVariety)
	Xvariety <- rep(Ryr[1:length(nVariety)], times = nVariety)

    RGSCacc <- resultL$predAcc[["RGSC"]]
    # add VDP pred acc here! for skipping gens
    # if(any(names(resultL$predAcc)) VDPacc <- 
    theorMax <- maxBv(resultL$SP)
	
	nVar = unique(nVariety)

    return(list(SP = resultL$SP, paramL = resultL$paramL, Rcyc = Rcyc, varMean = varMean, sdRGSC = SDgRGSC, 
				VgRGSC = VgRGSC, VgVDP = VgVDP, gvRGSC = gvRGSC, gvVDP = gvVDP,
    			sRGSC = sRGSC, iRGSC = iRGSC, sVDP = sVDP, iVDP = iVDP, sTotal = sTotal, iTotal = iTotal, 
    			nVar = nVar, vx = Xvariety, vy = Yvariety, RGSCyr = Ryr, RGSCacc = RGSCacc, theorMax = theorMax))
}


getYrange <- function(simR) range(c(simR$gv + simR$sdRGSC, simR$gv - simR$sdRGSC, simR$vy))

invertList <-  function(ll) {
    nms <- unique(unlist(lapply(ll, function(X) names(X))))
    ll <- lapply(ll, function(X) setNames(X[nms], nms))
    ll <- apply(do.call(rbind, ll), 2, as.list)
    lapply(ll, function(X) X[!sapply(X, is.null)])
}

plotPop <- function(simL, Rgen = RGSCgen, vLine = "none", popcol = "#000000", alpha = "0D", alphaMean = "0D", pch = 1) {
	polycol <- paste0(popcol, alphaMean)
	popcol <- paste0(popcol, alpha)
	# attach(simL)
	xpoly <- c(Rgen, rev(Rgen), Rgen[1])
	ypoly <- c(simL$gvRGSC + simL$sdRGSC, rev(simL$gvRGSC - simL$sdRGSC), simL$gvRGSC[1] + simL$sdRGSC[1])

	polygon(x = xpoly, y = ypoly, col = polycol, border = NA)
	lines(x = Rgen, y = simL$gvRGSC, type = "l", col = popcol, lwd = 2)
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


# getTotalIntensity <- function(x) {
#     S <- x$vy - x$gv[x$RGSCyr - 1]
# 	i <- S / x$Vg[x$RGSCyr - 1]
# 	list(S = S, i = i)
# }

# plotIntensity <- function(simL, x, popcol = "#000000"){

#     xInt <- getIntensity(simL)
#     lines(simL$Rcyc, simL$Vg, type = "l", lty = 1, col = popcol)
#     lines(simL$RGSCyr, xInt$S, type = "l", lty = 2, col = popcol)
#     legend("topright", legend = c("Vg", "S"), lty = c(1, 2))
# }

# popList <- list(p1 = simrun, p2 = simrun)

# popLabs = NULL; varLine = "connect"; meanVariety = TRUE; legendPos = "topleft"; plotReps = FALSE; plotVg = TRUE; plotSelInt = TRUE
simPlot <- function(popList, cols = "#000000", popLabs = NULL, varLine = "none", meanVariety = TRUE, legendPos = "topleft", plotReps = FALSE, plotVg = TRUE, plotSelInt = TRUE){

	# cols = c("#000000", "#000000")
    if (length(cols) != length(popList)) stop("cols must be same length as popList!")
    lineCol <- paste0(cols, "FF")
    ptCol <- paste0(cols, "FF")
    polyCol <- paste0(cols, "4D")

    if (is.null(popLabs)) popLabs <- names(popList)

  	# nVar <- tail(with(popList[[1]][[1]][["paramL"]], nFam * famSize * cumprod(selectTrials)), 1)
  	nVar <- unique(unlist(rlapply(popList, "[[", i = "nVar", level = 2, combine = c)))
    cyclePerYr <- unique(unlist(rlapply(popList, function(x) x[["paramL"]][["cyclePerYr"]], level = 2, combine = c)))
    if(length(nVar) > 1) warning("number of varieties differ between pops!")
    if(length(nVar) > 1) warning("cycles per year differ between pops!")

    simStats <- lapply(popList, function(x) lapply(x, function(xx) xx[!names(xx) %in% c("SP", "paramL", "VgVDP", "gvVDP")]))

    simStatsInv <- lapply(simStats, invertList)
    simReps <- rlapply(simStatsInv, level = 3, combine = rbind) 
    simAvg <- rlapply(simReps, f = colMeans, level = 2, na.rm = TRUE)

    RGSCyr <- simAvg[[1]][["RGSCyr"]]
    RGSCgen <- simAvg[[1]][["Rcyc"]]
    yr <- RGSCyr / cyclePerYr
    xlims <- range(c(0, RGSCgen))
    ylims <- range(sapply(simReps, getYrange)) * 1.1

    if (meanVariety) {
    	simAvg <- lapply(simAvg, function(x) {x[["vy"]] <- x[["varMean"]]; x[["vx"]] <- x[["RGSCyr"]]; x})
    }
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
		ylims2 <- c(0, max(sapply(simAvg, "[[", "VgRGSC")))
		plot(NA, xlim = xlims, ylim = ylims2, xaxt = "n", xlab = "generation", ylab = "Vg", main = "Genetic variance across generations")
	    axis(1, at = c(0, RGSCyr), labels = c(0, yr))
	    for (i in 1:length(simAvg)) {
	    	lines(simAvg[[i]][["Rcyc"]], simAvg[[i]][["Vg"]], type = "l", lwd = 2, lty = 1, col = cols[[i]])
	    }	
	    legend("topright", legend = popLabs, col=cols, lty = 1, lwd = 2, pch = 16)
	}



# THIS NEEDS TO BE CLEANED UP!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	if(plotSelInt){
		# for(j in c("sVDP", "iVDP", "sRGSC", "iRGSC"){
		# 	si <- lapply(simAvg, "[[", j)
		# 	ylimsSi <- range(c(0, unlist(si)))
	
		# 	plot(NA, xlim = xlims, ylim = ylims3, xaxt = "n", xlab = "generation", ylab = "S", main = "Selection differential across generations")
		#     axis(1, at = c(0, RGSCyr), labels = c(0, yr))
		    
		#     for (i in 1:length(sRGSC)) {
		#     	lines(RGSCyr, sVDP[[i]], type = "l", lwd = 2, lty = 1, col = cols[[i]])
		#     }	
		#     legend("topright", legend = popLabs, col=cols, lty = 1, lwd = 2, pch = 16)


		# }


		# selInt <- lapply(simAvg, getIntensity)
		sVDP <- lapply(simAvg, "[[", "sVDP")
		iVDP <- lapply(simAvg, "[[", "iVDP")
		sRGSC <- lapply(simAvg, "[[", "sRGSC")
		iRGSC <- lapply(simAvg, "[[", "iRGSC")

		ylims3 <- range(unlist(sVDP))
		# ylims4 <- range(unlist(int))
		ylims4 <- c(-2, max(unlist(iVDP)))

		ylims5 <- range(unlist(sRGSC))
		# ylims4 <- range(unlist(int))
		ylims6 <- c(-2, max(unlist(iRGSC)))


		plot(NA, xlim = xlims, ylim = ylims3, xaxt = "n", xlab = "generation", ylab = "S", main = "Selection differential across generations in VDP")
	    axis(1, at = c(0, RGSCyr), labels = c(0, yr))
	    
	    for (i in 1:length(sRGSC)) {
	    	lines(RGSCyr, sVDP[[i]], type = "l", lwd = 2, lty = 1, col = cols[[i]])
	    }	
	    legend("topright", legend = popLabs, col=cols, lty = 1, lwd = 2, pch = 16)

	   	plot(NA, xlim = xlims, ylim = ylims4, xaxt = "n", xlab = "generation", ylab = "S/Vg", main = "Selection intensity across generations in VDP")
	    axis(1, at = c(0, RGSCyr), labels = c(0, yr))
	    
	    for (i in 1:length(iVDP)) {
	    	lines(RGSCyr, iVDP[[i]], type = "l", lwd = 2, lty = 1, col = cols[[i]])
	    }

	    legend("topleft", legend = popLabs, col=cols, lty = 1, lwd = 2, pch = 16)
		plot(NA, xlim = xlims, ylim = ylims5, xaxt = "n", xlab = "generation", ylab = "S", main = "Selection differential across generations in RGSC")
	    axis(1, at = c(0, RGSCyr), labels = c(0, yr))
	    
	    for (i in 1:length(sRGSC)) {
	    	lines(RGSCgen[-1], sRGSC[[i]], type = "l", lwd = 2, lty = 1, col = cols[[i]])
	    }	
	    legend("topright", legend = popLabs, col=cols, lty = 1, lwd = 2, pch = 16)

	   	plot(NA, xlim = xlims, ylim = ylims6, xaxt = "n", xlab = "generation", ylab = "S/Vg", main = "Selection intensity across generations in RGSC")
	    axis(1, at = c(0, RGSCyr), labels = c(0, yr))
	    
	    for (i in 1:length(iRGSC)) {
	    	lines(RGSCgen[-1], iRGSC[[i]], type = "l", lwd = 2, lty = 1, col = cols[[i]])
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

