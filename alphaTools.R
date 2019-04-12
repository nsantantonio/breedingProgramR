# getArgs <- function(defaultArgs) {
#   args <- commandArgs(TRUE)
#   isAssn <- grepl("=", args)
#   userArgs <- args[isAssn]
#   needEval <- grepl("\\(|\\)|\\:", userArgs) 
#   argSplit <- strsplit(userArgs, "=")
#   argList <- lapply(argSplit, "[[", 2)
#   names(argList) <- lapply(argSplit, "[[", 1)
#   argList[needEval] <- lapply(argList[needEval], function(x) eval(parse(text = x)))
#   argList[!needEval] <- lapply(argList[!needEval], function(x) strsplit(x, ",")[[1]])
#   argList[!needEval] <- type.convert(argList[!needEval], as.is = TRUE)
#   print(argList)
#   defaultArgs[names(argList)] <- argList
#   defaultArgs
# }

getComArgs <- function(defaultArgs = NULL) {
  defaults <- !is.null(defaultArgs)
  args <- commandArgs(TRUE)
  isAssn <- grepl("=", args)
  userArgs <- args[isAssn]
  needEval <- grepl("\\(|\\)|\\:", userArgs) 
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

mergePopsRec <- function(popList) mergePops(lapply(popList, function(x) if (is.list(x)) mergePopsRec(x) else x))


h2toVe <- function(h2, Vg = 1) Vg * (1-h2) / h2
geth2 <- function(GSfit) GSfit@Vu / sum(GSfit@Vu, GSfit@Ve)
gen <- function(i) paste0("gen", i)
logSelInd <- function(pop, sel) pop@id %in% sel  
rSel <- function(sel) Reduce("&", sel) 
maxBv <- function(simParam, traits = 1) sapply(simParam$traits[traits], function(x) sum(abs(x@addEff)))
sdUnCor <- function(x) sqrt(mean(x^2) - mean(x)^2)
getRRh2 <- function(rrFit) solve(pop0pred@Vu + pop0pred@Ve) %*% pop0pred@Vu
getAcc <- function(pop) cor(gv(pop), ebv(pop))
dummyFunc <- function(x, retrn) retrn

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

estIntensity <- function(VDP, i, start = "trial1", end = "variety", estFunc = pheno, Gvar = varA) {
	S <- mean(pheno(VDP[[end]][[gen(i - nTrial)]])) - mean(pheno(VDP[[start]][[gen(i - nTrial)]]))
	i <- S / sqrt(Gvar(VDP[[start]][[gen(i - nTrial)]]))
}



getSel <- function(selCrit, n) {
	if (is.data.frame(selCrit)){
		selCrit <- selCrit[order(selCrit[["selCrit"]], decreasing = TRUE), ]
		sel <- as.matrix(selCrit[1:n, c("p1", "p2")])
	} else {
		sel <- names(sort(selCrit, decreasing = TRUE))[1:n]
	}
	sel
}

getTrueQTLeff <- function(simParam, trait = 1) simParam$traits[[trait]]@addEff 

getTrueBV <- function(pop, simParam, trait = 1) {
	M <- pullQtlGeno(pop)	
	u <- getTrueQTLeff(simParam)
	M %*% u
} 




# 	selCrit <- selFuncOut(RGSC[[lastRGSCgen]], GSmodel[[lastGSmodel]], quant = xInt[[gen(i)]])
# 	parSel <- getSel(selCrit, nFam)

# 	if (is.matrix(parSel)){
# 		if (ncol(parSel) == 2){
# 			# if (famSize > 1) parSel <- parSel[rep(1:nFam, each = famSize), ] # dont think this is necessary, you dont need to make the cropss more than once, just DH/self. 
# 			selGStoP <- makeCross(RGSC[[lastRGSCgen]], crossPlan = parSel) # perhaps this should be set to the previous RGSC generation??
# 		} else {
# 			stop("Something is wrong with parent selection. Expecting 2 columns, p1 and p2 indicating parent pairs.")
# 		}
# 	} else {
# 		selGStoP <- RGSC[[lastRGSCgen]][parSel]
# 		if (any(!selGStoP@id %in% parSel)) stop("parent selection out of RGSC failed!")
# 	}
# } else {
# }

# what do you want to do here?

getSelfVar <- function(M, u, fdiff = NULL) {
	H <- M == 1
	Hu <- H %*% u^2 
	if(!is.null(fdiff)) Hu <- Hu * (1 - 2^(-fdiff))
	Hu
}

# function to return expected quantiles from sampling DH individuals
simDHdist <- function(pop, GSfit, retQuant = FALSE, quant = 0.9, n = 200) {
	DHdist <- function(i){
		DH <- makeDH(pop[i], nDH = n)
		DH <- setEBV(DH, GSfit, simParam = simParam)
		DHebv <- ebv(DH)
		if(retQuant) quantile(DHebv, quant)[[1]] else var(DHebv)
	}
	sapply(pop@id, DHdist)
}

expDist <- function(pop, GSfit, quant, returnQuant = TRUE, sampleDH = FALSE, pullGeno = pullSnpGeno, Gvar = varA, updateEBV = FALSE, ...){
		expVar <- do.call(getSelfVar, getArgs(getSelfVar, M = pullGeno(pop), u = GSfit@markerEff, ...))
		if(returnQuant) {
			if(updateEBV) pop <- setEBV(pop, GSfit)
			parVal <- ebv(pop)
			x <- qnorm(quant, sd = sqrt(expVar))
			expVar <- parVal + x * expVar
		} # need top check why this seems to work so poorly...
	}
	expVar
}

expDistSel <- function(nSel, pop, GSfit, quant, distFunc = expDistInd, pullGeno = pullSnpGeno, Gvar = varA, ...) {
	# expVar <- do.call(distFunc, getArgs(distFunc, pop = pop, GSfit = GSfit, quant = quant, ...))
	expVar <- distFunc(pop = pop, GSfit = GSfit, quant = quant, ...))
	selection <- getSel(expVar, n = nSel)
	if(is.data.frame(selection)){
		selection <- makeCross(pop, crossPlan = selection) 
	}
	selection
}


expDistPairs <- function(pop, GSfit, quant, returnQuant = TRUE, pullGeno = pullSnpGeno, Gvar = varA) {
	parVal <- ebv(pop)
	rownames(parVal) <- pop@id
	M <- pullGeno(pop)
	K <- genCov() else K <- genCov(M, u = c(GSfit@markerEff))

	Vg <- Gvar(pop)
	if (prod(dim(Vg)) > 1) stop("can only handle a single trait!") else Vg <- Vg[[1]]
	parents <- do.call(rbind, combn(pop@id, 2, simplify = FALSE))
	colnames(parents) <- c("p1", "p2")
	pE <- getPopMeanVar(parVal, K, Vg)
	Eq <- pE$pbar + qnorm(quant, sd = sqrt(pE$crossvar)) # does this make sense?
	selCrit <- data.frame(parents, selCrit = Eq)
  # cor(pE$pbar, pE$pbar + qnorm(1 - eSelInt, sd = sqrt(pE$crossvar)))
}

truncSel <- function(nSel, pop, traits = 1, ...) selectInd(pop, nInd = nSel, trait = traits, use = selectOutRGSC, ...)

# 	# rownames(parVal) <- pop@id
# 	# M <- pullGeno(pop)
# 	# rownames(M) <- pop@id

# 	# expVar <- expDistInd(pop, GSfit, 0.9)
# 	expVar <- expDistInd(pop, GSfit, 0.9, returnQuant = FALSE)
# 	expVarDH <- expDistInd(pop, GSfit, 0.9, returnQuant = FALSE, sampleDH = TRUE)
# 	# expVarDH <- simDHdist(pop, GSfit, retQuant = FALSE)
# 	trueVar <- getSelfVar(pullQtlGeno(pop), getTrueQTLeff(SP), fdiff = NULL)

# 	parVal <- ebv(pop)


# 	cor(expVar, expVarDH)
# 	cor(parVal, expVar)
# 	cor(parVal, expVarDH)
# 	cor(trueVar, expVar)
# 	cor(trueVar, expVarDH) # not sampling DH is better than expected?

# 	K <- genCov(M)
# 	K <- genCov(M, u = c(u))
# 	K <- genCov(M, u = c(u), scaleD = TRUE, sumVar = FALSE, printMeanDiag = TRUE)

# 	Vg <- Gvar(pop)
# 	if (prod(dim(Vg)) > 1) stop("can only handle a single trait!") else Vg <- Vg[[1]]
# 	parents <- do.call(rbind, combn(pop@id, 2, simplify = FALSE))
# 	colnames(parents) <- c("p1", "p2")
# 	pE <- getPopMeanVar(parVal, K, Vg)
# 	Eq <- pE$pbar + qnorm(quant, sd = sqrt(pE$crossvar)) # does this make sense?
# 	selCrit <- data.frame(parents, selCrit = Eq)

# 	parSel <- getSel(selCrit, nFam)

#   # cor(pE$pbar, pE$pbar + qnorm(1 - eSelInt, sd = sqrt(pE$crossvar)))
# }




getPopMeanVar <- function(parVal, parCov, Vg){
	if (ncol(parVal) > 1) stop("cannot use more than 1 trait...") 
	pbar <- combn(parVal, 2, mean)
	pCovar <- parCov[lower.tri(parCov)] * Vg
	pVarSum <- combn(diag(parCov) * Vg, 2, sum) 
	crossvar <- max(0, pVarSum - 2 * pCovar)
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



randCross <- function(pop, nFam, famSize = 1){ # note this is just a random sampler, to illustrate how one might build a function to pick pairs. 
	allCrosses <- combn(pop@id, 2)
	resample <- if (nFam > ncol(allCrosses)) TRUE else FALSE 
	crosses <- allCrosses[, sample(1:ncol(allCrosses), nFam, replace = resample)]
	if (famSize > 1) crosses[, rep(1:nFam, each = famSize)]
	t(crosses)
}

getPopStats <- function(resultL, meanVariety = FALSE){
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
    } else if (vLine == "curve"){	
		smoothFit <- with(simL, loess(vy ~ vx)) # this is super annoying... why doesnt work with data = ?????
		smx <- seq(min(simL$vx), max(simL$vx), by = 0.1)
		lines(smx, predict(smoothFit, newdata = data.frame(vx = smx)), type = "l", col = popcol, lwd = 1)
    } else if (vLine == "connect"){
    	lines(simL$vx, simL$vy, col = popcol, lwd = 2)
    }
}

simPlot <- function(popList, baseCols = "#000000", popLabs = NULL, varLine = "none", meanVariety = TRUE, legendPos = "topleft"){
    avgInGen <- function(x, lx) lapply(x, function(xx, len = lx) tapply(xx, rep(1:(length(xx)/len), each = len), mean))

    if (length(baseCols) != length(popList)) stop("baseCols must be same length as popList!")
    lineCol <- paste0(baseCols, "FF")
    ptCol <- paste0(baseCols, "FF")
    polyCol <- paste0(baseCols, "4D")

    if (is.null(popLabs)) popLabs <- names(popList)

    GScylcePerYr <- popList[[1]][[1]][["paramL"]][["GScylcePerYr"]]
    simStats <- lapply(popList, function(x) lapply(x, function(xx) xx[!names(xx) %in% c("SP", "paramL")]))
	# this is fucking hacky.................... UGH!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    if (meanVariety) {
    	for (i in 1:length(simStats)) {
    		l <- length(simStats[[i]])
    		for (j in 1:l) simStats[[i]][[j]][c("vy", "vx")] <- avgInGen(simStats[[i]][[j]][c("vy", "vx")], l)
    	}
    }
    simStatsInv <- lapply(simStats, invertList)
    simReps <- rlapply(simStatsInv, level = 3, combine = rbind) 
    simAvg <- rlapply(simReps, f = colMeans, level = 2, na.rm = TRUE)

    RGSCyr <- simAvg[[1]][["RGSCyr"]]
    RGSCgen <- simAvg[[1]][["Rcyc"]]
    yr <- RGSCyr / GScylcePerYr
    xlims <- range(c(0, RGSCgen))
    ylims <- range(sapply(simReps, getYrange)) * 1.1

    plot(NA, xlim = xlims, ylim = ylims, xaxt = "n", xlab = "generation", ylab = "standardized genetic value")
    axis(1, at = c(0, RGSCyr), labels = c(0, yr))

    for (i in 1:length(popList)){
    	invisible(lapply(simStats[[i]], plotPop, Rgen = RGSCgen, popcol = baseCols[i]))
    	plotPop(simAvg[[i]], popcol = baseCols[i], alpha = "FF", alphaMean = "4D", Rgen = RGSCgen, vLine = varLine, pch = 16)
    }

    if (length(popList) > 1) {
        legend(legendPos, legend = popLabs, col=baseCols, lty = 1, lwd = 2)
      } else {
        legend(legendPos, legend = c("RGSC mean", expression(paste('RGSC ', sigma[g])), "Variety mean"), 
         bty = "n",
         col = c(lineCol, "black", ptCol),
         lty = c(1, 0, 0), lwd = c(2, 0, 0),
         pch = c(NA, 22, 16),
         pt.bg = c(NA, polyCol, ptCol),
        )
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

