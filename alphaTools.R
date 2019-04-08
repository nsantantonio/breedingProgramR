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
simDHdist <- function(pop, GSfit, quant = 0.9, DHsampleSize = 200) {
	qDH <- function(i){
		DH <- makeDH(pop[i], nDH = DHsampleSize)
		DH <- setEBV(DH, GSfit, simParam = simParam)
		DHebv <- ebv(DH)
		quantile(DHebv, quant)[[1]]
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

getExpDist <- function(pop, GSfit, quant, pullGeno = pullSnpGeno, Gvar = varA) {
	parVal <- ebv(pop)
	rownames(parVal) <- pop@id
	K <- vanRaden1(pullGeno(pop))
	Vg <- Gvar(pop)
	if(prod(dim(Vg)) > 1) stop("can only handle a single trait!") else Vg <- Vg[[1]]
	pE <- getPopMeanVar(parVal, K, Vg)
	do.call(rbind, combn(pop@id, 2, simplify = FALSE))
	pE$pbar + qnorm(quant, sd = sqrt(pE$crossvar)) # does this make sense?
	# cor(pE$pbar, pE$pbar + qnorm(1 - eSelInt, sd = sqrt(pE$crossvar)))
}


crossPlanFunc <- function(pop, nFam, famSize = 1){ # note this is just a random sampler, to illustrate how one might build a function to pick pairs. 
	allCrosses <- combn(pop@id, 2)
	resample <- if(nFam > ncol(allCrosses)) TRUE else FALSE 
	crosses <- allCrosses[, sample(1:ncol(allCrosses), nFam, replace = resample)]
	if(famSize > 1) crosses[, rep(1:nFam, each = famSize)]
	t(crosses)
}



getPopStats <- function(resultL, meanVariety = TRUE){
    VDPparam <- rlapply(resultL[["VDP"]], f = genParam, level = 2)
    RGSCparam <- lapply(resultL[["RGSC"]], genParam)
    # pop <- list(RGSC = RGSCparam, VDP = VDPparam)

    nYr <- length(resultL[["VDP"]][[1]])
    GScylcePerYr <- (length(resultL[["RGSC"]]) - 1) / nYr 
    yr <- 1:nYr
    Ryr <- yr * GScylcePerYr
    Rcyc <- 1:(GScylcePerYr * nYr)

    VgRGSC <- sapply(RGSCparam, "[[", "varG")
    # gvRGSC <- sapply(pop[["RGSC"]], function(x) mean(x$gv_a)) # this is misleading
    gvRGSC <- sapply(RGSCparam, function(x) mean(x$gv_a) + x$gv_mu) # this is correct
    gvVariety <- lapply(VDPparam[["variety"]], function(x) x$gv_a + x$gv_mu)
    SDgRGSC <- sqrt(VgRGSC)
   
    if(meanVariety){
      Yvariety <- sapply(gvVariety, mean)
      Xvariety <- Ryr[1:length(Yvariety)]
    } else {
      nVariety <- sapply(gvVariety, nrow)
      Yvariety <- unlist(gvVariety)
      Xvariety <- rep(Ryr[1:length(nVariety)], times = nVariety)
    }
    acc <- resultL$predAcc[["RGSC"]]
    theorMax <- maxBv(resultL$SP)
    return(list(SP = resultL$SP, Rcyc = Rcyc, Vg = VgRGSC, gv = gvRGSC, sd = SDgRGSC, vx = Xvariety, vy = Yvariety, RGSCyr = Ryr, acc = acc, theorMax = theorMax))
}



# getPopStats <- function(pop, Ryr = RGSCyr, meanVariety = TRUE){
#     VgRGSC <- sapply(pop[["RGSC"]], genicVarG)
#     gvRGSC <- sapply(pop[["RGSC"]], function(x) mean(gv(x)))
#     # gvVariety <- sapply(run1[["VDP"]][["variety"]], function(x) mean(gv(x)))
#     gvVariety <- lapply(pop[["VDP"]][["variety"]], function(x) gv(x))
#     SDgRGSC <- sqrt(VgRGSC)
   
#     # meanVariety <- TRUE
#     if(meanVariety){
#       Yvariety <- sapply(gvVariety, mean)
#       Xvariety <- Ryr[1:length(Yvariety)]
#     } else {
#       nVariety <- sapply(gvVariety, nrow)
#       Yvariety <- unlist(gvVariety)
#       Xvariety <- rep(Ryr[1:length(nVariety)], times = nVariety)
#     }

#     return(list(Vg = VgRGSC, gv = gvRGSC, sd = SDgRGSC, vx = Xvariety, vy = Yvariety))
# }

getYrange <- function(simR) range(c(simR$gv + simR$sd, simR$gv - simR$sd, simR$vy))

invertList <-  function(ll) {
    nms <- unique(unlist(lapply(ll, function(X) names(X))))
    ll <- lapply(ll, function(X) setNames(X[nms], nms))
    ll <- apply(do.call(rbind, ll), 2, as.list)
    lapply(ll, function(X) X[!sapply(X, is.null)])
}

plotPop <- function(simL, Rgen = RGSCgen, vLine = FALSE, popcol = "#000000", alpha = "0D", alphaMean = "0D"){
    polycol <- paste0(popcol, alphaMean)
    popcol <- paste0(popcol, alpha)
    # attach(simL)
    xpoly <- c(Rgen, rev(Rgen), Rgen[1])
    ypoly <- c(simL$gv + simL$sd, rev(simL$gv - simL$sd), simL$gv[1] + simL$sd[1])

    polygon(x = xpoly, y = ypoly, col = polycol, border = NA)
    lines(x = Rgen, y = simL$gv, type = "l", col = popcol, lwd = 2)
    points(simL$vx, simL$vy, col = popcol, pch = 1)
    if(vLine) abline(lm(simL$vy ~ simL$vx ), col = popcol)
  }

simPlot <- function(popList, baseCols = "#000000", popLabs = NULL, varLine = TRUE){
    if(class(popList[[1]][[1]][[1]]) == "Pop") popList <- list(popList) # checks that popList is a list of pop
    
    if(length(baseCols) != length(popList)) stop("baseCols must be same length as popList!")
    lineCol <- paste0(baseCols, "FF")
    ptCol <- paste0(baseCols, "FF")
    polyCol <- paste0(baseCols, "4D")

    if(!is.null(popLabs)) popLabs <- names(popList)
    yr <- 1:nYr
    RGSCyr <- yr * GScylcePerYr
    RGSCgen <- c(0, 1:(nYr * GScylcePerYr))
    
    simStats <- rlapply(popList, f=getPopStats, level = 2, Ryr = RGSCyr)
    simStatsInv <- lapply(simStats, invertList)

    simReps <- rlapply(simStatsInv, level = 3, combine = rbind) 
    simAvg <- rlapply(simReps, f = colMeans, level = 2)

    xlims <- range(c(0, RGSCgen))
    ylims <- range(sapply(simReps, getYrange)) * 1.1

    plot(NA, xlim = xlims, ylim = ylims, xaxt = "n", xlab = "generation", ylab = "standardized genetic value")
    axis(1, at = c(0, RGSCyr), labels = c(0, yr))

    for(i in 1:length(popList)){
      invisible(lapply(simStats[[i]], plotPop, Rgen = RGSCgen))
      plotPop(simAvg[[i]], popcol = baseCols[i], alpha = "FF", alphaMean = "4D", Rgen = RGSCgen, vLine = varLine)
    }

    if(length(popList) > 1) {
        legend("topright", legend = popLabs, col=baseCols, lty = 1)
      } else {
        legend("topright", legend = c("RGSC mean", expression(paste('RGSC ', sigma[g])), "Variety mean"), 
         # fill = c(NA, polyCol, NA), 
         bty = "n",
         col = c(lineCol, "black", ptCol),
         lty = c(1, 0, 0), lwd = c(2, 0, 0),
         pch = c(NA, 22, 16),
         pt.bg = c(NA, polyCol, ptCol),
        )
  }
}

# # popList <- list(simrun, simrun)
# popList <- list(simrun)




# sel <- if(RGSCselect) "sel" else "nosel" 
# f <- if(selF2) "F2" else "F1" 

# source("testparam.R")
# # pdf("RGSCselF2.pdf", width = 10, height = 7)
# pdf(paste0("RGSC", sel, f, ".pdf"), width = 10, height = 7)
# # simPlot(popList, c("#00FF00", "#0000FF"))
# simPlot(popList)
# dev.off()

