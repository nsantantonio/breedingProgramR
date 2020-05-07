#' solqp function
#'
#' function to (do something)
#'
#' @param pop [value]
#' @param GSfit [value]
#' @param use [value]
#' @param nCrosses [value]
#' @param simParam [value]
#' @param lambda [value]. Default is NULL
#' @param fthresh [value]. Default is NULL
#' @param gain [value]. Default is NULL
#' @param truncqp [value]. Default is NULL
#' @param allowSelf [value]. Default is FALSE
#' @param weightLoci [value]. Default is FALSE
#' @param pullGeno [value]. Default is pullSnpGeno
#' @param verbose [value]. Default is FALSE
#' @param nProgeny [value]. Default is 1
#' @param maxProp [value]. Default is 1
#' @return [value]
#' @details [fill in details here]
#' @examples none
#' @export
solqp <- function(pop, GSfit, use, nCrosses, simParam, lambda = NULL, fthresh = NULL, gain = NULL, truncqp = NULL, allowSelf = FALSE, weightLoci = FALSE, pullGeno = pullSnpGeno, verbose = FALSE, nProgeny = 1, maxProp = 1, ...){
	suppressMessages(require(LowRankQP))
	inbreedingCoef <- function(cee, Kmat) 1/2 * crossprod(cee, Kmat) %*% cee
	expectedGain <- function(cee, gebvs) crossprod(cee, gebvs)
	betterSample <- function(x, ...) x[sample(length(x), ...)]
	if(is.character(use)) use <- match.fun(use)

	n <- nInd(pop)
	
	M <- pullGeno(pop, simParam = simParam)
	K <- if(weightLoci) genCov(M, u = c(GSfit@markerEff)) else genCov(M)	
	gebvs <- use(pop)
	popTruth <- genParam(pop, simParam = simParam)

	Vg <- popTruth$varG
	if(var(gebvs) == 0){
		msg(2, "No variance in GEBVs! random intermating progressing...\n")
		selection <- randomCross(pop, nFam = nCrosses, nProgeny = nProgeny)
	} else {
		if(is.null(lambda) & is.null(fthresh) & is.null(gain) & is.null(truncqp)) stop("lambda {0, 1}, fthresh {>0}, gain {>0} or truncqp {>0} must be provided a value")
		if(!is.null(gain)) lambda <- 1 else if(is.null(lambda)) lambda <- 0:100 * 0.01 else if (lambda > 1 | lambda < 0) stop("Supplied lambda values must be between 0 and 1.")

		f <- NULL
		g <- NULL
		A <- matrix(1, nrow = 1, ncol = n)
		u <- matrix(1, ncol = 1, nrow = n) * maxProp
		b <- 1
		log <- list()
		cee <- list()
		for(k in 1:length(lambda)){
			H <- 2 * lambda[k] * K # the one half is included in the optimization. 
			d <- if(is.null(gain)) -(1 - lambda[k]) * gebvs else rep(0, length(gebvs))
			if(!is.null(gain)) {
				b <- c(1, gain)
				A <- rbind(A, c(gebvs))
			}
			log[[k]] <- capture.output({solutionqp <- suppressMessages(invisible(LowRankQP(Vmat = H, dvec = d, Amat = A, bvec = b, uvec = u, method = "LU", verbose = FALSE)))})
			cee[[k]] <- solutionqp$alpha
			f <- c(f, c(inbreedingCoef(solutionqp$alpha, K)))
			g <- c(g, c(expectedGain(solutionqp$alpha, gebvs)))
		}
		if (!is.null(fthresh) & is.null(gain)) {
			whichLambda <- if(fthresh < min(f)) which.min(f) else which(f == max(f[f <= fthresh]))
		} else if (!is.null(truncqp)) {
			zero <- 1 / (2 * nCrosses)
			nPars <- sapply(cee, function(x) sum(x >= (zero)))
			whichLambda <- if(truncqp > max(nPars)) which.max(nPars) else which(nPars == min(nPars[nPars >= truncqp]))
			if(length(whichLambda) > 1) whichLambda <- whichLambda[1]
		} else {
			whichLambda <- 1
		}
		msg(2, "lambda:", lambda[whichLambda])
		propPar <- round(cee[[whichLambda]] * 2 * nCrosses)
		rownames(propPar) <- pop@id
		pars <- rep(pop@id, times = propPar)
		if (length(unique(pars)) == 1){
			selection <- do.call(rbind, rep(list(rep(unique(pars), 2)), nCrosses))
		} else {
			parList <- list()
			index <- 1:length(pars)
			for(i in 1:min(nCrosses, floor(length(pars) / 2))) {
				p1 <- betterSample(index, 1)
				samplep2 <- index[pars[index] != pars[p1]]
				p2 <- if(allowSelf) betterSample(index[-p1], 1) else betterSample(samplep2, 1)
				if(length(p2) == 0) next
				if(is.na(pars[p1]) | is.na(pars[p1])) print(paste(pars[p1], pars[p2], sep = "  :  "))
				if(pars[p1] == pars[p2]) {msg(2, "oops\n"); break}
				crossi <- c(p1, p2)
				parList[[i]] <- pars[crossi]
				index <- index[!index %in% crossi]
			}
			selection <- do.call(rbind, parList)
		}
	}
	if(nProgeny > 1) selection <- selection[rep(1:nrow(selection), each = nProgeny), ] 
	list(pop = makeCross(pop, crossPlan = selection, simParam = simParam), lambda = lambda[whichLambda])
}
