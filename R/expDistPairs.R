#' expDistPairs function
#'
#' function to (do something)
#'
#' @param pop [value]
#' @param GSfit [value]
#' @param nSel [value]
#' @param quant [value]
#' @param nCrosses [value]
#' @param use [value]
#' @param returnQuant [value]. Default is TRUE
#' @param weightLoci [value]. Default is FALSE
#' @param pullGeno [value]. Default is pullSnpGeno
#' @param maxCrossPerParent [value]. Default is 0
#' @param Gvar [value]. Default is estVg
#' @param weight [value]. Default is 0.5
#' @param nProgeny [value]. Default is 1
#' @param verbose [value]. Default is FALSE
#' @return [value]
#' @details [fill in details here]
#' @examples none
#' @export
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
			msg(2, "Not enough possible crosses with maxP =", maxP, "! Increasing maxP to", maxP + 1, "and retrying...\n")
			maxP <- maxP + 1
		}
		selection <- getSel(selCrit, n = nCrosses, high = TRUE, maxP = maxP) # rerun with high = TRUE??
		lenSel <- nrow(selection)
	}
	if(verbose) print(table(selection))
	
	selection <- getSel(selCrit, n = nCrosses)
	if(nEx > 1) selection <- selection[rep(1:nrow(selection), times = nEx)[1:nCrosses], ]
	if(nProgeny > 1) selection <- selection[rep(1:nrow(selection), each = nProgeny), ] 
	makeCross(pop, crossPlan = selection)
}
