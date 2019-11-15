#' maxVar function
#'
#' function to (do something)
#'
#' @param pop [value]
#' @param GSfit [value]
#' @param nSel [value]
#' @param nCrosses [value]
#' @param use [value]
#' @param weightLoci [value]. Default is FALSE
#' @param pullGeno [value]. Default is pullSnpGeno
#' @param maxCrossPerParent [value]. Default is 0
#' @param verbose [value]. Default is FALSE
#' @param nProgeny [value]. Default is 1
#' @return [value]
#' @details [fill in details here]
#' @examples none
#' @export
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
			msg(2, "Not enough possible crosses with maxP =", maxP, "! Increasing maxP to", maxP + 1, "and retrying...\n")
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
