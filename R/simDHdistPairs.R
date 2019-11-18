#' simDHdistPairs function
#'
#' function to (do something)
#'
#' @param nSel [value]
#' @param pop [value]
#' @param GSfit [value]
#' @param nCrosses [value]
#' @param use [value]
#' @param retQuant [value]. Default is FALSE
#' @param quant [value]. Default is 0.9
#' @param nDH [value]. Default is 200
#' @param weight [value]. Default is 0.5
#' @param maxCrossPerParent [value]. Default is 0
#' @param nSimCrosses [value]. Default is 1
#' @param nProgeny [value]. Default is 1
#' @param verbose [value]. Default is FALSE
#' @return [value]
#' @details [fill in details here]
#' @examples none
#' @export
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

	if(verbose) msg(2, "simulating distribution of", nDH, "DH for", nSimCrosses, "crosses for each of", nCombos, "parental pairs\n")
	simVar <- do.call(simDHdist, getArgs(simDHdist, nSel = nInd(popX), pop = popX, GSfit = GSfit, retQuant = retQuant, quant = quant, nDH = nDH, weight = weight, returnPop = FALSE, ...))
	if(nSimCrosses > 1) simVar <- tapply(simVar, crosses, mean)

	selCrit <- data.frame(parents, selCrit = simVar)
	lenSel <- 0
	while(lenSel < nCrosses / nEx){
		if(lenSel > 0) {
			msg(2, "Not enough possible crosses with maxP =", maxP, "! Increasing maxP to", maxP + 1, "and retrying...\n")
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
