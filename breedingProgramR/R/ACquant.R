#' ACquant function
#'
#' function to (do something)
#'
#' @param pop [value]
#' @param GSfit [value]
#' @param nSel [value]
#' @param nCrosses [value]
#' @param use [value]
#' @param acTrunc [value]. Default is 1
#' @param evapRate [value]. Default is 0.05
#' @param nAnts [value]. Default is 500
#' @param pherPower [value]. Default is 1.5
#' @param verbose [value]. Default is FALSE
#' @param nProgeny [value]. Default is 1
#' @return [value]
#' @details [fill in details here]
#' @examples none
#' @export
ACquant <- function(pop, GSfit, nSel, nCrosses, use, acTrunc = 1, evapRate = 0.05, nAnts = 500, pherPower = 1.5, verbose = FALSE, nProgeny = 1, ...){
	if(is.character(use)) use <- match.fun(use)
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
		selectedParents <- acOpt(use(pop), n = nSel, xAt0 = TRUE, targetFunc = popQuant, pherFunc = pherFuncMax, evapRate = evapRate, nAnts = nAnts, pherPower = pherPower)
	}
	selection <- randomCross(pop[selectedParents], nFam = nCrosses, nProgeny = nProgeny)

	if(nEx > 1) selection <- selection[rep(1:nrow(selection), times = nEx)[1:nCrosses], ]
	if(nProgeny > 1) selection <- selection[rep(1:nrow(selection), each = nProgeny), ] 
	newpop <- makeCross(pop, crossPlan = selection) 
	if(verbose){
		msg(2, "Selected population variance diff:", {varA(newpop) - popVar})
		msg(2, "Selected population mean diff:", {mean(gv(newpop)) - popMean})
	}
	newpop
}
