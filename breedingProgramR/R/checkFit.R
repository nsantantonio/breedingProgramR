#' checkFit function
#'
#' function to (do something)
#'
#' @param pop [value]
#' @return [value]
#' @details [fill in details here]
#' @examples none
#' @export
checkFit <- function(pop){
	GSfit <- GSfunc(pop, traits = 1, use = "pheno", snpChip = 1, simParam = simParam)
	pop <- setEBV(pop, GSfit, simParam = simParam)
	M <- pullSnpGeno(pop)

	msg(2, "alphaSimR RRBLUP iterations:", GSfit@iter)
	
	msg(2, "alphaSimR RR-BLUP Vu:",  GSfit@Vu, ", Ve:", GSfit@Ve)

	require(EMMREML)
	rr <- emmreml(y = pheno(pop), X = matrix(1, nInd(pop), 1), Z = pullSnpGeno(pop), K = diag(simParam$snpChips[[1]]@nLoci))
	af <- getAF(pop)
	msg(2, "emreml RR-BLUP Vu:",  rr$Vu, ", Ve:", rr$Ve)
	msg(2, "Correlation of marker effect estimates:", cor(rr$u, GSfit@markerEff))
	msg(2, "emreml RR-BLUP Vg:", rr$Vu * sum(2 * af * (1-af)))
	M %*% rr$u - ebv(pop) 

	K <- vanRaden1()
	mean(diag(K))
	gblup <- emmreml(y = pheno(pop), X = matrix(1, nInd(pop), 1), Z = diag(nInd(pop)), K = K)
	msg(2, "emreml GBLUP Vg:",  gblup$Vu, ", Ve:", gblup$Ve)

	msg(2, "Prediction accuracy alphaSimR:",  getAcc(pop))
	msg(2, "Prediction accuracy emmreml:",  cor(gblup$u, bv(pop)))
	msg(2, "Correlation of emmreml and alphaSimR genetic effect estimates:",  cor(gblup$u, ebv(pop)))
}
