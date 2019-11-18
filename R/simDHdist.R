#' simDHdist function
#'
#' function to (do something)
#'
#' @param nSel [value]
#' @param pop [value]
#' @param GSfit [value]
#' @param retQuant [value]. Default is FALSE
#' @param quant [value]. Default is 0.9
#' @param nDH [value]. Default is 200
#' @param weight [value]. Default is 0.5
#' @param nProgeny [value]. Default is 1
#' @param returnPop [value]. Default is TRUE
#' @param bigmem [value]. Default is FALSE
#' @return [value]
#' @details [fill in details here]
#' @examples none
#' @export
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
