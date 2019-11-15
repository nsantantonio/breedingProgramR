#' expDist function
#'
#' function to (do something)
#'
#' @param nSel [value]
#' @param pop [value]
#' @param GSfit [value]
#' @param use [value]
#' @param quant [value]
#' @param returnQuant [value]. Default is TRUE
#' @param pullGeno [value]. Default is pullSnpGeno
#' @param updateEBV [value]. Default is FALSE
#' @param weight [value]. Default is 0.5
#' @param nProgeny [value]. Default is 1
#' @return [value]
#' @details [fill in details here]
#' @examples none
#' @export
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
