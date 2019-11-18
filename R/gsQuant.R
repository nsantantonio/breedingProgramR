#' gsQuant function
#'
#' function to (do something)
#'
#' @param sel [value]
#' @param ebvs [value]
#' @param sigma [value]. Default is 2
#' @param w [value]. Default is 0.5
#' @param nCrosses [value]. Default is 100
#' @param nProgenyPerCross [value]. Default is 1
#' @return [value]
#' @details [fill in details here]
#' @examples none
#' @export
gsQuant <- function(sel, ebvs, sigma = 2, w = 0.5, nCrosses = 100, nProgenyPerCross = 1) {
	if(!all(ebv(pop[sel]) == ebvs[sel])) stop("ebvs of pop dont match those of ant!")
	simCrosses <- randomCross(pop[sel], nFam = nCrosses, nProgeny = nProgenyPerCross)
	simPop <- makeCross(pop[sel], simCrosses)
	simPop <- setEBV(simPop, GSfit)
	quant(ebv(simPop), sigma = sigma, w = w)
}
