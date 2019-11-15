#' withinFamSel function
#'
#' function to (do something)
#'
#' @param pop [value]
#' @param GSfit [value]
#' @param use [value]
#' @param int [value]
#' @param nProgeny [value]
#' @param nGenInbr [value]
#' @param nGenThisYr [value]
#' @param lambdaInbr [value]. Default is NULL
#' @param fthreshInbr [value]. Default is NULL
#' @param gainInbr [value]. Default is NULL
#' @param truncqpInbr [value]. Default is NULL
#' @param ssd [value]
#' @param simParam [value]
#' @return [value]
#' @details [fill in details here]
#' @examples none
#' @export
withinFamSel <- function(pop, GSfit, use, int, nProgeny, nGenInbr, nGenThisYr, lambdaInbr = NULL, fthreshInbr = NULL, gainInbr = NULL, truncqpInbr = NULL, ssd, simParam, ...){
	nFam <- nInd(pop)
	nProgPerFam <- nProgeny / int
	if (nProgPerFam %% 1 != 0) {
		nProgPerFam <- round(nProgPerFam)
		nSelToTrial <- round(nProgPerFam * withinFamInt)
		msg(2, "\nNOTE: Selection intensities within familiy have been rounded to the nearest integer resulting in", nSelToTrial, "progeny per family selected from", nProgPerFam, "progeny per family\n")
	}
	if (ssd) {
		i <- 0
		while(i <= nGenInbr) {
			nP <- if(i == 0) nProgPerFam else 1 
			pop <- self(pop, nProgeny = nP)
			i <- i + 1
		}
	} else {
		pop <- makeDH(pop, nDH = nProgPerFam)
	}

	pop <- setEBV(pop, GSfit, simParam = simParam)
	fams <- split(1:(nProgeny * nFam), rep(1:nFam, each = nProgeny))
	faml <- list()
	for (j in names(fams)){
		faml[[j]] <- selectInd(pop[fams[[j]]], nInd = nProgeny, trait = 1, use = use) 
	}
	mergePops(faml)
}
