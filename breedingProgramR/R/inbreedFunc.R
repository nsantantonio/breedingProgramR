#' inbreedFunc function
#'
#' function to (do something)
#'
#' @param pop [value]
#' @param GSfit [value]
#' @param nProgPerFam [value]
#' @return [value]
#' @details [fill in details here]
#' @examples none
#' @export
inbreedFunc <- function(pop, GSfit, nProgPerFam, ...) {
	pop <- if(ssh) self(pop, nProgeny = nProgPerFam) else makeDH(pop, nDH = nProgPerFam)
	if (withinFamInt < 1) {
		pop <- setEBV(pop, GSmodel[[lastGSmodel]], simParam = simParam)
		fams <- split(1:(nProgPerFam * nFam), rep(1:nFam, each = nProgPerFam))
		faml <- list()
		for (j in names(fams)){
			faml[[j]] <- selectInd(pop[fams[[j]]], nInd = famSize, trait = 1, use = selectOut) 
		}
		pop <- mergePops(faml)		
	}
	# # print mean genotypic value of DH 
	# if (verbose) print(sapply(VDP[[trials[1]]], function(x) mean(gv(x))))
	if(nInd(pop) != nFam * famSize) stop("selToP is wrong...")
}
