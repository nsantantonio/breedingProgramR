#' pullLoci function
#'
#' function to (do something)
#'
#' @param simParam [value]
#' @param snpChip [value]. Default is 1
#' @param asList [value]. Default is FALSE
#' @return [value]
#' @details [fill in details here]
#' @examples none
#' @export
pullLoci <- function(simParam, snpChip = 1, asList = FALSE) {
	splitIndex <- function(lpc, loc) split(loc, rep(1:length(lpc), times = lpc))
	getIndex <- function(sites, nPerChr, loc){
		indexList <- splitIndex(sites, 1:sum(sites))
		locList <- splitIndex(nPerChr, loc)
		unlist(lapply(1:length(locList), function(i) indexList[[i]][locList[[i]]]))
	}
	sites <- SP$segSites
	QTLnPerChr <- simParam$traits[[snpChip]]@lociPerChr
	SNPnPerChr <- simParam$snpChips[[snpChip]]@lociPerChr
	QTLloc <- simParam$traits[[snpChip]]@lociLoc 
	SNPloc <- simParam$snpChips[[snpChip]]@lociLoc
	QTLsites <- getIndex(sites, QTLnPerChr, QTLloc)
	SNPsites <- getIndex(sites, SNPnPerChr, SNPloc)
	intersect(QTLsites, SNPsites)
	list(QTL = QTLsites, SNP = SNPsites)
}
