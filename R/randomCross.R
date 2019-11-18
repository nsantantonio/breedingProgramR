#' randomCross function
#'
#' function to (do something)
#'
#' @param pop [value]
#' @param nFam [value]
#' @param nProgeny [value]. Default is 1
#' @return [value]
#' @details [fill in details here]
#' @examples none
#' @export
randomCross <- function(pop, nFam, nProgeny = 1){ # note this is just a random sampler, to illustrate how one might build a function to pick pairs. 
	allCrosses <- combn(pop@id, 2)
	resample <- if (nFam > ncol(allCrosses)) TRUE else FALSE 
	crosses <- allCrosses[, sample(1:ncol(allCrosses), nFam, replace = resample)]
	if (nProgeny > 1) crosses[, rep(1:nFam, each = nProgeny)]
	t(crosses)
}
