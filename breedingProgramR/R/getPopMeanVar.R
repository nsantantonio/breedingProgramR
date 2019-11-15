#' getPopMeanVar function
#'
#' function to (do something)
#'
#' @param parVal [value]
#' @param parCov [value]
#' @param Vg [value]
#' @return [value]
#' @details [fill in details here]
#' @examples none
#' @export
getPopMeanVar <- function(parVal, parCov, Vg){
	if (ncol(parVal) > 1) stop("cannot use more than 1 trait...") 
	pbar <- combn(parVal, 2, mean)
	pCovar <- parCov[lower.tri(parCov)] * Vg
	pVarSum <- combn(diag(parCov) * Vg, 2, sum) 
	crossvar <- pVarSum - 2 * pCovar
	crossvar[crossvar < 0] <- 0
	list(pbar = pbar, crossvar = crossvar, pcovar = pCovar)
}
