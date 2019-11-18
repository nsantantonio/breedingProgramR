#' sampleFounderPop function
#'
#' function to (do something)
#'
#' @param founderPop [value]
#' @param size [value]
#' @param n [value]
#' @param savePop [value]. Default is FALSE
#' @param seed [value]. Default is NULL
#' @return [value]
#' @details [fill in details here]
#' @examples none
#' @export
sampleFounderPop <- function(founderPop, size, n, savePop = FALSE, seed = NULL) {
	if(!is.null(seed)) set.seed(seed)
	popSamples <- list()
	for(i in 1:n){
		samplei <- sample(1:nInd(founderPop), nFounder)
		popSamples[[i]] <- if(savePop) founderPop[samplei] else samplei 
	}
	popSamples
}
