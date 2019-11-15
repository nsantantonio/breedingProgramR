#' getAcc function
#'
#' function to (do something)
#'
#' @param pop [value]
#' @param simParam [value]
#' @return [value]
#' @details [fill in details here]
#' @examples none
#' @export
getAcc <- function(pop, simParam) {
	popTrue <- gv(pop)
	popPred <- ebv(pop)
	if(nrow(popTrue) == 1 | var(c(popTrue)) == 0 | var(c(popPred)) == 0) NA else cor(popTrue, popPred)
}
