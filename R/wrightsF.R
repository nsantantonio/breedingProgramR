#' wrightsF function
#'
#' function to (do something)
#'
#' @param M [value]
#' @param returnNA [value]. Default is TRUE
#' @return [value]
#' @details [fill in details here]
#' @examples none
#' @export
wrightsF <- function(M, returnNA = TRUE){
	n <- nrow(M)
	p <- colMeans(M) / 2
	d <- colSums(M == 1) / n # assumes M in 0, 1, 2
	v <- 2 * p *(1-p)
	wF <- 1 - d/v
	wF[is.na(wF)] <- if(returnNA) NA else 0
}
