#' estVg function
#'
#' function to (do something)
#'
#' @param pop [value]
#' @param GSfit [value]. Default is NULL
#' @param GSfunc [value]. Default is RRBLUP
#' @return [value]
#' @details [fill in details here]
#' @examples none
#' @export
estVg <- function(pop, GSfit = NULL, GSfunc = RRBLUP) {
	if(is.null(GSfit)) GSfit <- GSfunc(pop)
	af <- getAF(pop)
	Vg <- GSfit@Vu
	Vg <- if(is.matrix(Vg) & prod(dim(Vg)) == 1) Vg[[1]]
	Vg * sum(2 * af * (1-af))
}
