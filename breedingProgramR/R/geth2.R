#' geth2 function
#'
#' function to (do something)
#'
#' @param pop [value]
#' @param GSfit [value]
#' @return [value]
#' @details [fill in details here]
#' @examples none
#' @export
geth2 <- function(pop, GSfit) {
	Vg <- estVg(pop, GSfit)
	Vg / sum(Vg, GSfit@Ve)
}
