#' getTotalIntensity function
#'
#' function to (do something)
#'
#' @param x [value]
#' @return [value]
#' @details [fill in details here]
#' @examples none
#' @export
getTotalIntensity <- function(x) {
    S <- x$vy - x$gv[x$RCRSyr - 1]
	i <- S / x$Vg[x$RCRSyr - 1]
	list(S = S, i = i)
}
