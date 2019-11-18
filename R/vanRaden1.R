#' vanRaden1 function
#'
#' function to (do something)
#'
#' @param M [value]
#' @return [value]
#' @details [fill in details here]
#' @examples none
#' @export
vanRaden1 <- function(M){
	Z <- scale(M, scale = FALSE)
	p <- attributes(Z)[["scaled:center"]] / 2
	ZZt <- tcrossprod(Z)
	ZZt / (2 * crossprod(p, 1-p)[[1]])
}
