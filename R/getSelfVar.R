#' getSelfVar function
#'
#' function to (do something)
#'
#' @param M [value]
#' @param u [value]. Default is NULL
#' @param fdiff [value]. Default is NULL
#' @return [value]
#' @details [fill in details here]
#' @examples none
#' @export
getSelfVar <- function(M, u = NULL, fdiff = NULL) {
	if(is.null(u)) u <- rep(1, ncol(M))
	H <- M == 1
	Hu <- H %*% u^2 
	if(!is.null(fdiff)) Hu <- Hu * (1 - 2^(-fdiff))
	# if(ncol(Hu) == 1) Hu <- c(Hu)
	Hu
}
