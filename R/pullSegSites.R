#' pullSegSites function
#'
#' function to (do something)
#'
#' @param pop [value]
#' @param returnMatrix [value]. Default is TRUE
#' @return [value]
#' @details [fill in details here]
#' @examples none
#' @export
pullSegSites <- function(pop, returnMatrix = TRUE){
	rawToSum <- function(xk) {
		class(xk) <- "integer"
		rowSums(xk)
	}
	geno <- lapply(pop@geno, function(x) t(apply(x, 3, rawToSum)))
	if (returnMatrix) geno <- do.call(cbind, geno)
	geno
}
