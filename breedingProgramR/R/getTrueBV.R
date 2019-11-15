#' getTrueBV function
#'
#' function to (do something)
#'
#' @param pop [value]
#' @param simParam [value]
#' @param trait [value]. Default is 1
#' @return [value]
#' @details [fill in details here]
#' @examples none
#' @export
getTrueBV <- function(pop, simParam, trait = 1) {
	M <- pullQtlGeno(pop)	
	u <- getTrueQTLeff(simParam)
	M %*% u
} 
