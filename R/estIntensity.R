#' estIntensity function
#'
#' function to (do something)
#'
#' @param VDP [value]
#' @param i [value]
#' @param nT [value]. Default is nTrial
#' @param start [value]. Default is "trial1"
#' @param end [value]. Default is "variety"
#' @param estFunc [value]. Default is pheno
#' @param Gvar [value]. Default is estVg
#' @return [value]
#' @details [fill in details here]
#' @examples none
#' @export
estIntensity <- function(VDP, i, nT = nTrial, start = "trial1", end = "variety", estFunc = pheno, Gvar = estVg, simParam = NULL) {
	S <- mean(pheno(VDP[[end]][[gen(i - nT)]])) - mean(pheno(VDP[[start]][[gen(i - nT)]]))
	i <- S / sqrt(Gvar(VDP[[start]][[gen(i - nT)]], simParam = simParam))
	if(is.matrix(i) & prod(dim(i)) == 1) i <- i[[1]] else if(is.matrix(i)) msg(2, "intensity has dimensions:", dim(i))
	i
}
