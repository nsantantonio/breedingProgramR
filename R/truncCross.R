#' truncCross function
#'
#' function to (do something)
#'
#' @param pop [value]
#' @param nSel [value]
#' @param nCrosses [value]
#' @param use [value]
#' @param nProgeny [value]. Default is 1
#' @param crossFunc [value]. Default is randomCross
#' @param traits [value]. Default is 1
#' @return [value]
#' @details [fill in details here]
#' @examples none
#' @export
truncCross <- function(pop, nSel, nCrosses, use, nProgeny = 1, crossFunc = randomCross, traits = 1, ...) {
	if(nSel < nInd(pop)) selPop <- do.call(selectInd2, getArgs(selectInd2, pop = pop, nSel = nSel, use = use, trait = traits, ...)) else selPop <- pop
	selection <- do.call(crossFunc, getArgs(crossFunc, pop = selPop, nFam = nCrosses, nProgeny = nProgeny, ...))
	makeCross(pop, selection)
}
