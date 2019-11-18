#' selectInd2 function
#'
#' function to (do something)
#'
#' @param pop [value]
#' @param nSel [value]
#' @param use [value]
#' @param trait [value]. Default is 1
#' @param selFunc [value]. Default is identity
#' @return [value]
#' @details [fill in details here]
#' @examples none
#' @export
selectInd2 <- function(pop, nSel, use, trait = 1, selFunc = identity){
	if(is.character(use)) use <- match.fun(use)
	sel <- use(pop)
	if(ncol(sel) < 1) stop("Something is wrong! I dont appear to have a phenotype/GEBV! perhaps I was never predicted?")
	names(sel) <- pop@id
	selection <- getSel(selFunc(sel), n = nSel)
	pop[selection]
}
