#' getF function
#'
#' function to (do something)
#'
#' @param pop [value]
#' @param pullGeno [value]. Default is pullSnpGeno
#' @return [value]
#' @details [fill in details here]
#' @examples none
#' @export
getF <- function(pop, pullGeno = pullSnpGeno) {1 + (0.5 - rowMeans(pullGeno(pop) == 1)) * 2} # not sue this is right...
