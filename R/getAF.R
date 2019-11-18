#' getAF function
#'
#' function to (do something)
#'
#' @param pop [value]
#' @param pullGeno [value]. Default is pullSnpGeno
#' @return [value]
#' @details [fill in details here]
#' @examples none
#' @export
getAF <- function(pop, pullGeno = pullSnpGeno) { colMeans(pullGeno(pop)) / 2 }
