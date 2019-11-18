#' maxBv function
#'
#' function to (do something)
#'
#' @param x [value]
#' @return [value]
#' @details [fill in details here]
#' @examples none
#' @export
maxBv <- function(simParam, traits = 1) { sapply(simParam$traits[traits], function(x) { sum(abs(x@addEff)) }) }
