#' truncSel function
#'
#' function to (do something)
#'
#' @param pop [value]
#' @param nSel [value]
#' @param use [value]
#' @param traits [value]. Default is 1
#' @return [value]
#' @details [fill in details here]
#' @examples none
#' @export
truncSel <- function(pop, nSel, use, traits = 1, ...) { do.call(selectInd2, getArgs(selectInd2, pop = pop, nSel = nSel, use = use, trait = traits, ...)) }
