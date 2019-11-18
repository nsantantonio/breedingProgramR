#' mergePopsRec function
#'
#' function to (do something)
#'
#' @param x [value]
#' @return [value]
#' @details [fill in details here]
#' @examples none
#' @export
mergePopsRec <- function(popList) { mergePops(lapply(popList, function(x) {if (is.list(x)) mergePopsRec(x) else x})) }
