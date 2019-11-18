#' quant function
#'
#' function to (do something)
#'
#' @param x [value]
#' @param sigma [value]
#' @param w [value]
#' @return [value]
#' @details [fill in details here]
#' @examples none
#' @export
quant <- function(x, sigma, w) {w * mean(x) + (1 - w) * sigma * sd(x)}
