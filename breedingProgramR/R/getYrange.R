#' getYrange function
#'
#' function to (do something)
#'
#' @param simR [value]
#' @return [value]
#' @details [fill in details here]
#' @examples none
#' @export
getYrange <- function(simR) { range(c(simR$gv + simR$sdRCRS, simR$gv - simR$sdRCRS, simR$vy)) }
