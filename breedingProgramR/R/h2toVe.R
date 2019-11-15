#' h2toVe function
#'
#' function to (do something)
#'
#' @param h2 [value]
#' @param Vg [value]. Default is 1
#' @return [value]
#' @details [fill in details here]
#' @examples none
#' @export
h2toVe <- function(h2, Vg = 1) { Vg * (1-h2) / h2 }
