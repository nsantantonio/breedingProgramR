#' getRRh2 function
#'
#' function to (do something)
#'
#' @param rrFit [value]
#' @return [value]
#' @details [fill in details here]
#' @examples none
#' @export
getRRh2 <- function(rrFit) { solve(pop0pred@Vu + pop0pred@Ve) %*% pop0pred@Vu }
