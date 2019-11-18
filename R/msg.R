#' msg function
#'
#' function to (do something)
#'
#' @param n [value]
#' @return [value]
#' @details [fill in details here]
#' @examples none
#' @export
msg <- function(n, ...) { cat(do.call(paste, c(rep(list("    "), n), ..., collapse = "", sep = " ")), "\n", sep = "") }
