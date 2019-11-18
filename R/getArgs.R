#' getArgs function
#'
#' function to (do something)
#'
#' @param f [value]
#' @return [value]
#' @details [fill in details here]
#' @examples none
#' @export
getArgs <- function(f, ...) { list(...)[names(list(...)) %in% names(formals(f))] }
