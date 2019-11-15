#' rlapply function
#'
#' function to (do something)
#'
#' @param l [value]
#' @param f [value]. Default is identity
#' @param level [value]. Default is 1
#' @param combine [value]. Default is list
#' @param counter [value]. Default is 1
#' @return [value]
#' @details [fill in details here]
#' @examples none
#' @export
rlapply <- function(l, f = identity, level = 1, combine = list, counter = 1, ...){
	args <- list(...)
	if (counter < level){
		do.call(lapply, c(list(X = l, FUN = rlapply, f = f, level = level, combine = combine, counter = counter + 1), args))
	} else {
		result <- do.call(lapply, c(list(X = l, FUN = f), args))
		if (identical(combine, list)) return(result) else return(do.call(combine, result))
	}
}
