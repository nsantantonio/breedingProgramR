#' argumentChange function
#'
#' function to (do something)
#'
#' @param defaultArgs [value]
#' @param userArgs [value]
#' @return [value]
#' @details [fill in details here]
#' @examples none
#' @export
argumentChange <- function(defaultArgs, userArgs){
	defaultArgs[names(userArgs)] <- userArgs
	return(defaultArgs)
}
