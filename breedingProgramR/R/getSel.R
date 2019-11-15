#' getSel function
#'
#' function to (do something)
#'
#' @param selCrit [value]
#' @param n [value]
#' @param high [value]. Default is TRUE
#' @param variable [value]. Default is "selCrit"
#' @param parentCols [value]. Default is c("p1", "p2")
#' @param maxP [value]. Default is NULL
#' @return [value]
#' @details [fill in details here]
#' @examples none
#' @export
getSel <- function(selCrit, n, high = TRUE, variable = "selCrit", parentCols = c("p1", "p2"), maxP = NULL) {
	len <- if(is.data.frame(selCrit)) nrow(selCrit) else length(selCrit)
	if(len < n) n <- nrow(selCrit)
	if (is.data.frame(selCrit)){
		selCrit <- selCrit[order(selCrit[[variable]], decreasing = high), ]
		if(!is.null(maxP)){
			sel <- dFSel(selCrit, limit = maxP, val = variable, parentCols = parentCols)
			sel <- sel[1:min(nrow(sel), n), , drop = FALSE] 
		} else {
			sel <- as.matrix(selCrit[1:n, parentCols])
		}
	} else {
		sel <- names(sort(selCrit, decreasing = high))[1:n]
	}
	sel
}
