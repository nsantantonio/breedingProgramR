#' dFSel function
#'
#' function to (do something)
#'
#' @param dF [value]
#' @param limit [value]. Default is 1
#' @param val [value]. Default is "selCrit"
#' @param parentCols [value]. Default is c("p1", "p2")
#' @param returnPar [value]. Default is TRUE
#' @return [value]
#' @details [fill in details here]
#' @examples none
#' @export
dFSel <- function(dF, limit = 1, val = "selCrit", parentCols = c("p1", "p2"), returnPar = TRUE) {
	dFord <- list(fwd = order(dF[[val]]), rev = order(-dF[[val]]))
	dup <- duplicated(dF[[val]])
	dFord <- lapply(dFord, function(x) x[!dup])
	if(!(list({1:nrow(dF)}[!dup]) %in% dFord | list({nrow(dF):1}[!rev(dup)]) %in% dFord)) {
		stop("dF must be sorted in order to select!")
	}
	
	parMat <- as.matrix(dF[, parentCols])
	parents <- sort(unique(c(parMat)))
	parCount <- rep(0, length(parents))
	names(parCount) <- parents
	rows <- NULL

	for (i in 1:nrow(parMat)) {
		parenti <- parMat[i, ]
		if(all(parCount[parenti] < limit)) {
			parCount[parenti] <- parCount[parenti] + 1
			rows <- c(rows, i)
		}
	}
	if(returnPar) parMat[rows, , drop = FALSE] else rows
}
