#' genCov function
#'
#' function to (do something)
#'
#' @param M [value]
#' @param u [value]. Default is NULL
#' @param absU [value]. Default is TRUE
#' @param sumVar [value]. Default is TRUE
#' @param scaleD [value]. Default is TRUE
#' @param inclm [value]. Default is TRUE
#' @param p [value]. Default is NULL
#' @return [value]
#' @details [fill in details here]
#' @examples none
#' @export
genCov <- function(M, u = NULL, absU = TRUE, sumVar = TRUE, scaleD = TRUE, inclm = TRUE, p = NULL){
	if(is.matrix(u)) u <- c(u)
	Z <- scale(M, scale = FALSE)
	m <- if(inclm) ncol(Z) else 1
	if(is.null(p)) p <- attributes(Z)[["scaled:center"]] / 2
	v <- 2 * p * (1 - p)
	if(is.null(u) & sumVar){
		ZDZt <- tcrossprod(Z) / sum(v)
	} else {
		if(all(u == 1) & length(u) == 1) u <- rep(1, ncol(Z))
		d <- u
		if(absU) d <- abs(d)
		if(!sumVar) {
			seg <- v != 0
			d <- d[seg] / (v*m)[seg]
			Z <- Z[, seg]
		}
		if(scaleD) d <- d / mean(d)
		ZDZt <- tcrossprod(Z %*% diag(d), Z)
		if(sumVar) ZDZt <- ZDZt / sum(v)
	}
	ZDZt
}
