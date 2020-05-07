#' invertList function
#'
#' function to (do something)
#'
#' @param ll [value]
#' @return [value]
#' @details [fill in details here]
#' @examples none
#' @export
invertList <-  function(ll) {
	allIsNull <- function(x) all(sapply(x, is.null))
	if(allIsNull(ll)) return(NULL)
    nms <- unique(unlist(lapply(ll, function(X) names(X))))
    ll <- lapply(ll, function(X) setNames(X[nms], nms))
    ll <- apply(do.call(rbind, ll), 2, as.list)
    lapply(ll, function(X) X[!sapply(X, is.null)])
}
