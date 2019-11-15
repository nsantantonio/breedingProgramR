#' weightedQuantile function
#'
#' function to (do something)
#'
#' @param mu [value]
#' @param sigmasq [value]
#' @param quant [value]. Default is 0.9
#' @param w [value]. Default is 0.5
#' @return [value]
#' @details [fill in details here]
#' @examples none
#' @export
weightedQuantile <- function(mu, sigmasq, quant = 0.9, w = 0.5) {
	msg(2, "            weight parameter has value:", w)
	w * mu + (1-w) * qnorm(quant, sd = sqrt(sigmasq))
}
