#' plotPop function
#'
#' function to (do something)
#'
#' @param simL [value]
#' @param Rgen [value]. Default is RCRSgen
#' @param vLine [value]. Default is "none"
#' @param popcol [value]. Default is "#000000"
#' @param alpha [value]. Default is "0D"
#' @param alphaMean [value]. Default is "0D"
#' @param pch [value]. Default is 1
#' @return [value]
#' @details [fill in details here]
#' @examples none
#' @export
plotPop <- function(simL, Rgen = RCRSgen, vLine = "none", popcol = "#000000", alpha = "0D", alphaMean = "0D", pch = 1) {
	polycol <- paste0(popcol, alphaMean)
	popcol <- paste0(popcol, alpha)
	xpoly <- c(Rgen, rev(Rgen), Rgen[1])
	ypoly <- c(simL$gvRCRS + simL$sdRCRS, rev(simL$gvRCRS - simL$sdRCRS), simL$gvRCRS[1] + simL$sdRCRS[1])

	polygon(x = xpoly, y = ypoly, col = polycol, border = NA)
	lines(x = Rgen, y = simL$gvRCRS, type = "l", col = popcol, lwd = 2)
	points(simL$vx, simL$vy, col = popcol, pch = pch)
    if (vLine == "linear") {
		abline(with(simL, lm(vy ~ vx)), col = popcol, lwd = 2)
    } else if (vLine %in% c("poly", "curve")){	    	
    	fit <- if(vLine == "poly") with(simL, lm(vy ~ poly(vx, 2))) else with(simL, loess(vy ~ vx))
		smx <- seq(min(simL$vx), max(simL$vx), by = 0.1)
		lines(smx, predict(fit, newdata = data.frame(vx = smx)), type = "l", col = popcol, lwd = 1)
    } else if (vLine == "connect"){
    	lines(simL$vx, simL$vy, col = popcol, lwd = 2)
    }
}
