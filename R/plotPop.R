#' plotPop function
#'
#' function to (do something)
#'
#' @param simL [value]
#' @param vLine [value]. Default is "none"
#' @param popcol [value]. Default is "#000000"
#' @param alpha [value]. Default is "0D"
#' @param alphaMean [value]. Default is "0D"
#' @param pch [value]. Default is 1
#' @return [value]
#' @details [fill in details here]
#' @examples none
#' @export
plotPop <- function(simL, vLine = "none", popcol = "#000000", alpha = "80", alphaMean = "0D", ...) {
	polycolor <- paste0(popcol, alphaMean)
	popcolor <- paste0(popcol, alpha)
	Rgen <- simL$Rcyc / simL$paramL$cyclePerYr
	xpoly <- c(Rgen, rev(Rgen), Rgen[1])
	ypoly <- c(simL$gvRCRS + simL$sdRCRS, rev(simL$gvRCRS - simL$sdRCRS), simL$gvRCRS[1] + simL$sdRCRS[1])

	# plot(NA, ylim = c(0, max(simL$vy)), xlim = c(0, simL$paramL$nYr))

	polygon(x = xpoly, y = ypoly, col = polycolor, border = NA)
	lines(x = Rgen, y = simL$gvRCRS, type = "l", col = popcolor, lwd = 2)
	points(simL$vx, simL$vy, col = popcolor, ...)
    if (vLine == "linear") {
		abline(with(simL, lm(vy ~ vx)), col = popcolor, lwd = 2)
    } else if (vLine %in% c("poly", "curve")){	    	
    	fit <- if(vLine == "poly") with(simL, lm(vy ~ poly(vx, 2))) else with(simL, loess(vy ~ vx))
		smx <- seq(min(simL$vx), max(simL$vx), by = 0.1)
		lines(smx, predict(fit, newdata = data.frame(vx = smx)), type = "l", col = popcolor, lwd = 2)
    } else if (vLine == "connect"){
    	lines(1:simL$paramL$nYr, simL$gvVDP$variety, col = popcolor, lwd = 2)
    }
}



# traditionalTruncation$Rcyc
# rapidCycleTruncation$Rcyc

# traditionalTruncation$RCRSyr
# rapidCycleTruncation$RCRSyr