#' simPlot function
#'
#' function to (do something)
#'
#' @param popList [value]
#' @param cols [value]. Default is "#000000"
#' @param popLabs [value]. Default is NULL
#' @param varLine [value]. Default is "none"
#' @param meanVariety [value]. Default is TRUE
#' @param legendPos [value]. Default is "topleft"
#' @param plotReps [value]. Default is FALSE
#' @param plotVg [value]. Default is TRUE
#' @param plotSelInt [value]. Default is TRUE
#' @return [value]
#' @details [fill in details here]
#' @examples none
#' @export
simPlot <- function(popList, cols = "#000000", popLabs = NULL, varLine = "none", meanVariety = TRUE, legendPos = "topleft", plotReps = FALSE, plotVg = TRUE, plotSelInt = TRUE){

    if (length(cols) != length(popList)) stop("cols must be same length as popList!")
    lineCol <- paste0(cols, "FF")
    ptCol <- paste0(cols, "FF")
    polyCol <- paste0(cols, "4D")

    if (is.null(popLabs)) popLabs <- names(popList)

  	nVar <- unique(unlist(rlapply(popList, "[[", i = "nVar", level = 2, combine = c)))
    cyclePerYr <- unique(unlist(rlapply(popList, function(x) x[["paramL"]][["cyclePerYr"]], level = 2, combine = c)))
    if(length(nVar) > 1) warning("number of varieties differ between pops!")
    if(length(nVar) > 1) warning("cycles per year differ between pops!")

    simStats <- lapply(popList, function(x) lapply(x, function(xx) xx[!names(xx) %in% c("SP", "paramL", "VgVDP", "gvVDP")]))

    simStatsInv <- lapply(simStats, invertList)
    simReps <- rlapply(simStatsInv, level = 3, combine = rbind) 
    simAvg <- rlapply(simReps, f = colMeans, level = 2, na.rm = TRUE)

    RCRSyr <- simAvg[[1]][["RCRSyr"]]
    RCRSgen <- simAvg[[1]][["Rcyc"]]
    yr <- RCRSyr / cyclePerYr
    xlims <- range(c(0, RCRSgen))
    ylims <- range(sapply(simReps, getYrange)) * 1.1

    if (meanVariety) {
    	simAvg <- lapply(simAvg, function(x) {x[["vy"]] <- x[["varMean"]]; x[["vx"]] <- x[["RCRSyr"]]; x})
    }
    plot(NA, xlim = xlims, ylim = ylims, xaxt = "n", xlab = "generation", ylab = "standardized genetic value")
    axis(1, at = c(0, RCRSyr), labels = c(0, yr))

    for (i in 1:length(popList)){
    	if(plotReps) invisible(lapply(simStats[[i]], plotPop, Rgen = RCRSgen, popcol = cols[i]))
    	plotPop(simAvg[[i]], popcol = cols[i], alpha = "FF", alphaMean = "4D", Rgen = RCRSgen, vLine = varLine, pch = 16)
    }

    if (length(popList) > 1) {
        legend(legendPos, legend = popLabs, col=cols, lty = 1, lwd = 2, pch = 16)
      } else {
        legend(legendPos, legend = c("RCRS mean", expression(paste('RCRS ', sigma[g])), "Variety mean"), 
         bty = "n",
         col = c(lineCol, "black", ptCol),
         lty = c(1, 0, 0), lwd = c(2, 0, 0),
         pch = c(NA, 22, 16),
         pt.bg = c(NA, polyCol, ptCol),
        )
	}
	if(plotVg){
		ylims2 <- c(0, max(sapply(simAvg, "[[", "VgRCRS")))
		plot(NA, xlim = xlims, ylim = ylims2, xaxt = "n", xlab = "generation", ylab = "Vg", main = "Genetic variance across generations")
	    axis(1, at = c(0, RCRSyr), labels = c(0, yr))
	    for (i in 1:length(simAvg)) {
	    	lines(simAvg[[i]][["Rcyc"]], simAvg[[i]][["VgRCRS"]], type = "l", lwd = 2, lty = 1, col = cols[[i]])
	    }	
	    legend("topright", legend = popLabs, col=cols, lty = 1, lwd = 2, pch = 16)
	}



# THIS NEEDS TO BE CLEANED UP!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	if(plotSelInt){
	
		sVDP <- lapply(simAvg, "[[", "sVDP")
		iVDP <- lapply(simAvg, "[[", "iVDP")
		sRCRS <- lapply(simAvg, "[[", "sRCRS")
		iRCRS <- lapply(simAvg, "[[", "iRCRS")

		ylims3 <- range(unlist(sVDP))
		# ylims4 <- range(unlist(int))
		ylims4 <- c(-2, max(unlist(iVDP)))

		ylims5 <- range(unlist(sRCRS))
		# ylims4 <- range(unlist(int))
		ylims6 <- c(-2, max(unlist(iRCRS)))


		plot(NA, xlim = xlims, ylim = ylims3, xaxt = "n", xlab = "generation", ylab = "S", main = "Selection differential across generations in VDP")
	    axis(1, at = c(0, RCRSyr), labels = c(0, yr))
	    
	    for (i in 1:length(sRCRS)) {
	    	lines(RCRSyr, sVDP[[i]], type = "l", lwd = 2, lty = 1, col = cols[[i]])
	    }	
	    legend("topright", legend = popLabs, col=cols, lty = 1, lwd = 2, pch = 16)

	   	plot(NA, xlim = xlims, ylim = ylims4, xaxt = "n", xlab = "generation", ylab = "S/Vg", main = "Selection intensity across generations in VDP")
	    axis(1, at = c(0, RCRSyr), labels = c(0, yr))
	    
	    for (i in 1:length(iVDP)) {
	    	lines(RCRSyr, iVDP[[i]], type = "l", lwd = 2, lty = 1, col = cols[[i]])
	    }

	    legend("topleft", legend = popLabs, col=cols, lty = 1, lwd = 2, pch = 16)
		plot(NA, xlim = xlims, ylim = ylims5, xaxt = "n", xlab = "generation", ylab = "S", main = "Selection differential across generations in RCRS")
	    axis(1, at = c(0, RCRSyr), labels = c(0, yr))
	    
	    for (i in 1:length(sRCRS)) {
	    	lines(RCRSgen[-1], sRCRS[[i]], type = "l", lwd = 2, lty = 1, col = cols[[i]])
	    }	
	    legend("topright", legend = popLabs, col=cols, lty = 1, lwd = 2, pch = 16)

	   	plot(NA, xlim = xlims, ylim = ylims6, xaxt = "n", xlab = "generation", ylab = "S/Vg", main = "Selection intensity across generations in RCRS")
	    axis(1, at = c(0, RCRSyr), labels = c(0, yr))
	    
	    for (i in 1:length(iRCRS)) {
	    	lines(RCRSgen[-1], iRCRS[[i]], type = "l", lwd = 2, lty = 1, col = cols[[i]])
	    }
	    legend("topleft", legend = popLabs, col=cols, lty = 1, lwd = 2, pch = 16)
	}
}
