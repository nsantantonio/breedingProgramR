#' getPopStats function
#'
#' function to (do something)
#'
#' @param resultL [value]
#' @param meanVariety [value]. Default is TRUE
#' @param verbose [value]. Default is FALSE
#' @return [value]
#' @details [fill in details here]
#' @examples none
#' @export
getPopStats <- function(resultL, meanVariety = TRUE, verbose = FALSE){
    VDPparam <- rlapply(resultL[["VDP"]], f = genParam, level = 2)
    RCRSparam <- lapply(resultL[["RCRS"]], genParam)

    nYr <- length(resultL[["VDP"]][[1]])
    GScylcePerYr <- (length(resultL[["RCRS"]]) - 1) / nYr 
    yr <- 1:nYr
    Ryr <- yr * GScylcePerYr
    Rcyc <- c(0, 1:(GScylcePerYr * nYr))

    VgRCRS <- sapply(RCRSparam, "[[", "varG")
    gvRCRS <- sapply(RCRSparam, function(x) mean(x$gv_a) + x$gv_mu) # this is correct
    sRCRS <- gvRCRS[-1] - gvRCRS[-length(gvRCRS)]
    iRCRS <- sRCRS / sqrt(VgRCRS[-length(VgRCRS)])

    VgVDP <- rlapply(VDPparam, "[[", i = "varG", level = 2, combine = c)
    gvVDP <- rlapply(VDPparam, function(x) mean(x$gv_a) + x$gv_mu, level = 2, combine = c)
    sVDP <- gvVDP$variety - gvVDP$trial1
    iVDP <- sVDP / sqrt(VgVDP$trial1)
  
  	RyrIndex <- Ryr - GScylcePerYr + 1
  	sTotal <- gvVDP$variety - gvRCRS[RyrIndex]
	iTotal <- sTotal / sqrt(VgRCRS[RyrIndex])

    gvVariety <- lapply(VDPparam[["variety"]], function(x) x$gv_a + x$gv_mu)
    SDgRCRS <- sqrt(VgRCRS)

    varMean <- gvVDP[["variety"]]
	nVariety <- sapply(gvVariety, nrow)
	Yvariety <- unlist(gvVariety)
	Xvariety <- rep(Ryr[1:length(nVariety)], times = nVariety)


    RCRSacc <- resultL$predAcc[["RCRS"]]
    VDPacc <- resultL$predAcc[["VDP"]]
    RCRSoutAcc <- resultL$predAcc[["RCRSout"]]
    VDPinAcc <- resultL$predAcc[["VDPin"]]

    theorMax <- maxBv(resultL$SP)
	
	nVar = unique(nVariety)

    return(list(SP = resultL$SP, paramL = resultL$paramL, Rcyc = Rcyc, varMean = varMean, sdRCRS = SDgRCRS, 
				VgRCRS = VgRCRS, VgVDP = VgVDP, gvRCRS = gvRCRS, gvVDP = gvVDP,
    			sRCRS = sRCRS, iRCRS = iRCRS, sVDP = sVDP, iVDP = iVDP, sTotal = sTotal, iTotal = iTotal, 
    			nVar = nVar, vx = Xvariety, vy = Yvariety, RCRSyr = Ryr, RCRSacc = RCRSacc, VDPacc = VDPacc,
    			RCRSoutAcc = RCRSoutAcc, VDPinAcc = VDPinAcc, theorMax = theorMax))
}
