
getPopStats <- function(pop, meanVariety = TRUE){
    VgRGSC <- sapply(pop[["RGSC"]], genicVarG)
    gvRGSC <- sapply(pop[["RGSC"]], function(x) mean(gv(x)))
    # gvVariety <- sapply(run1[["VDP"]][["variety"]], function(x) mean(gv(x)))
    gvVariety <- lapply(pop[["VDP"]][["variety"]], function(x) gv(x))
    SDgRGSC <- sqrt(VgRGSC)
   
    # meanVariety <- TRUE
    if(meanVariety){
      Yvariety <- sapply(gvVariety, mean)
      Xvariety <- RGSCyr[1:length(Yvariety)]
    } else {
      nVariety <- sapply(gvVariety, nrow)
      Yvariety <- unlist(gvVariety)
      Xvariety <- rep(RGSCyr[1:length(nVariety)], times = nVariety)
    }

    return(list(Vg = VgRGSC, gv = gvRGSC, sd = SDgRGSC, vx = Xvariety, vy = Yvariety))
}

getYrange <- function(simR) range(c(simR$gv + simR$sd, simR$gv - simR$sd, simR$vy))

invertList <-  function(ll) {
    nms <- unique(unlist(lapply(ll, function(X) names(X))))
    ll <- lapply(ll, function(X) setNames(X[nms], nms))
    ll <- apply(do.call(rbind, ll), 2, as.list)
    lapply(ll, function(X) X[!sapply(X, is.null)])
}

plotPop <- function(simL, popcol = "#000000", alpha = "0D", alphaMean = "0D"){
    polycol <- paste0(popcol, alphaMean)
    popcol <- paste0(popcol, alpha)
    # attach(simL)
    xpoly <- c(RGSCgen, rev(RGSCgen), RGSCgen[1])
    ypoly <- c(simL$gv + simL$sd, rev(simL$gv - simL$sd), simL$gv[1] + simL$sd[1])

    polygon(x = xpoly, y = ypoly, col = polycol, border = NA)
    lines(x = RGSCgen, y = simL$gv, type = "l", col = popcol, lwd = 2)
    points(simL$vx, simL$vy, col = popcol, pch = 1)
  }

simPlot <- function(popList, baseCols){
    if(class(popList[[1]][[1]][[1]]) == "Pop") popList <- list(popList) # checks that popList is a list of pop

    simStats <- rlapply(popList, f=getPopStats, level = 2)
    simStatsInv <- lapply(simStats, invertList)

    simReps <- rlapply(simStatsInv, level = 3, combine = rbind) 
    simAvg <- rlapply(simReps, f = colMeans, level = 2)

    yr <- 1:nYr
    RGSCgen <- c(0, 1:(tail(yr * GScylcePerYr, 1)))

    xlims <- range(c(0, RGSCgen))
    ylims <- range(sapply(simReps, getYrange)) * 1.1

    plot(NA, xlim = xlims, ylim = ylims, xaxt = "n", xlab = "generation", ylab = "standardized genetic value")
    axis(1, at = c(0, RGSCyr), labels = c(0, RGSCyrLab))

    for(i in 1:length(popList)){
      invisible(lapply(simStats[[i]], plotPop))
      plotPop(simAvg[[i]], popcol = baseCols[i], alpha = "FF", alphaMean = "4D")
    }
  legend("topright", legend = c("RGSC mean", expression(paste('RGSC ', sigma[g])), "Variety mean"), 
         # fill = c(NA, polyCol, NA), 
         bty = "n",
         col = c(lineCol, "black", ptCol),
         lty = c(1, 0, 0), lwd = c(2, 0, 0),
         pch = c(NA, 22, 16),
         pt.bg = c(NA, polyCol, ptCol),
  )
}

popList <- list(simrun, simrun)

pdf("genValByGenTEST.pdf", width = 10, height = 7)
# 
simPlot(popList, c("#00FF00", "#0000FF"))
dev.off()

