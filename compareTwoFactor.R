parDir <- getwd()
source(paste0(parDir, "/alphaTools.R"))

defArgs <- list(f1 = NULL, f2 = NULL, simName = NULL, slice = "yr", projDir = "")
defArgs <- list(f1 = c("rgsc0.5", "rgsc0.4", "rgsc0.3", "rgsc0.2", "rgsc0.1", "rgsc0.05"), 
				f2 = c("vdp4x25", "vdp10x50", "vdp20x100", "vdp40x200"), 
				simName = "VDPvsRGSCintensity", 
				factColor = "f2", projDir = "VDPvsRGSCintensity/")

stdArgNames <- names(defArgs)
defArgs <- getComArgs(defArgs)
altArgs <- names(defArgs)[!names(defArgs) %in% stdArgNames] 
attach(defArgs)

# argList <- as.list(commandArgs(TRUE))
# argList <- c("truncSel_truncCross" "truncSel_expDistPairs" "truncSel_simDHdistPairs")

popList <- list()
for(i in f1){
	for(j in f2){
		simNameij <- paste(simName, i, j, sep = "_")
		load(paste0(parDir, "/results/", projDir, simNameij, "/", simNameij, ".RData"))
		if(factColor == "f1") popList[[i]][[j]] <- simrun else  popList[[j]][[i]] <- simrun 
	}
}

cols <- c("#006600", "#00008C", "#660000", "#8C8C00", "#660066", "#006666")
cols <- cols[1:length(popList)]
# cols <- if(length(argList) == 1) "#000000" else if(length(argList) <= 6) cols[1:length(argList)] else stop("I can only do up to 4 sims!")




twoFactorPlot <- function(popList, cols = "#000000", fLabs = NULL, varLine = "none", meanVariety = FALSE, legendPos = "topleft", plotReps = FALSE, plotVg = TRUE, plotSelInt = TRUE){

    if (length(cols) != length(popList)) stop("cols must be same length as popList!")
    lineCol <- paste0(cols, "FF")
    ptCol <- paste0(cols, "FF")
    polyCol <- paste0(cols, "4D")

    if (is.null(fLabs)) fLabs <- names(popList)

    nVar <- unique(round(unlist(rlapply(popList, level = 3, f = function(x) {tail(with(x[["paramL"]], nFam * famSize * cumprod(selectTrials)), 1)}, combine = c))))
    cyclePerYr <- unique(round(unlist(rlapply(popList, level = 3, f = function(x) {x[["paramL"]][["cyclePerYr"]]}, combine = c))))

    simStats <- rlapply(popList, function(x) {x[!names(x) %in% c("SP", "paramL", "VgVDP", "gvVDP")]}, level = 3)
    simStatsInv <- rlapply(simStats, invertList, level = 2)
    simReps <- rlapply(simStatsInv, level = 4, combine = rbind) 
    simAvg <- rlapply(simReps, f = colMeans, level = 3, na.rm = TRUE)

    RGSCyr <- simAvg[[1]][[1]][["RGSCyr"]]
    RGSCgen <- simAvg[[1]][[1]][["Rcyc"]]
    yr <- RGSCyr / cyclePerYr
    xlims <- range(c(0, RGSCgen))
    ylims <- range(do.call(rbind, rlapply(simReps, getYrange, level = 2, combine = rbind))) * 1.1

    if (meanVariety) {
    	simAvg <- rlapply(simAvg, function(x) {x[["vy"]] <- x[["varMean"]]; x[["vx"]] <- x[["RGSCyr"]]; x}, level = 2)
    }

	varL <- rlapply(simAvg, "[[", level = 2, i = "vy", rbind)
	gsL <- rlapply(simAvg, "[[", level = 2, i = "gvRGSC", rbind)
	SL <- rlapply(simAvg, "[[", level = 2, i = "sVDP", rbind)
	iL <- rlapply(simAvg, "[[", level = 2, i = "iVDP", rbind)

	# lapply(iL, rowMeans)
	# getIntensity

	varCoord <- list(x = 1:length(f1names), z = yr, y= varL[[i]])
	RGSCcoord <- list(x = 1:length(f1names), z = RGSCgen, y= gsL[[i]])

	ylimits <- range(0, unlist(c(varL, gsL)))
	min(sqrt())
	
	pdf("forLabMeeting1.pdf")
	par(mfrow = c(4, 5))
	par(mar = c(0, 0, 0, 0), oma = c(4, 5, 0.5, 0.5))
	par(tcl = -0.25)
	par(mgp = c(2, 0.6, 0))
	for(k in yr){
		# if(k %in% 16:20) prx <- 's' else prx <- 'n'
		if(k %in% c(1, 6, 11, 16)) pry <- 's' else pry <- 'n'
		plot(1:length(f1names), varL[[1]][, k], type = "l", ylim = ylimits, xaxt = "n", col = cols[1], yaxt = pry)
		text(2, 14, paste0("year ", k))
		if(k %in% 16:20) axis(1, at = 1:length(f1names), labels = gsub("rgsc", "", f1names))
		for(i in 2:length(varL)){
		lines(1:length(f1names), varL[[i]][, k], col = cols[i])
		}
	}
	mtext("Recurrent GS intensity", 1,  cex = 1.7, padj = 1.6, outer = TRUE)
	# mtext(expression(r^2), 2,  cex = 2, padj = -1.5, outer = TRUE, las = 2)
	mtext("Genetic Value of Varieties", 2,  cex = 1.5, adj = 2.1, padj = 1.5, outer = TRUE)

	# dev.off()


	# pdf("forLabMeeting2.pdf")
	par(mfrow = c(4, 5))
	par(mar = c(0, 0, 0, 0), oma = c(4, 5, 0.5, 0.5))
	par(tcl = -0.25)
	par(mgp = c(2, 0.6, 0))
	for(k in RGSCyr){
		# if(k %in% 16:20) prx <- 's' else prx <- 'n'
		if(k %in% c(2, 12, 22, 32)) pry <- 's' else pry <- 'n'
		plot(1:length(f1names), gsL[[1]][, k], type = "l", ylim = ylimits, xaxt = "n", col = cols[1], yaxt = pry, lty = 2)
		text(2, 14, paste0("year ", k/2))
		if(k %in% c(32, 34, 36, 38, 40)) axis(1, at = 1:length(f1names), labels = gsub("rgsc", "", f1names))
		for(i in 2:length(gsL)){
		lines(1:length(f1names), gsL[[i]][, k], col = cols[i], lty = 2)
		}
	}
	mtext("Recurrent GS intensity", 1,  cex = 1.7, padj = 1.6, outer = TRUE)
	# mtext(expression(r^2), 2,  cex = 2, padj = -1.5, outer = TRUE, las = 2)
	mtext("Genetic Value of Varieties", 2,  cex = 1.5, adj = 2.1, padj = 1.5, outer = TRUE)

	dev.off()




	f1names <- names(popList[[1]]) 
	alpha = "0D"
	alphaMean = "4D"
	angle <- 180
	yangle <- 20

	# for(i in 1:length(varL)){
	# 	if(i == 1) {
	# 		persp(x = 1:length(f1names), y = yr, z = varL[[i]], col = paste0(cols[[i]], alpha), border = paste0(cols[[i]], alphaMean), theta = angle, phi = yangle, 
	# 			xlab = simName, ylab = "Year", zlab = "genetic value")
	# 	} else {
	# 		persp(varL[[i]], col = paste0(cols[[i]], alpha), border = paste0(cols[[i]], alphaMean), theta = angle, phi = yangle, box = FALSE)
	# 	}
	# 	if(i != length(varL)) par(new=TRUE)
	# }

	maxy <- max(unlist(varL))
	library(rgl)
	rgl.open()
	rgl.bg(color = "white") 

	rgl.lines(c(0, length(f1names)), c(0, 0), c(0, 0), color = "black")
	rgl.lines(c(0, 0), c(0, maxy), c(0, 0), color = "black")
	rgl.lines(c(0, 0), c(0, 0), c(0,max(yr)), color = "black")
	# rgl.texts(1:length(f1names), text = f1names, color = "black", adj = c(0.5, -0.8), size = 2)
	for(i in 1:length(varL)){
		rgl.surface(x = 1:length(f1names), z = yr, y= varL[[i]], col = rep(paste0(cols[[i]]), prod(dim(varL[[i]]))), alpha = 0.25)#, border = paste0(cols[[i]], alphaMean),)
	}

	library(rgl)
	rgl.open()
	rgl.bg(color = "white") 

	rgl.lines(c(0, length(f1names)), c(0, 0), c(0, 0), color = "black")
	rgl.lines(c(0, 0), c(0, maxy), c(0, 0), color = "black")
	rgl.lines(c(0, 0), c(0, 0), c(0,max(yr)), color = "black")
# rgl.texts(1:length(f1names), text = f1names, color = "black", adj = c(0.5, -0.8), size = 2)
	for(i in 1:length(gsL)){
		rgl.surface(x = 1:length(f1names), z = RGSCgen, y= gsL[[i]], col = rep(paste0(cols[[i]]), prod(dim(gsL[[i]]))), alpha = 0.25)#, border = paste0(cols[[i]], alphaMean),)
	}



Slice on Year

rgl.clear()


    # plot means of RGSC and VDP output
    plot(NA, xlim = xlims, ylim = ylims, xaxt = "n", xlab = "generation", ylab = "standardized genetic value")
    axis(1, at = c(0, RGSCyr), labels = c(0, yr))

    for (i in 1:length(popList)){
    	if(plotReps) invisible(lapply(simStats[[i]], plotPop, Rgen = RGSCgen, popcol = cols[i]))
    	plotPop(simAvg[[i]], popcol = cols[i], alpha = "FF", alphaMean = "4D", Rgen = RGSCgen, vLine = varLine, pch = 16)
    }

    if (length(popList) > 1) {
        legend(legendPos, legend = fLabs, col=cols, lty = 1, lwd = 2, pch = 16)
      } else {
        legend(legendPos, legend = c("RGSC mean", expression(paste('RGSC ', sigma[g])), "Variety mean"), 
         bty = "n",
         col = c(lineCol, "black", ptCol),
         lty = c(1, 0, 0), lwd = c(2, 0, 0),
         pch = c(NA, 22, 16),
         pt.bg = c(NA, polyCol, ptCol),
        )
	}
	if(plotVg){
		ylims2 <- c(0, max(sapply(simAvg, "[[", "Vg")))
		plot(NA, xlim = xlims, ylim = ylims2, xaxt = "n", xlab = "generation", ylab = "Vg", main = "Genetic variance across generations")
	    axis(1, at = c(0, RGSCyr), labels = c(0, yr))
	    for (i in 1:length(simAvg)) {
	    	lines(simAvg[[i]][["Rcyc"]], simAvg[[i]][["Vg"]], type = "l", lwd = 2, lty = 1, col = cols[[i]])
	    }	
	    legend("topright", legend = fLabs, col=cols, lty = 1, lwd = 2, pch = 16)
	}


# check for selection intensity off by one
	if(plotSelInt){
		selInt <- lapply(simAvg, getIntensity)
		S <- lapply(selInt, "[[", "S")
		int <- lapply(selInt, "[[", "i")
		ylims3 <- range(unlist(S))
		# ylims4 <- range(unlist(int))
		ylims4 <- c(-2, max(unlist(int)))

		plot(NA, xlim = xlims, ylim = ylims3, xaxt = "n", xlab = "generation", ylab = "S", main = "Selection differential across generations")
	    axis(1, at = c(0, RGSCyr), labels = c(0, yr))
	    
	    for (i in 1:length(S)) {
	    	lines(RGSCyr, S[[i]], type = "l", lwd = 2, lty = 1, col = cols[[i]])
	    }	
	    legend("topright", legend = fLabs, col=cols, lty = 1, lwd = 2, pch = 16)

	   	plot(NA, xlim = xlims, ylim = ylims4, xaxt = "n", xlab = "generation", ylab = "S/Vg", main = "Selection intensity across generations")
	    axis(1, at = c(0, RGSCyr), labels = c(0, yr))
	    
	    for (i in 1:length(int)) {
	    	lines(RGSCyr, int[[i]], type = "l", lwd = 2, lty = 1, col = cols[[i]])
	    }
	    legend("topleft", legend = fLabs, col=cols, lty = 1, lwd = 2, pch = 16)
	}
}

















params <- lapply(popList, function(x) x[[1]][["paramL"]])
llunion <- function(ll){
	lunion <- function(l1, l2) l1[union(names(l1)[l1 %in% l2], names(l2)[l2 %in% l1])]
	inCom <- Reduce(lunion, ll)
	notInCom <- lapply(ll, function(x) x[!x %in% inCom])
	unique(unlist(lapply(notInCom, names)))
}

# lapply(params, "[[", "useTrue")

ldiff <- llunion(params)
print(ldiff)

pdfName <- paste0(parDir, "/figures/", paste(argList, collapse = "_"), ".pdf") 
pdf(pdfName, width = 12, height = 7)
simPlot(popList, cols, varLine = "poly")
dev.off()

if(system("hostname", intern = TRUE) == "cbsurobbins.biohpc.cornell.edu") system(paste0("scp ", pdfName, " Bender:~/Dropbox/optibreedSim/figures/"))