
formatPop <- function(popL, depth = 1, meanVariety = TRUE){
		get1 <- function(l, what, depth, counter = 0) if(depth == counter) return(l[[what]]) else get1(l[[1]], what, depth, counter = counter + 1)
		get1names <- function(depth, l, counter = 0) if(depth == counter) return(names(l)) else get1names(depth, l[[1]], counter = counter + 1)

	    nVar <- unique(round(unlist(rlapply(popL, level = depth, f = function(x) {tail(with(x[["paramL"]], nFam * famSize * cumprod(selectTrials)), 1)}, combine = c))))
	    cyclePerYr <- unique(round(unlist(rlapply(popL, level = depth, f = function(x) {x[["paramL"]][["cyclePerYr"]]}, combine = c))))

	    simStats <- rlapply(popL, function(x) {x[!names(x) %in% c("SP", "paramL", "VgVDP", "gvVDP")]}, level = depth)
	    simStatsInv <- rlapply(simStats, invertList, level = depth - 1)
	    simReps <- rlapply(simStatsInv, level = depth + 1, combine = rbind) 
	    simReps <- rlapply(simReps, level = depth - 1, f = function(x) x[!sapply(x, is.null)]) 
	    simAvg <- rlapply(simReps, f = colMeans, level = depth, na.rm = TRUE)

	    RGSCyr <- get1(simAvg, "RGSCyr", depth - 1)
	    RGSCgen <- get1(simAvg, "Rcyc", depth - 1)
	    yr <- RGSCyr / cyclePerYr
	    xlims <- range(c(0, RGSCgen))
	    ylims <- range(unlist(rlapply(simReps, getYrange, level = depth - 1))) * 1.1

	    if (meanVariety) {
	    	simAvg <- rlapply(simAvg, function(x) {x[["vy"]] <- x[["varMean"]]; x[["vx"]] <- x[["RGSCyr"]]; x}, level = depth - 1)
	    }

		varL <- rlapply(simAvg, "[[", level = depth - 1, i = "vy", rbind)
		gsL <- rlapply(simAvg, "[[", level = depth - 1, i = "gvRGSC", rbind)
		Vg <- rlapply(simAvg, "[[", level = depth - 1, i = "VgRGSC", rbind)
		SL <- rlapply(simAvg, "[[", level = depth - 1, i = "sVDP", rbind)
		iL <- rlapply(simAvg, "[[", level = depth - 1, i = "iVDP", rbind)

		# lapply(iL, rowMeans)
		# getIntensity
		factorNames <- lapply(0:(depth - 2), get1names, l = simAvg)
		names(factorNames) <- names(factors)

		list(nVar = nVar, cyclePerYr = cyclePerYr, RGSCyr = RGSCyr, RGSCgen = RGSCgen, yr = yr, xlims = xlims, ylims = ylims, 
			 varL = varL, Vg = Vg, gsL = gsL, SL = SL, iL = iL, factorNames = factorNames)
	}

# popL <- popList[[1]]; depth = 6; fLabs = NULL; plotReps = FALSE; plotVg = TRUE; plotSelInt = TRUE


yearPanes <- function(statList, val, flist, byF, x, xnames, colCat, title, xlabel, ylabel, cols = c("#006600", "#00008C", "#660000", "#8C8C00", "#660066", "#006666"), meanVariety = TRUE){
	if(length(val) == 1) val <- rep(val, 2)
    # if (length(cols) != length(popL)) stop("cols must be same length as popL!")
    # lineCol <- paste0(cols, "FF")
    # ptCol <- paste0(cols, "FF")
    # polyCol <- paste0(cols, "4D")

    # if (is.null(fLabs)) fLabs <- names(popL)

	
	
	# statList <- formatPop(popList[[1]], depth = 5)
	# statList <- formatPop(popList, depth = 6)

	# names(statList$factorNames) <- names(factors)
	# statList$Vg

	# byF <- factors$select
	# colCat <- factors$VDPint
	cols <- c("#006600", "#00008C", "#660000", "#8C8C00", "#660066", "#006666")
	cols <- cols[1:length(colCat)]

	nF <- sapply(flist, length) 
	if(sum(nF > 1) > 1) stop("I can only handle two factor levels for one factor!")
	justOne <- if(all(nF == 1)) TRUE else FALSE

	if(!justOne) flist[nF == 1] <- lapply(flist[nF == 1], function(x) rep(x, 2)) 

	# val2 <- "gsL"

	# title <- "title"
	f1 <- flist[[1]]
	f2 <- flist[[2]]
	f3 <- flist[[3]]
	y1 <- statList[[val[1]]][[f1[1]]][[f2[1]]][[f3[1]]]
	yl <- list(y1)
	if(!justOne) {
		y2 <- statList[[val[2]]][[f1[2]]][[f2[2]]][[f3[1]]]
		yl <- c(yl, list(y2))
	}
	# z1 <- statList[[val2]][["1000QTL"]][["trunc"]][["rgsc"]]
	# z2 <- statList[[val2]][["1000QTL"]][["ACquant"]][["rgsc"]]
	# yl <- list(y1, y2, z1, z2)


	ylimits <- range(0, unlist(yl))
	# xnames <- statList$factorNames$crossInt
	# pdf("truncSel.pdf")
	par(mfrow = c(5, 6))
	par(mar = c(0, 0, 0, 0), oma = c(4, 5, 5, 0.5))
	par(tcl = -0.25)
	par(mgp = c(2, 0.6, 0))
	for(k in statList$yr){
		# if(k %in% 16:20) prx <- 's' else prx <- 'n'
		if(k %in% c(1, 7, 13, 19, 25)) pry <- 's' else pry <- 'n'
		# x <- 1 - as.numeric(xnames)
		# plot(x, y1[[1]][, k], type = "l", ylim = ylimits, xaxt = "n", col = cols[1], yaxt = pry)
		plot(NA, type = "l", xlim = range(x), ylim = ylimits, xaxt = "n", col = cols[1], yaxt = pry)
		text(2, 14, paste0("year ", k))
		if(k %in% 25:30) axis(1, at = x, labels = gsub("rgsc", "", xnames))
	
		for(j in 1:length(yl)){
			y <- yl[[j]]
			for(i in 1:length(y)){
				yik <- y[[i]][, k]
				best <- which.max(yik)
				lines(x, yik, col = cols[i], lty = j)
				points(x[best], yik[best], col = cols[i], pch = c(16, 1)[j])
			}
		}
	}
	mtext(xlabel, 1,  cex = 1, padj = 2.5, outer = TRUE)
	mtext(ylabel, 2,  cex = 1, adj = 0.5, padj = -2, outer = TRUE)
	mtext(title, 3,  cex = 1, padj = -1, outer = TRUE)
}




intensityPanes <- function(statList, val, flist, panes, x, xnames, colCat, title, xlabel, ylabel, cols = c("#006600", "#00008C", "#660000", "#8C8C00", "#660066", "#006666"), meanVariety = TRUE){
	if(length(val) == 1) val <- rep(val, 2)

	cols <- c("#006600", "#00008C", "#660000", "#8C8C00", "#660066", "#006666")
	cols <- cols[1:length(colCat)]

	nF <- sapply(flist, length) 
	if(sum(nF > 1) > 1) stop("I can only handle two factor levels for one factor!")
	justOne <- if(all(nF == 1)) TRUE else FALSE

	if(!justOne) flist[nF == 1] <- lapply(flist[nF == 1], function(x) rep(x, 2)) 

	f1 <- flist[[1]]
	f2 <- flist[[2]]
	f3 <- flist[[3]]
	y1 <- statList[[val[1]]][[f1[1]]][[f2[1]]][[f3[1]]]
	yl <- list(y1)
	if(!justOne) {
		y2 <- statList[[val[2]]][[f1[2]]][[f2[2]]][[f3[1]]]
		yl <- c(yl, list(y2))
	}

	bestL <- list()
	for (i in 1:length(yl)) {
		bestL[[i]] <- list()
		for (j in 1:length(yl[[i]])) {
			bestL[[i]][[j]]<- apply(yl[[i]][[j]], 2, which.max)
		}
	}

	ylimits <- range(0, unlist(yl))

	par(mfrow = c(length(panes), 1))
	par(mar = c(0, 0, 0, 0), oma = c(4, 5, 5, 0.5))
	par(tcl = -0.25)
	par(mgp = c(2, 0.6, 0))
	for(k in 1:length(panes)){
		# if(k %in% 16:20) prx <- 's' else prx <- 'n'
		if(k %in% c(1:length(x))) pry <- 's' else pry <- 'n'
		# x <- 1 - as.numeric(xnames)
		# plot(x, y1[[1]][, k], type = "l", ylim = ylimits, xaxt = "n", col = cols[1], yaxt = pry)
		plot(NA, type = "l", xlim = range(x), ylim = ylimits, xaxt = "n", col = cols[1], yaxt = pry)
		text(3, 0.9 * ylimits[2], paste0("RGSC intensity: ", panes[k]))
		if(k %in% length(panes)) axis(1, at = x, labels = xnames)
	
		for(j in 1:length(yl)){
			y <- yl[[j]]
			for(i in 1:length(y)){
				yik <- y[[i]][k, ]
				best <- bestL[[j]][[i]] == k
				lines(x, yik, col = cols[i], lty = j)
				points(x[best], yik[best], col = cols[i], pch = c(16, 1)[j])
			}
		}
	}
	mtext(xlabel, 1,  cex = 1, padj = 2.5, outer = TRUE)
	mtext(ylabel, 2,  cex = 1, adj = 0.5, padj = -2, outer = TRUE)
	mtext(title, 3,  cex = 1, padj = -1, outer = TRUE)
}





# }

# VDPvsRGSCintensity10QTLACquant_rgsc0.1_vdp4x25

# VDPvsRGSCintensityACquant_rgsc0.05_vdp10x50.RData
parDir <- getwd()
source(paste0(parDir, "/alphaTools.R"))

# defArgs <- list(f1 = NULL, f2 = NULL, simName = NULL, slice = "yr", projDir = "")
defArgs <- list(factors = list(QTL = c("1000QTL", "100QTL", "10QTL"),
							   select = c("", "ACquant"),
							   crossSel = c("rgsc", "trad"), 
							   crossInt = c("0.5", "0.4", "0.3", "0.2", "0.1", "0.05"), 
							   VDPint = c("vdp4x25", "vdp10x50", "vdp20x100", "vdp40x200")
							   ), 
				simName = "VDPvsRGSCintensity30yr", 
				projDir = "VDPvsRGSCintensityQTL30yr/")

byVDP <- TRUE
defArgs <- getComArgs(defArgs)
attach(defArgs)

# argList <- as.list(commandArgs(TRUE))
# argList <- c("truncSel_truncCross" "truncSel_expDistPairs" "truncSel_simDHdistPairs")

popList <- list()
for(i in factors$QTL){
	for(j in factors$select){
		for(k in factors$crossSel){
			for(l in factors$crossInt){
				for(m in factors$VDPint){
					simNameijklm <- paste0(c(simName, c(i, j, k, l, m)), c("", "", "_", "", "_", ""), collapse = "")
					load(paste0(parDir, "/results/", projDir, simNameijklm, "/", simNameijklm, ".RData"))
					sel <- if(j == "") "trunc" else j
					if(byVDP) popList[[i]][[sel]][[k]][[m]][[l]] <- simrun else popList[[i]][[sel]][[k]][[l]][[m]] <- simrun
					nullsim <- sapply(simrun, is.null)
					if(any(nullsim)) cat(simNameijklm, "has", sum(nullsim), "missing replicates!!!\n")
				}
			}
		}
	}
}
if(byVDP) factors <- factors[c(1:3, 5:4)] # reverse order for plotting

cols <- c("#006600", "#00008C", "#660000", "#8C8C00", "#660066", "#006666")
cols <- cols[1:length(popList[[1]])]



# cols <- if(length(argList) == 1) "#000000" else if(length(argList) <= 6) cols[1:length(argList)] else stop("I can only do up to 4 sims!")

# listDepth <- function(l, counter = 0){
# 	if(all(sapply(l, class) == "list")) {
# 	# if(is.list(l[[1]])) {
# 		counter <- listDepth(l[[1]], counter = counter + 1)
# 	}
# 	return(counter)
# }
# listDepth(popList)

# names(popList[[1]][[1]][[1]][[1]][[1]][[1]])


statList <- formatPop(popList, depth = 6)
val <- "varL"

byF <- factors$select
colCat <- factors$VDPint
flist <- list(f1 = "1000QTL", f2 = c("trunc", "ACquant"), f3 = "rgsc") 

xnames = statList$yr
x <- statList$yr

# x <- 1 - as.numeric(xnames)




# yearPanes(statList, val = "varL", flist = list(f1 = "1000QTL", f2 = c("trunc"), f3 = "rgsc"), 
# 		  byF = factors$select, x = 1- as.numeric(statList$factorNames$crossInt), xnames = statList$factorNames$crossInt,
# 		  colCat = factors$VDPint, title = "Truncation", xlabel = "RGSC intensity", ylabel = "Variety mean")



# intensityPanes(statList, val = "varL", flist = list(f1 = "1000QTL", f2 = c("trunc"), f3 = "rgsc"), 
# 		  panes = factors$crossInt, x = statList$yr, xnames = statList$yr,
# 		  colCat = factors$VDPint, title = "Truncation", xlabel = "Year", ylabel = "Variety mean")


# yearPanes(statList, val = "varL", flist = list(f1 = "1000QTL", f2 = c("ACquant"), f3 = "rgsc"), 
# 		  byF = factors$select, x = 1- as.numeric(statList$factorNames$crossInt), xnames = statList$factorNames$crossInt,
# 		  colCat = factors$VDPint, title = "Truncation vs. Ant Colony Quantile", xlabel = "RGSC intensity", ylabel = "Variety mean")


# intensityPanes(statList, val = "varL", flist = list(f1 = "1000QTL", f2 = c("ACquant"), f3 = "rgsc"), 
# 		  panes = factors$crossInt, x = statList$yr, xnames = statList$yr,
# 		  colCat = factors$VDPint, title = "Ant Colony Quantile", xlabel = "Year", ylabel = "Variety mean")


for(i in c("10QTL", "100QTL", "1000QTL")){
	pdf(paste0("figures/VDPvsRGSCintensityTrucVsACquant_", i, "byYear.RData"))
	yearPanes(statList, val = "varL", flist = list(f1 = i, f2 = c("trunc", "ACquant"), f3 = "rgsc"), 
		  byF = factors$select, x = 1- as.numeric(statList$factorNames$crossInt), xnames = statList$factorNames$crossInt,
		  colCat = factors$VDPint, title = "Truncation vs. Ant Colony Quantile", xlabel = "RGSC intensity", ylabel = "Variety mean")

	yearPanes(statList, val = "gsL", flist = list(f1 = i, f2 = c("trunc", "ACquant"), f3 = "rgsc"), 
		  byF = factors$select, x = 1- as.numeric(statList$factorNames$crossInt), xnames = statList$factorNames$crossInt,
		  colCat = factors$VDPint, title = "Truncation vs. Ant Colony Quantile", xlabel = "RGSC intensity", ylabel = "RGSC population mean")

	yearPanes(statList, val = "Vg", flist = list(f1 = i, f2 = c("trunc", "ACquant"), f3 = "rgsc"), 
			  byF = factors$select, x = 1- as.numeric(statList$factorNames$crossInt), xnames = statList$factorNames$crossInt,
			  colCat = factors$VDPint, title = "Truncation vs. Ant Colony Quantile", xlabel = "RGSC intensity", ylabel = "Genetic Variance in RGSC")
	dev.off()
	
	pdf(paste0("figures/VDPvsRGSCintensityTrucVsACquant_", i, "byInt.RData"), width = 5)

	intensityPanes(statList, val = "varL", flist = list(f1 = i, f2 = c("trunc", "ACquant"), f3 = "rgsc"), 
		  panes = factors$crossInt, x = statList$yr, xnames = statList$yr,
		  colCat = factors$VDPint, title = "Truncation vs. Ant Colony Quantile", xlabel = "Year", ylabel = "Variety mean")

	intensityPanes(statList, val = "gsL", flist = list(f1 = i, f2 = c("trunc", "ACquant"), f3 = "rgsc"), 
		  panes = factors$crossInt, x = statList$RGSCgen, xnames = statList$RGSCgen,
		  colCat = factors$VDPint, title = "Truncation vs. Ant Colony Quantile", xlabel = "Year", ylabel = "RGSC population mean")

	intensityPanes(statList, val = "Vg", flist = list(f1 = i, f2 = c("trunc", "ACquant"), f3 = "rgsc"), 
			  panes = factors$crossInt, x = statList$RGSCgen, xnames = statList$RGSCgen,
			  colCat = factors$VDPint, title = "Truncation vs. Ant Colony Quantile", xlabel = "RGSC generation", ylabel = "Genetic Variance in RGSC")
	dev.off()

}


# yearPanes(statList, val = "Vg", flist = list(f1 = "1000QTL", f2 = c("trunc", "ACquant"), f3 = "rgsc"), 
# 		  byF = factors$select, x = 1- as.numeric(statList$factorNames$crossInt), xnames = statList$factorNames$crossInt,
# 		  colCat = factors$VDPint, title = "Truncation vs. Ant Colony Quantile", xlabel = "RGSC intensity", ylabel = "Genetic Variance in RGSC")

# intensityPanes(statList, val = "Vg", flist = list(f1 = "1000QTL", f2 = c("trunc", "ACquant"), f3 = "rgsc"), 
# 		  panes = factors$crossInt, x = statList$RGSCgen, xnames = statList$RGSCgen,
# 		  colCat = factors$VDPint, title = "Truncation vs. Ant Colony Quantile", xlabel = "RGSC generation", ylabel = "Genetic Variance in RGSC")




# yearPanes(statList, val = "varL", flist = list(f1 = "100QTL", f2 = c("trunc", "ACquant"), f3 = "rgsc"), 
# 		  byF = factors$select, x = 1- as.numeric(statList$factorNames$crossInt), xnames = statList$factorNames$crossInt,
# 		  colCat = factors$VDPint, title = "Truncation vs. Ant Colony Quantile - 100 QTL", xlabel = "RGSC intensity", ylabel = "Variety mean")

# intensityPanes(statList, val = "varL", flist = list(f1 = "100QTL", f2 = c("trunc", "ACquant"), f3 = "rgsc"), 
# 		  panes = factors$crossInt, x = statList$yr, xnames = statList$yr,
# 		  colCat = factors$VDPint, title = "Truncation vs. Ant Colony Quantile - 100 QTL", xlabel = "Year", ylabel = "Variety mean")


# yearPanes(statList, val = "Vg", flist = list(f1 = "100QTL", f2 = c("trunc", "ACquant"), f3 = "rgsc"), 
# 		  byF = factors$select, x = 1- as.numeric(statList$factorNames$crossInt), xnames = statList$factorNames$crossInt,
# 		  colCat = factors$VDPint, title = "Truncation vs. Ant Colony Quantile - 100 QTL", xlabel = "RGSC intensity", ylabel = "Genetic Variance in RGSC")

# intensityPanes(statList, val = "Vg", flist = list(f1 = "100QTL", f2 = c("trunc", "ACquant"), f3 = "rgsc"), 
# 		  panes = factors$crossInt, x = statList$RGSCgen, xnames = statList$RGSCgen,
# 		  colCat = factors$VDPint, title = "Truncation vs. Ant Colony Quantile - 100 QTL", xlabel = "RGSC generation", ylabel = "Genetic Variance in RGSC")


# yearPanes(statList, val = "Vg", flist = list(f1 = "10QTL", f2 = c("trunc", "ACquant"), f3 = "rgsc"), 
# 		  byF = factors$select, x = 1- as.numeric(statList$factorNames$crossInt), xnames = statList$factorNames$crossInt,
# 		  colCat = factors$VDPint, title = "Truncation vs. Ant Colony Quantile", xlabel = "RGSC intensity", ylabel = "Genetic Variance in RGSC")

# intensityPanes(statList, val = "Vg", flist = list(f1 = "10QTL", f2 = c("trunc", "ACquant"), f3 = "rgsc"), 
# 		  panes = factors$crossInt, x = statList$RGSCgen, xnames = statList$RGSCgen,
# 		  colCat = factors$VDPint, title = "Truncation vs. Ant Colony Quantile", xlabel = "RGSC generation", ylabel = "Genetic Variance in RGSC")


# yearPanes(statList, val = "Vg", flist = list(f1 = "10QTL", f2 = c("trunc", "ACquant"), f3 = "rgsc"), 
# 		  byF = factors$select, x = 1- as.numeric(statList$factorNames$crossInt), xnames = statList$factorNames$crossInt,
# 		  colCat = factors$VDPint, title = "Truncation vs. Ant Colony Quantile", xlabel = "RGSC intensity", ylabel = "Genetic Variance in RGSC")

# intensityPanes(statList, val = "Vg", flist = list(f1 = "10QTL", f2 = c("trunc", "ACquant"), f3 = "rgsc"), 
# 		  panes = factors$crossInt, x = statList$RGSCgen, xnames = statList$RGSCgen,
# 		  colCat = factors$VDPint, title = "Truncation vs. Ant Colony Quantile", xlabel = "RGSC generation", ylabel = "Genetic Variance in RGSC")










statList <- formatPop(popList, depth = 6)
val <- "varL"

byF <- factors$select
colCat <- factors$VDPint
panes <- factors$crossInt
flist <- list(f1 = "1000QTL", f2 = c("trunc", "ACquant"), f3 = "rgsc") 

xnames = statList$yr
x <- statList$yr














pdf("truncSel.pdf")
twoFactorPlot(popListTrunc, cols = cols)
dev.off()



# VDPvsRGSCintensityACquant_rgsc0.05_vdp10x50.RData
parDir <- getwd()
source(paste0(parDir, "/alphaTools.R"))

defArgs <- list(f1 = NULL, f2 = NULL, simName = NULL, slice = "yr", projDir = "")
defArgs <- list(f1 = c("rgsc0.5", "rgsc0.4", "rgsc0.3", "rgsc0.2", "rgsc0.1", "rgsc0.05"), 
				f2 = c("vdp4x25", "vdp10x50", "vdp20x100", "vdp40x200"), 
				simName = "VDPvsRGSCintensityACquant", 
				factColor = "f2", projDir = "VDPvsRGSCintensityACquant/")

stdArgNames <- names(defArgs)
defArgs <- getComArgs(defArgs)
altArgs <- names(defArgs)[!names(defArgs) %in% stdArgNames] 
attach(defArgs)

# argList <- as.list(commandArgs(TRUE))
# argList <- c("truncSel_truncCross" "truncSel_expDistPairs" "truncSel_simDHdistPairs")

popListACquant <- list()
for(i in f1){
	for(j in f2){
		simNameij <- paste(simName, i, j, sep = "_")
		load(paste0(parDir, "/results/", projDir, simNameij, "/", simNameij, ".RData"))
		if(factColor == "f1") popListACquant[[i]][[j]] <- simrun else  popListACquant[[j]][[i]] <- simrun 
	}
}

cols <- c("#006600", "#00008C", "#660000", "#8C8C00", "#660066", "#006666")
cols <- cols[1:length(popListACquant)]
# cols <- if(length(argList) == 1) "#000000" else if(length(argList) <= 6) cols[1:length(argList)] else stop("I can only do up to 4 sims!")




pdf("ACquantSel.pdf")
twoFactorPlot(popListACquant, cols = cols)
dev.off()



threeFactorPlot <- function(popList, cols = "#000000", fLabs = NULL, varLine = "none", meanVariety = TRUE, legendPos = "topleft", plotReps = FALSE, plotVg = TRUE, plotSelInt = TRUE){

	formatPop <- function(popL){
	    if (length(cols) != length(popL)) stop("cols must be same length as popL!")
	    lineCol <- paste0(cols, "FF")
	    ptCol <- paste0(cols, "FF")
	    polyCol <- paste0(cols, "4D")

	    if (is.null(fLabs)) fLabs <- names(popL)

	    nVar <- unique(round(unlist(rlapply(popL, level = 3, f = function(x) {tail(with(x[["paramL"]], nFam * famSize * cumprod(selectTrials)), 1)}, combine = c))))
	    cyclePerYr <- unique(round(unlist(rlapply(popL, level = 3, f = function(x) {x[["paramL"]][["cyclePerYr"]]}, combine = c))))

	    simStats <- rlapply(popL, function(x) {x[!names(x) %in% c("SP", "paramL", "VgVDP", "gvVDP")]}, level = 3)
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
		f1names <- names(popL[[1]]) 

		# varCoord <- list(x = 1:length(f1names), z = yr, y= varL[[i]])
		# RGSCcoord <- list(x = 1:length(f1names), z = RGSCgen, y= gsL[[i]])

		# ouput list here. then plot / lines...
	}

	lapply(popList, formatPop)

	ylimits <- range(0, unlist(c(varL, gsL)))

	
	# pdf("truncSel.pdf")
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
	mtext("Recurrent GS intensity", 1,  cex = 1, padj = 2.5, outer = TRUE)
	mtext("Genetic Value of Varieties", 2,  cex = 1, adj = 0.5, padj = -2, outer = TRUE)

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
	mtext("Recurrent GS intensity", 1,  cex = 1, padj = 2.5, outer = TRUE)
	mtext("Genetic Value of Recurrent Population", 2,  cex = 1, adj = 0.5, padj = -2, outer = TRUE)

	# dev.off()

}





























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