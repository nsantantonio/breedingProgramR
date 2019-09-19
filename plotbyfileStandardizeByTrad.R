# popL = popList; depth = 2; meanVariety = TRUE; removeErr = TRUE

formatPop <- function(popL, depth = 1, meanVariety = TRUE, removeErr = TRUE){
	get1 <- function(l, what, depth, counter = 0) if(depth == counter) return(l[[what]]) else get1(l[[1]], what, depth, counter = counter + 1)
	get1names <- function(depth, l, counter = 0) if(depth == counter) return(names(l)) else get1names(depth, l[[1]], counter = counter + 1)
	rmErr <- function(x) {
		classx <- lapply(x, class)
		err <- sapply(classx, function(x) "error" %in% x)
		lapply(x[err], print)
		x[!err]
	}

	getErr <- function(x) {
		classx <- sapply(x, class)
		err <- sapply(classx, function(x) "error" %in% x)
		x[err]
	}
	whichErr <- function(x) {
		classx <- lapply(x, class)
		err <- which(sapply(classx, function(x) "error" %in% x))
	}

	rmErrList <- function(x) {
		classx <- lapply(x, class)
		noterr <- sapply(classx, function(x) "list" %in% x)
		x[noterr]
	}


	rmShort <- function(x, len = 10) {
		lenxnotok <- sapply(x, length) < len
		if(any(lenxnotok)) msg(1, "WARNING! ", sum(lenxnotok), " simulations failed! Returning", length(x) - sum(lenxnotok), "sucessful simulations")
		x[!lenxnotok]
	}

	if(removeErr) errs <- rlapply(popL, level = depth - 1, f= getErr)
	if(removeErr) errsIndex <- rlapply(popL, level = depth - 1, f= whichErr)
	if(removeErr) popL <- rlapply(popL, level = depth - 1, f= rmErrList)
	if(removeErr) popL <- rlapply(popL, level = depth - 1, f= rmShort)

	# z <- NULL
    # for(i in 1:length(popL[[1]])) z[i] <- popL[[1]][[i]][["paramL"]][["nFam"]]
    nVar <- unique(round(unlist(rlapply(popL, level = depth, f = function(x) {tail(with(x[["paramL"]], nFam * famSize * cumprod(selectTrials)), 1)}, combine = c))))
    cyclePerYr <- unique(round(unlist(rlapply(popL, level = depth, f = function(x) {x[["paramL"]][["cyclePerYr"]]}, combine = c))))
    nTrial <- unique(round(unlist(rlapply(popL, level = depth, f = function(x) {length(x[["paramL"]][["selectTrials"]])}, combine = c))))
    nYr <- unique(round(unlist(rlapply(popL, level = depth, f = function(x) {x[["paramL"]][["nYr"]]}, combine = c))))
    trad <- sapply(rlapply(popL, level = depth, f = function(x) {x[["paramL"]][["traditional"]]}, combine = c), unique)

    if (length(nYr) > 1) stop("number of years dont match!")

    simStats <- rlapply(popL, function(x) {x[!names(x) %in% c("SP", "paramL", "VgVDP", "gvVDP", "VDPacc")]}, level = depth)
   
    # invert list is failing here. dont know why....
    simStatsInv <- rlapply(simStats, invertList, level = depth - 1)
    simReps <- rlapply(simStatsInv, level = depth + 1, combine = rbind) 
    simReps <- rlapply(simReps, level = depth - 1, f = function(x) x[!sapply(x, is.null)]) 

    simAvg <- rlapply(simReps, f = colMeans, level = depth, na.rm = TRUE)

	# RGSCyr <- lapply(rlapply(simReps, f = "[[", level = depth - 1, i = "RGSCyr"), unique)
    
    # RGSCyr <- get1(simAvg, "RGSCyr", depth - 1)
    # if(all(trad > 0)) 
    # RGSCgen <- get1(simAvg, "Rcyc", depth - 1)
    # yr <- RGSCyr / cyclePerYr
    yr <- 1:nYr
    # xlims <- range(c(0, RGSCgen))
    # ylims <- range(unlist(rlapply(simReps, getYrange, level = depth - 1))) * 1.1

    if (meanVariety) {
    	simAvg <- rlapply(simAvg, function(x) {x[["vy"]] <- x[["varMean"]]; x[["vx"]] <- x[["RGSCyr"]]; x}, level = depth - 1)
    }

	varL <- rlapply(simAvg, "[[", level = depth - 1, i = "vy")
	RGSCyr <- rlapply(simAvg, "[[", level = depth - 1, i = "vx")
	RGSCgen <- rlapply(simAvg, "[[", level = depth - 1, i = "Rcyc")
	RGSCacc <- rlapply(simAvg, "[[", level = depth - 1, i = "RGSCacc")
	VDPinAcc <- rlapply(simAvg, "[[", level = depth - 1, i = "VDPinAcc")
	RGSCoutAcc <- rlapply(simAvg, "[[", level = depth - 1, i = "RGSCoutAcc")
	gsL <- rlapply(simAvg, "[[", level = depth - 1, i = "gvRGSC")
	Vg <- rlapply(simAvg, "[[", level = depth - 1, i = "VgRGSC")
	SL <- rlapply(simAvg, "[[", level = depth - 1, i = "sVDP")
	iL <- rlapply(simAvg, "[[", level = depth - 1, i = "iVDP")

	list(nVar = nVar, cyclePerYr = cyclePerYr, nTrial = nTrial, nYr = nYr, RGSCyr = RGSCyr, RGSCgen = RGSCgen, yr = yr, 
		# xlims = xlims, ylims = ylims, 
		VDPinAcc = VDPinAcc, RGSCoutAcc = RGSCoutAcc,
		 varL = varL, Vg = Vg, gsL = gsL, RGSCacc = RGSCacc, SL = SL, iL = iL, errors = errs)
}

plotPopVar <- function(x, y, s, vLine = "none", popcol = "#000000", alpha = "", alphaMean = "33", ...) {
	polycol <- paste0(popcol, alphaMean)
	popcol <- paste0(popcol, alpha)
	xpoly <- c(x, rev(x), x[1])
	ypoly <- c(y + s, rev(y - s), y[1] + s[1])

	polygon(x = xpoly, y = ypoly, col = polycol, border = NA)
	lines(x = x, y = y, type = "l", col = popcol, ...)
}

parDir <- getwd()
source(paste0(parDir, "/alphaTools.R"))


defArgs <- list(vdp = c("10x50", "20x75", "40x100"), 
	tradfilenames = c('results/traditionalBest/traditionalBest30yr1000QTL_trad2_intWithin0.2_intAcross0.5_truth0_vdp'),
	filenames=c( 'results/truncation/truncation30yr1000QTL_trunc_truth0_rgsc0.3_vdp',
                        'results/quadprog/quadprog30yr1000QTL_fin0.005_pull3_truth0_rgsc0.2_vdp', 
                        'results/quadprog/quadprog30yr1000QTL_fin0.005_fout0.1_N1_pull0_truth0_rgsc0.2_vdp',
                        'results/phenoRGSC/phenoRGSC30yr1000QTL_quadprog_fin0.005_fout0.1_N1_pull0_phRS0.6_sepTrn0_truth0_rgsc0.2_vdp'),
    figDir="figures/bestModels",
    figName="best.pdf",
    labels = c('truncation', 'optimalContribution', 'optContBranch', 'optContBranch_PhenoRS'),
    whichfiles = NULL,
    printYears = seq(5, 30, 5), 
    cols = NULL,
    fromSummary = TRUE,
    lineType = NULL,
    applyColor = TRUE,
    lineWeight = 1.5
    )

# defArgs <- list(vdp = c("10x50", "20x75", "40x100"), 
# 	tradfilenames = c('results/traditional/traditional30yr1000QTL_trad2_intWithin0.2_intAcross0.5_truth0_vdp'),
# 	filenames = c('resultSummary/truncation/truncation30yr1000QTL_trunc_truth0_rgsc0.1_vdp'),
#     figDir = "figures/test",
#     figName = "test.pdf",
#     labels = c('trunc_phRS0.4'),
#     whichfiles = NULL,
#     printYears = seq(5, 30, 5), 
#     cols = c("#000000", "#FF0000"),
#     fromSummary = TRUE)

defArgs <- list(vdp = NULL, tradfilenames = NULL, filenames = NULL, figDir = "figures/test", figName = "test.pdf", labels = NULL, whichfiles = NULL, printYears = seq(5, 30, 5), cols = NULL, lineType = 1, applyColor = FALSE, fromSummary = TRUE)
defArgs <- getComArgs(defArgs)


defCols <- c("#70a845", "#c45ca2", "#4aad92", "#cb584c", "#7878cd", "#b88f40")[c(4, 5, 6, 1, 2, 3)]

if(is.null(defArgs[["filenames"]])) stop("please provide filenames!")
if(is.null(defArgs[["lineType"]])) defArgs[["lineType"]] <- 0:length(defArgs[["filenames"]]) + 1 
if(is.null(defArgs[["labels"]])) defArgs[["labels"]] <- 1:length(defArgs[["filenames"]])
# if(is.null(defArgs[["whichfiles"]])) defArgs[["whichfiles"]] <- 1:length(defArgs[["filenames"]])
if(is.null(defArgs[["printYears"]])) defArgs[["printYears"]] <- FALSE

attach(defArgs)


if(fromSummary)	{
	filenames <- gsub("results/", "resultSummary/", filenames)
	tradfilenames <- gsub("results/", "resultSummary/", tradfilenames)
} else {
	if(all(gsub("/.*", "", filenames) == "resultSummary") & all(gsub("/.*", "", tradfilenames) == "resultSummary")) fromSummary <- TRUE
}

if(!is.null(defArgs[["vdp"]])) {
	tradfilenames <- paste0(tradfilenames, vdp, ".RData")
	filenames <- c(sapply(filenames, paste0, vdp, ".RData"))
}

if(is.null(whichfiles)) whichfiles <- 1:length(defArgs[["filenames"]])

tradvdp <- gsub(".*vdp|\\.RData", "", tradfilenames)
vdp <- gsub(".*vdp|\\.RData", "", filenames)

nVdp <- length(tradvdp)
nExp <- length(vdp) / nVdp
nf <- length(filenames)
uniqLabs <- labels
if(length(labels) != nf & length(labels) %% 1 == 0) labels <- rep(labels, each = nVdp)
if (nExp %% 1 != 0) stop("wrong number of filenames to match trad")
if (!all(rep(tradvdp, nExp) == vdp)) stop("vdps dont match!")
labelsLong <- paste0(labels, "_", vdp)

popList <- list()
for(i in 1:length(filenames)){
	# load(paste0(parDir, "/results/", projDir, simNameij, "/", simNameij, ".RData"))
	load(filenames[i])
	if(fromSummary) {
		popList[[i]] <- statList
	} else {
		cat(filenames[i], "has", length(simrun), "founder pops with", length(simrun[[1]]), "reps...\n")
		popList[[i]] <- simrun 
		nullsim <- sapply(simrun, is.null)
		if(any(nullsim)) cat(simNameij, "has", sum(nullsim), "missing replicates!!!\n")
	}
}

tradList <- list()
for(i in 1:length(tradfilenames)){
	# load(paste0(parDir, "/results/", projDir, simNameij, "/", simNameij, ".RData"))
	load(tradfilenames[i])
	if(fromSummary) {
		tradList[[i]] <- statList
	} else {
		cat(tradfilenames[i], "has", length(simrun), "founder pops with", length(simrun[[1]]), "reps...\n")
		tradList[[i]] <- simrun 
		nullsim <- sapply(simrun, is.null)
		if(any(nullsim)) cat(simNameij, "has", sum(nullsim), "missing replicates!!!\n")
	}
}


# length(popList)
# length(popList[[1]])
# length(popList[[1]][[1]])
# length(popList[[1]][[1]][[1]])

if(fromSummary){
	tradStatList <- rlapply(invertList(rep(tradList, nExp)), level = 1, f = unlist, recursive = FALSE)
	statList <- rlapply(invertList(popList), level = 1, f = unlist, recursive = FALSE)
} else {
	tradList <- rlapply(tradList, level = 2, f = unlist, recursive = FALSE)
	popList <- rlapply(popList, level = 2, f = unlist, recursive = FALSE)
	statList <- formatPop(popList, depth = 2)
	tradStatList <- formatPop(rep(tradList, nExp), depth = 2)
}

# for(i in 1:length(popList)){
# "SP"         "paramL"     "Rcyc"       "varMean"    "sdRGSC"    
# "VgRGSC"     "VgVDP"      "gvRGSC"     "gvVDP"      "sRGSC"     
# "iRGSC"      "sVDP"       "iVDP"       "sTotal"     "iTotal"    
# "nVar"       "vx"         "vy"         "RGSCyr"     "RGSCacc"   
# "VDPacc"     "RGSCoutAcc" "VDPinAcc"   "theorMax" 
# }

getYranges <- function(x) range(unlist(x)) * c(0.95, 1.1)


figDir <- paste0(parDir, "/", figDir)
system(paste0("simdir=", figDir, "
if [ ! -d $simdir ]; then
	mkdir -p $simdir
else 
	echo 'Directory already exists! Any results therein will be overwritten without further notice...'
fi"))

# select <- names(statList[["varL"]])
# select = c("trad2", "trunc", names(statList[["varL"]])[!names(statList[["varL"]]) %in% c("trad2","trunc")])
# select = c("trad", names(statList[["varL"]])[!names(statList[["varL"]]) %in% c("trad")])
# select = c("quadprog_fin0.01_pull3", "quadprog_fin0.01_fout0.1_N1_pull0", "quadprog_fin0.01_fout0.2_N1_pull0", "quadprog_fin0.01_fout0.3_N1_pull0")

# for(i in names(statList$varL[[1]])){ # "trunc" "quadprog_fin0.005_fout0.2" "quadprog_fin0.01_fout0.2" 


if(!is.null(cols) & is.numeric(cols)) cols <- defCols[cols]  # to reorder / select default colors
if(applyColor & is.null(cols)) {
	cols <- defCols[0:nExp + 1]
	vdpCols <- defCols[1:nVdp]
} else {
	cols <- rep("#000000", nExp + 1)
	vdpCols <- rep("#000000", nVdp)
}

varStd <- list()
for(i in 1:length(statList$varL)) varStd[[i]] <- statList$varL[[i]] / tradStatList$varL[[i]] 

xvar <- statList$RGSCyr[[which.max(sapply(statList$RGSCyr, max))]]
ylimsvar <- getYranges(varStd)
ylimsvar2 <- getYranges(c(tradStatList$varL, statList$varL))
xrgsc <- statList$RGSCgen
gens <- sapply(xrgsc, length)
xrgsc[gens == min(gens)] <- lapply(xrgsc[gens == min(gens)], function(x) x * max(gens - 1) / min(gens - 1))
ylimsgs <- getYranges(statList$gsL)
varShift <- statList$nTrial*statList$cyclePerYr
nYr <- statList$nYr
ylimsVg <- getYranges(statList$Vg)
ylimsGS <- getYranges(statList$gsL)

printYr <- if(printYears[[1]]) printYears * unique(statList$cyclePerYr) else print("boo")

if(length(lineType) == 1) ltyVDP <- lineType else ltyVDP <- lineType[1:nVdp]

# whichfiles <- 7:9

# # varieties relative to traditional
# plot(NA, type = "l", ylim = ylimsvar, xlim = range(xvar), xlab = "year", ylab = "Variety Mean", main = "Variety Means", xaxt = "n", col = cols[1])
# abline(1, 0, lty = 1, col = "gray80")
# axis(1, at = c(0, xvar), labels = c(0, statList$yr))
# l <- 1
# for(i in whichfiles) lines(xvar, varStd[[i]], lty = i); l <- l + 1
# legend("topleft", legend = labels[whichfiles], lty = 1:length(whichfiles))

# # varieties against traditional

# plot(NA, type = "l", ylim = ylimsvar2, xlim = ylimsvar2, xlab = "", ylab = "", col = cols[1])
# abline(0, 1, col = "gray80")
# # axis(1, at = statList$yr, labels = c(0, statList$yr))
# l <- 1
# for(i in whichfiles) {
# 	lines(tradStatList$varL[[i]], statList$varL[[i]], lty = i)
# 	l <- l + 1
# }

# legend("topleft", legend = labels[whichfiles], lty = 1:nf)

index <- split(1:nf, rep(1:{nf / nVdp}, each = nVdp))


pdfName <- paste0(figDir, "/", figName) 
pdf(pdfName, width = 15, height = 5)


par(mfrow = c(1, nVdp), mar = c(0, 0, 0, 0), oma = c(4, 5, 0.5, 0.5), tcl = -0.25, mgp = c(2, 0.6, 0))
v <- "varL"
ylims <- getYranges(c(tradStatList[[v]], statList[[v]]))

for(i in 1:nVdp){
	pry <- if(i == 1) 's' else 'n'
	plot(xvar, tradStatList[[v]][[index[[1]][[i]]]], type = "l", ylim = ylims, xlim = range(xvar), xlab = "", ylab = "", xaxt = "n", yaxt = pry, col = cols[1], lwd = lineWeight)
	# axis(1, at = c(0, xvar), labels = c(0, statList$yr))
	axis(1, at = c(0, printYr), labels = c(0, printYears))
	for(k in 1:nExp) lines(xvar, statList[[v]][[index[[k]][i]]], lty = lineType[k + 1], col = cols[k + 1], lwd = lineWeight)
	if(i == 1) legend("topleft", legend = c("traditional", uniqLabs), lty = lineType, col = cols, lwd = lineWeight)
	legend("topright", legend = vdp[i], bty = 'n')
}
mtext("Variety Means", side = 2, outer = TRUE, padj = -2)
mtext("Year", side = 1, outer = TRUE, padj = 2)

par(mfrow = c(1, nVdp), mar = c(0, 0, 0, 0), oma = c(4, 5, 0.5, 0.5), tcl = -0.25, mgp = c(2, 0.6, 0))
v <- "sVDPL"
ylims <- getYranges(c(tradStatList[[v]], statList[[v]]))

for(i in 1:nVdp){
	pry <- if(i == 1) 's' else 'n'
	plot(xvar, tradStatList[[v]][[index[[1]][[i]]]], type = "l", ylim = ylims, xlim = range(xvar), xlab = "", ylab = "", xaxt = "n", yaxt = pry, col = cols[1], lwd = lineWeight)
	# axis(1, at = c(0, xvar), labels = c(0, statList$yr))
	axis(1, at = c(0, printYr), labels = c(0, printYears))
	for(k in 1:nExp) lines(xvar, statList[[v]][[index[[k]][i]]], lty = lineType[k + 1], col = cols[k + 1], lwd = lineWeight)
	if(i == 1) legend("topleft", legend = c("traditional", uniqLabs), lty = lineType, col = cols, lwd = lineWeight)
	legend("topright", legend = vdp[i], bty = 'n')
}
mtext("Selection differential in VDP", side = 2, outer = TRUE, padj = -2)
mtext("Year", side = 1, outer = TRUE, padj = 2)


par(mfrow = c(1, nVdp), mar = c(0, 0, 0, 0), oma = c(4, 5, 0.5, 0.5), tcl = -0.25, mgp = c(2, 0.6, 0))
v <- "sTotalL"
ylims <- getYranges(c(tradStatList[[v]], statList[[v]]))

for(i in 1:nVdp){
	pry <- if(i == 1) 's' else 'n'
	plot(xvar, tradStatList[[v]][[index[[1]][[i]]]], type = "l", ylim = ylims, xlim = range(xvar), xlab = "", ylab = "", xaxt = "n", yaxt = pry, col = cols[1], lwd = lineWeight)
	# axis(1, at = c(0, xvar), labels = c(0, statList$yr))
	axis(1, at = c(0, printYr), labels = c(0, printYears))
	for(k in 1:nExp) lines(xvar, statList[[v]][[index[[k]][i]]], lty = lineType[k + 1], col = cols[k + 1], lwd = lineWeight)
	if(i == 1) legend("topleft", legend = c("traditional", uniqLabs), lty = lineType, col = cols, lwd = lineWeight)
	legend("topright", legend = vdp[i], bty = 'n')
}
mtext("Total selection differential", side = 2, outer = TRUE, padj = -2)
mtext("Year", side = 1, outer = TRUE, padj = 2)


par(mfrow = c(1, nVdp), mar = c(0, 0, 0, 0), oma = c(4, 5, 0.5, 0.5), tcl = -0.25, mgp = c(2, 0.6, 0))
v <- "iVDPL"
ylims <- getYranges(c(tradStatList[[v]], statList[[v]]))

for(i in 1:nVdp){
	pry <- if(i == 1) 's' else 'n'
	plot(xvar, tradStatList[[v]][[index[[1]][[i]]]], type = "l", ylim = ylims, xlim = range(xvar), xlab = "", ylab = "", xaxt = "n", yaxt = pry, col = cols[1], lwd = lineWeight)
	# axis(1, at = c(0, xvar), labels = c(0, statList$yr))
	axis(1, at = c(0, printYr), labels = c(0, printYears))
	for(k in 1:nExp) lines(xvar, statList[[v]][[index[[k]][i]]], lty = lineType[k + 1], col = cols[k + 1], lwd = lineWeight)
	if(i == 1) legend("topleft", legend = c("traditional", uniqLabs), lty = lineType, col = cols, lwd = lineWeight)
	legend("topright", legend = vdp[i], bty = 'n')
}
mtext("Selection intensity in VDP", side = 2, outer = TRUE, padj = -2)
mtext("Year", side = 1, outer = TRUE, padj = 2)


par(mfrow = c(1, nVdp), mar = c(0, 0, 0, 0), oma = c(4, 5, 0.5, 0.5), tcl = -0.25, mgp = c(2, 0.6, 0))
v <- "iTotalL"
ylims <- getYranges(c(tradStatList[[v]], statList[[v]]))

for(i in 1:nVdp){
	pry <- if(i == 1) 's' else 'n'
	plot(xvar, tradStatList[[v]][[index[[1]][[i]]]], type = "l", ylim = ylims, xlim = range(xvar), xlab = "", ylab = "", xaxt = "n", yaxt = pry, col = cols[1], lwd = lineWeight)
	# axis(1, at = c(0, xvar), labels = c(0, statList$yr))
	axis(1, at = c(0, printYr), labels = c(0, printYears))
	for(k in 1:nExp) lines(xvar, statList[[v]][[index[[k]][i]]], lty = lineType[k + 1], col = cols[k + 1], lwd = lineWeight)
	if(i == 1) legend("topleft", legend = c("traditional", uniqLabs), lty = lineType, col = cols, lwd = lineWeight)
	legend("topright", legend = vdp[i], bty = 'n')
}
mtext("Total selection intensity", side = 2, outer = TRUE, padj = -2)
mtext("Year", side = 1, outer = TRUE, padj = 2)


par(mfrow = c(1, nVdp), mar = c(0, 0, 0, 0), oma = c(4, 5, 0.5, 0.5), tcl = -0.25, mgp = c(2, 0.6, 0))
v <- "sVDPtoSTotal"
ylims <- getYranges(c(tradStatList[[v]], statList[[v]]))

for(i in 1:nVdp){
	pry <- if(i == 1) 's' else 'n'
	plot(xvar, tradStatList[[v]][[index[[1]][[i]]]], type = "l", ylim = ylims, xlim = range(xvar), xlab = "", ylab = "", xaxt = "n", yaxt = pry, col = cols[1], lwd = lineWeight)
	# axis(1, at = c(0, xvar), labels = c(0, statList$yr))
	axis(1, at = c(0, printYr), labels = c(0, printYears))
	for(k in 1:nExp) lines(xvar, statList[[v]][[index[[k]][i]]], lty = lineType[k + 1], col = cols[k + 1], lwd = lineWeight)
	if(i == 1) legend("topleft", legend = c("traditional", uniqLabs), lty = lineType, col = cols, lwd = lineWeight)
	legend("topright", legend = vdp[i], bty = 'n')
}
mtext("Ratio of VDP selection differential to total selection differential", side = 2, outer = TRUE, padj = -2)
mtext("Year", side = 1, outer = TRUE, padj = 2)


par(mfrow = c(1, nVdp), mar = c(0, 0, 0, 0), oma = c(4, 5, 0.5, 0.5), tcl = -0.25, mgp = c(2, 0.6, 0))
v <- "gsL"
ylims <- getYranges(statList[[v]])

for(i in 1:nVdp){
	pry <- if(i == 1) 's' else 'n'
	plot(NA, type = "l", ylim = ylims, xlim = range(unlist(xrgsc)), xlab = "", ylab = "", xaxt = "n", yaxt = pry, col = cols[1], lwd = lineWeight)
	axis(1, at = c(0, printYr), labels = c(0, printYears))
	for(k in 1:nExp) plotPopVar(x = xrgsc[[index[[k]][i]]], y = statList[[v]][[index[[k]][i]]], sqrt(statList[["Vg"]][[index[[k]][i]]]), popcol = cols[k + 1], lty = lineType[k + 1], lwd = lineWeight)
	if(i == 1) legend("topleft", legend = uniqLabs, lty = 1:nExp + 1)
	legend("topright", legend = vdp[i], bty = 'n')
}
mtext("Recurrent population mean", side = 2, outer = TRUE, padj = -2)
mtext("Year", side = 1, outer = TRUE, padj = 2)


par(mfrow = c(1, nVdp), mar = c(0, 0, 0, 0), oma = c(4, 5, 0.5, 0.5), tcl = -0.25, mgp = c(2, 0.6, 0))
v <- "Vg"
ylims <- getYranges(c(tradStatList[[v]], statList[[v]]))

for(i in 1:nVdp){
	pry <- if(i == 1) 's' else 'n'
	plot(c(0, xvar), tradStatList[[v]][[index[[1]][[i]]]], type = "l", ylim = ylims, xlim = range(unlist(xrgsc)), xlab = "", ylab = "", xaxt = "n", yaxt = pry, col = cols[1], lwd = lineWeight)
	axis(1, at = c(0, printYr), labels = c(0, printYears))
	for(k in 1:nExp) lines(x = xrgsc[[index[[k]][i]]], y = statList[[v]][[index[[k]][i]]], lty = lineType[k + 1], col = cols[k + 1], lwd = lineWeight)
	if(i == 1) legend("topright", legend = c("traditional", uniqLabs), lty = lineType, col = cols, lwd = lineWeight)
	legend("topleft", legend = vdp[i], bty = 'n')
}
mtext("Recurrent population genetic variance", side = 2, outer = TRUE, padj = -2)
mtext("Year", side = 1, outer = TRUE, padj = 2)



par(mfrow = c(1, nVdp), mar = c(0, 0, 0, 0), oma = c(4, 5, 0.5, 0.5), tcl = -0.25, mgp = c(2, 0.6, 0))
v <- "VgVDPtrial1"
ylims <- getYranges(c(tradStatList[[v]], statList[[v]]))

for(i in 1:nVdp){
	pry <- if(i == 1) 's' else 'n'
	plot(xvar, tradStatList[[v]][[index[[1]][[i]]]], type = "l", ylim = ylims, xlim = range(xvar), xlab = "", ylab = "", xaxt = "n", yaxt = pry, col = cols[1], lwd = lineWeight)
	# axis(1, at = c(0, xvar), labels = c(0, statList$yr))
	axis(1, at = c(0, printYr), labels = c(0, printYears))
	for(k in 1:nExp) lines(xvar, statList[[v]][[index[[k]][i]]], lty = lineType[k + 1], col = cols[k + 1], lwd = lineWeight)
	if(i == 1) legend("topleft", legend = c("traditional", uniqLabs), lty = lineType, col = cols, lwd = lineWeight)
	legend("topright", legend = vdp[i], bty = 'n')
}
mtext("Ratio of VDP selection differential to total selection differential", side = 2, outer = TRUE, padj = -2)
mtext("Year", side = 1, outer = TRUE, padj = 2)


par(mfrow = c(1, nVdp), mar = c(0, 0, 0, 0), oma = c(4, 5, 0.5, 0.5), tcl = -0.25, mgp = c(2, 0.6, 0))
v <- "RGSCacc"
ylims <- c(min(c(0, unlist(statList$RGSCacc)), na.rm = TRUE), max(c(1, unlist(statList$RGSCacc)), na.rm = TRUE))

for(i in 1:nVdp){
	pry <- if(i == 1) 's' else 'n'
	plot(NA, type = "l", ylim = ylims, xlim = range(unlist(xrgsc)), xlab = "", ylab = "", xaxt = "n", yaxt = pry, col = cols[1], lwd = lineWeight)
	axis(1, at = c(0, printYr), labels = c(0, printYears))
	for(k in 1:nExp) lines(statList$RGSCgen[[index[[k]][i]]][gen(statList$RGSCgen[[index[[k]][i]]]) %in% names(statList$RGSCacc[[index[[k]][i]]])], statList$RGSCacc[[index[[k]][i]]], lty = lineType[k + 1], col = cols[k + 1], lwd = lineWeight)
	if(i == 1) legend("topright", legend = c("traditional", uniqLabs), lty = lineType, col = cols, lwd = lineWeight)
	legend("topleft", legend = vdp[i], bty = 'n')
}
mtext("Genomic prediction accuracy in recurrent population", side = 2, outer = TRUE, padj = -2)
mtext("Year", side = 1, outer = TRUE, padj = 2)


par(mfrow = c(1, nExp), mar = c(0, 0, 0, 0), oma = c(4, 5, 0.5, 0.5), tcl = -0.25, mgp = c(2, 0.6, 0))
v <- "varStd"
ylims <- getYranges(varStd)
paneLabs <- if(length(uniqLabs) == nExp) uniqLabs else unique(labels)
legLabs <- unique(vdp)


for(i in 1:nExp){
	# prx <- if(i == nExp) 's' else 'n'
	pry <- if(i == 1) 's' else 'n'
	plot(NA, type = "l", ylim = ylimsvar, xlim = range(c(0, xvar)), xlab = "year", xaxt = "n", yaxt = pry)
	abline(1, 0, lty = 1, col = "gray80")
	axis(1, at = c(0, printYr), labels = c(0, printYears))
	# axis(1, at = c(0, xvar), labels = c(0, statList$yr))
	for(k in 1:nVdp) lines(xvar, varStd[index[[i]]][[k]], lty = ltyVDP[k], col = vdpCols[k], lwd = lineWeight)
	if(i == nExp) legend("topright", legend = legLabs,  lty = ltyVDP, col = vdpCols, lwd = lineWeight)
	legend("topleft", legend = paneLabs[i], bty = 'n')
}
mtext("Year", 1, padj = 2.5, outer = TRUE)
mtext("Proportion of Traditional", 2, padj = -2, outer = TRUE)


dev.off()


# pdfName <- paste0(figDir, "/pairs_", figName) 
# pdf(pdfName, width = 15, height = 15)

# par(mfrow = c(nExp, nExp), mar = c(0, 0, 0, 0), oma = c(4, 5, 0.5, 0.5), tcl = -0.25, mgp = c(2, 0.6, 0))
# dVar <- "varL"
# lVar <- "varL"
# uVar <- "gsL"

# for(i in 1:nExp){
# 	for(j in 1:nExp){
# 			prx <- if(i == nExp) 's' else 'n'
# 			pry <- if(j == 1) 's' else 'n'
# 		if(i == j){
# 			plot(NA, type = "l", ylim = ylimsvar, xlim = range(xvar), xlab = "year", xaxt = "n", yaxt = pry)
# 			abline(1, 0, lty = 1, col = "gray80")
# 			if(i == nExp) axis(1, at = c(0, xvar), labels = c(0, statList$yr))
# 			for(k in 1:nVdp) lines(xvar, varStd[index[[i]]][[k]], lty = k)
# 			# legend("topleft", legend = labels[index[[i]]], lty = 1:nVdp)
# 		} else if(i > j) {
# 			li <- statList[[lVar]][index[[i]]]
# 			lj <- statList[[lVar]][index[[j]]]
# 			plot(NA, type = "l", ylim = ylimsvar2, xlim = ylimsvar2, xlab = "", ylab = "",xaxt = prx, yaxt = pry, col = cols[1], lwd = lineWeight)
# 			abline(0, 1, col = "gray80")
# 			# axis(1, at = statList$yr, labels = c(0, statList$yr))
# 			for(k in 1:nVdp) lines(lj[[k]], li[[k]], lty = k)
# 		} else {
# 			ui <- statList[[uVar]][index[[i]]]
# 			uj <- statList[[uVar]][index[[j]]]
# 			plot(NA, type = "l", ylim = ylimsGS, xlim = ylimsGS, xlab = "", ylab = "", xaxt = prx, yaxt = pry, col = cols[1], lwd = lineWeight)
# 			abline(0, 1, col = "gray80")
# 			# axis(1, at = statList$yr, labels = c(0, statList$yr))
# 			for(k in 1:nVdp) lines(uj[[k]], ui[[k]], lty = k)
# 		}
# 	}
# }

# dev.off()



