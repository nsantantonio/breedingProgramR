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
    trad <- unique(round(unlist(rlapply(popL, level = depth, f = function(x) {x[["paramL"]][["traditional"]]}, combine = c))))

    simStats <- rlapply(popL, function(x) {x[!names(x) %in% c("SP", "paramL", "VgVDP", "gvVDP", "VDPacc")]}, level = depth)
   
    # invert list is failing here. dont know why....
    simStatsInv <- rlapply(simStats, invertList, level = depth - 1)
    simReps <- rlapply(simStatsInv, level = depth + 1, combine = rbind) 
    simReps <- rlapply(simReps, level = depth - 1, f = function(x) x[!sapply(x, is.null)]) 

    simAvg <- rlapply(simReps, f = colMeans, level = depth, na.rm = TRUE)

    RGSCyr <- get1(simAvg, "RGSCyr", depth - 1)
    # if(all(trad > 0)) 
    RGSCgen <- get1(simAvg, "Rcyc", depth - 1)
    yr <- RGSCyr / cyclePerYr
    xlims <- range(c(0, RGSCgen))
    ylims <- range(unlist(rlapply(simReps, getYrange, level = depth - 1))) * 1.1

    if (meanVariety) {
    	simAvg <- rlapply(simAvg, function(x) {x[["vy"]] <- x[["varMean"]]; x[["vx"]] <- x[["RGSCyr"]]; x}, level = depth - 1)
    }

	varL <- rlapply(simAvg, "[[", level = depth - 1, i = "vy")
	RGSCacc <- rlapply(simAvg, "[[", level = depth - 1, i = "RGSCacc")
	VDPinAcc <- rlapply(simAvg, "[[", level = depth - 1, i = "VDPinAcc")
	RGSCoutAcc <- rlapply(simAvg, "[[", level = depth - 1, i = "RGSCoutAcc")
	gsL <- rlapply(simAvg, "[[", level = depth - 1, i = "gvRGSC")
	Vg <- rlapply(simAvg, "[[", level = depth - 1, i = "VgRGSC")
	SL <- rlapply(simAvg, "[[", level = depth - 1, i = "sVDP")
	iL <- rlapply(simAvg, "[[", level = depth - 1, i = "iVDP")

	list(nVar = nVar, cyclePerYr = cyclePerYr, RGSCyr = RGSCyr, RGSCgen = RGSCgen, yr = yr, xlims = xlims, ylims = ylims, VDPinAcc = VDPinAcc, RGSCoutAcc = RGSCoutAcc,
		 varL = varL, Vg = Vg, gsL = gsL, RGSCacc = RGSCacc, SL = SL, iL = iL, errors = errs)
}

parDir <- getwd()
source(paste0(parDir, "/alphaTools.R"))


defArgs <- list(filenames = c('results/traditional/traditional30yr1000QTL_trad2_intWithin1_intAcross1_truth0_vdp20x75.RData', 'results/traditional/traditional30yr1000QTL_trad3_intWithin1_intAcross1_truth0_vdp20x75.RData'),
figDir = "figures/test", figName = "test.pdf", labels = NULL)
defArgs <- getComArgs(defArgs)

if(is.null(defArgs[["labels"]])) defArgs[["labels"]] <- 1:length(defArgs[["filenames"]])
attach(defArgs)

popList <- list()
for(i in 1:length(filenames)){
	# load(paste0(parDir, "/results/", projDir, simNameij, "/", simNameij, ".RData"))
	load(filenames[i])
	cat(filenames[i], "has", length(simrun), "founder pops with", length(simrun[[1]]), "reps...\n")
	popList[[i]] <- simrun 
	nullsim <- sapply(simrun, is.null)
	if(any(nullsim)) cat(simNameij, "has", sum(nullsim), "missing replicates!!!\n")
}

length(popList)
length(popList[[1]])
length(popList[[1]][[1]])
length(popList[[1]][[1]][[1]])
popList <- rlapply(popList, level = 2, f = unlist, recursive = FALSE)

statList <- formatPop(popList, depth = 2)

getYranges <- function(x) range(unlist(x)) * c(0.9, 1.1)


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
pdfName <- paste0(figDir, "/", figName) 
pdf(pdfName, width = 12, height = 7)

ylims <- getYranges(statList$varL)
plot(statList$RGSCyr, statList$varL[[1]], type = "l", ylim = ylims, xlab = "year", ylab = "Variety Mean", main = "Variety Means", xaxt = "n")
axis(1, at = c(0, statList$RGSCyr), labels = c(0, statList$yr))
for(s in 2:length(filenames)) lines(statList$RGSCyr, statList$varL[[s]], lty = s)
legend("topleft", legend = labels, lty = 1:length(filenames))
dev.off()	

# ylims <- getYranges(statList$gsL, select, j, k)
# plot(statList$RGSCgen, statList$gsL[[select[1]]][[j]][k,], type = "l", ylim = ylims, xlab = "year", ylab = "Recurrent Population Mean", main = "Recurrent Population Mean", xaxt = "n")
# axis(1, at = c(0, statList$RGSCyr), labels = c(0, statList$yr))
# for(s in 2:length(select)) lines(statList$RGSCgen, statList$gsL[[select[s]]][[j]][k,], lty = s)
# legend("topleft", legend = select, lty = 1:length(select))

# ylims <- getYranges(statList$Vg, select, j, k)
# plot(statList$RGSCgen, statList$Vg[[select[1]]][[j]][k,], type = "l", ylim = ylims, xlab = "year", main = "Recurrent Population Variance", ylab = "Variance", xaxt = "n")
# axis(1, at = c(0, statList$RGSCyr), labels = c(0, statList$yr))
# for(s in 2:length(select)) lines(statList$RGSCgen, statList$Vg[[select[s]]][[j]][k,], lty = s)
# legend("topright", legend = select, lty = 1:length(select))

# ylims <- c(0, 1)
# plot(statList$RGSCgen[-length(statList$RGSCgen)], statList$RGSCacc[[select[1]]][[j]][k,], type = "l", ylim = ylims, xlab = "year",  main = "Recurrent Population Prediction Accuracy", ylab = "Prediction Accuracy", xaxt = "n")
# axis(1, at = c(0, statList$RGSCyr), labels = c(0, statList$yr))
# for(s in 2:length(select)) lines(statList$RGSCgen[-length(statList$RGSCgen)], statList$RGSCacc[[select[s]]][[j]][k,], lty = s)
# legend("topright", legend = select, lty = 1:length(select))

# if(k %in% rownames(statList$RGSCoutAcc[[select[1]]][[j]])){
# 	ylims <- c(0, 1)
# 	plot(statList$RGSCyr, statList$RGSCoutAcc[[select[1]]][[j]][k,], type = "l", ylim = ylims, xlab = "year",  main = "Prediction Accuracy of material out of Recurrent Population", ylab = "Prediction Accuracy", xaxt = "n")
#     axis(1, at = c(0, statList$RGSCyr), labels = c(0, statList$yr))
# 	for(s in 2:length(select)) lines(statList$RGSCyr, statList$RGSCoutAcc[[select[s]]][[j]][k,], lty = s)
# 	legend("topright", legend = select, lty = 1:length(select))
# }
# if(k %in% rownames(statList$VDPinAcc[[select[1]]][[j]])){
# 	ylims <- c(0, 1)
# 	plot(statList$RGSCyr, statList$VDPinAcc[[select[1]]][[j]][k,], type = "l", ylim = ylims, xlab = "year",  main = "Prediction Accuracy of material into VDP", ylab = "Prediction Accuracy", xaxt = "n")
#     axis(1, at = c(0, statList$RGSCyr), labels = c(0, statList$yr))
# 	for(s in 2:length(select)) lines(statList$RGSCyr, statList$VDPinAcc[[select[s]]][[j]][k,], lty = s)
# 	legend("topright", legend = select, lty = 1:length(select))

# }
# 	}
# }


