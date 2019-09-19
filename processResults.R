# popL = popList; depth = 2; meanVariety = TRUE; removeErr = TRUE

formatPop <- function(popL, depth = 1, meanVariety = TRUE, removeErr = TRUE){
	get1 <- function(l, what, depth, counter = 0) if(depth == counter) return(l[[what]]) else get1(l[[1]], what, depth, counter = counter + 1)
	get1names <- function(depth, l, counter = 0) if(depth == counter) return(names(l)) else get1names(depth, l[[1]], counter = counter + 1)
	getVgRatio<- function(x){
		varg <- if(length(x[["VgRGSC"]]) > x[["paramL"]][["nYr"]] + 1) x[["VgRGSC"]][seq(1, x[["paramL"]][["nYr"]]*x[["paramL"]][["cyclePerYr"]], x[["paramL"]][["cyclePerYr"]])] else x[["VgRGSC"]][-length(x[["VgRGSC"]])]
		x[["VgVDP"]][["trial1"]] / varg
	}

	getSratio2 <- function(x){
		gvRS <- if(length(x[["gvRGSC"]]) > x[["paramL"]][["nYr"]] + 1) x[["gvRGSC"]][seq(1, x[["paramL"]][["nYr"]]*x[["paramL"]][["cyclePerYr"]], x[["paramL"]][["cyclePerYr"]])] else x[["gvRGSC"]][-length(x[["gvRGSC"]])]
		s <- {x[["gvVDP"]][["variety"]] - x[["gvVDP"]][["trial1"]]} / {x[["gvVDP"]][["variety"]] - gvRS}
		s[s == Inf] <- NA
		s
	}

	getSratio<- function(x)	{
		s <- x[["sVDP"]] / x[["sTotal"]]
		s[s == Inf] <- NA
		s
	}

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

    nVar <- unique(round(unlist(rlapply(popL, level = depth, f = function(x) {tail(with(x[["paramL"]], nFam * famSize * cumprod(selectTrials)), 1)}, combine = c))))
    cyclePerYr <- unique(round(unlist(rlapply(popL, level = depth, f = function(x) {x[["paramL"]][["cyclePerYr"]]}, combine = c))))
    pull <- unique(unlist(rlapply(popL, level = depth, f = function(x) {x[["paramL"]][["pullCycle"]]}, combine = c)))
    if(is.null(pull)) pull <- cyclePerYr
    nTrial <- unique(round(unlist(rlapply(popL, level = depth, f = function(x) {length(x[["paramL"]][["selectTrials"]])}, combine = c))))
    nYr <- unique(round(unlist(rlapply(popL, level = depth, f = function(x) {x[["paramL"]][["nYr"]]}, combine = c))))
    trad <- sapply(rlapply(popL, level = depth, f = function(x) {x[["paramL"]][["traditional"]]}, combine = c), unique)

    pulldiff <- cyclePerYr - pull

    if (length(nYr) > 1) stop("number of years dont match!")

    vdpStats <- rlapply(popL, function(x) {x[names(x) %in% c("VgVDP", "gvVDP", "VDPacc")]}, level = depth)
    VgVDPtrial1 <- rlapply( vdpStats, function(x) {x$VgVDP$trial1}, level = depth, combine = rbind)
	VgVDPtrial1 <- lapply(VgVDPtrial1, colMeans)

    simStats <- rlapply(popL, function(x) {x[!names(x) %in% c("SP", "paramL", "VgVDP", "gvVDP", "VDPacc")]}, level = depth)
    simStatsInv <- rlapply(simStats, invertList, level = depth - 1)
    simReps <- rlapply(simStatsInv, level = depth + 1, combine = rbind) 
    simReps <- rlapply(simReps, level = depth - 1, f = function(x) x[!sapply(x, is.null)]) 
    simAvg <- rlapply(simReps, f = colMeans, level = depth, na.rm = TRUE)

    simSd <- rlapply(simReps, f = function(x) {apply(x, 2, sd, na.rm = TRUE)}, level = depth)
    # simSd[[1]]$varMean

    # if (meanVariety) {
    # 	simAvg <- rlapply(simAvg, function(x) {x[["vy"]] <- x[["varMean"]]; x[["vx"]] <- x[["RGSCyr"]]; x}, level = depth - 1)
    # }

    yr <- 1:nYr
	varL <- rlapply(simAvg, "[[", level = depth - 1, i = "varMean")
	varSd <- rlapply(simSd, "[[", level = depth - 1, i = "varMean")
	RGSCyr <- rlapply(simAvg, "[[", level = depth - 1, i = "RGSCyr")
	RGSCgen <- rlapply(simAvg, "[[", level = depth - 1, i = "Rcyc")
	RGSCacc <- rlapply(simAvg, "[[", level = depth - 1, i = "RGSCacc")
	VDPinAcc <- rlapply(simAvg, "[[", level = depth - 1, i = "VDPinAcc")
	RGSCoutAcc <- rlapply(simAvg, "[[", level = depth - 1, i = "RGSCoutAcc")
	gsL <- rlapply(simAvg, "[[", level = depth - 1, i = "gvRGSC")
	Vg <- rlapply(simAvg, "[[", level = depth - 1, i = "VgRGSC")
	sVDPL <- rlapply(simAvg, "[[", level = depth - 1, i = "sVDP")
	iVDPL <- rlapply(simAvg, "[[", level = depth - 1, i = "iVDP")
	sTotalL <- rlapply(simAvg, "[[", level = depth - 1, i = "sTotal")
	iTotalL <- rlapply(simAvg, "[[", level = depth - 1, i = "iTotal")

	sVDPtoSTotal
	for(i in 1:length(popL)) iVDPL[[i]] / iTotal[[i]]

	VgVDPtoVgRGSC <- lapply(rlapply(popL, getVgRatio, level = depth, combine = rbind), colMeans, na.rm = TRUE)
	# sVDPtoSTotal2 <- lapply(rlapply(popL, getSratio2, level = depth, combine = rbind), colMeans, na.rm = TRUE)
	sVDPtoSTotal <- lapply(rlapply(popL, getSratio, level = depth, combine = rbind), colMeans, na.rm = TRUE)

	# for(i in 1:length(varL)){
	# 	varg <- if(length(Vg[[i]]) > nYr + 1) Vg[[i]][seq(1, nYr*cyclePerYr, cyclePerYr)] else Vg[[i]][-length(Vg[[i]])]
	# 	VgVDPtoVgRGSC[[i]] <- VgVDPtrial1[[i]] / varg
	# }


	# sVDPtoSTotal <- list()
	# for(i in 1:length(varL)){
	# 	sVDPtoSTotal[[i]] <- sVDPL[[i]] / sTotalL[[i]]
	# }

	list(nVar = nVar, cyclePerYr = cyclePerYr, nTrial = nTrial, nYr = nYr, RGSCyr = RGSCyr, RGSCgen = RGSCgen, yr = yr, 
		# xlims = xlims, ylims = ylims, 
		VDPinAcc = VDPinAcc, RGSCoutAcc = RGSCoutAcc, VgVDPtrial1 = VgVDPtrial1, VgVDPtoVgRGSC = VgVDPtoVgRGSC, sVDPtoSTotal = sVDPtoSTotal, 
		 varL = varL, varSd = varSd, Vg = Vg, gsL = gsL, RGSCacc = RGSCacc, sVDPL = sVDPL, iVDPL = iVDPL, sTotalL = sTotalL, iTotalL = iTotalL, 
		 sdVar = errors = errs)
}

parDir <- getwd()
source(paste0(parDir, "/alphaTools.R"))

# defArgs <- list(filenames = c('results/traditional/traditional30yr1000QTL_trad2_intWithin0.2_intAcross0.5_truth0_vdp10x50.RData','results/truncation/truncation30yr1000QTL_trunc_truth0_rgsc0.1_vdp10x50.RData', 'results/phenoRGSC/phenoRGSC30yr1000QTL_quadprog_fin0.005_fout0.1_N1_pull0_phRS0.6_sepTrn0_truth0_rgsc0.2_vdp10x50.RData'), newfilename = "test")
# defArgs <- list(filenames = c('results/traditional/traditional30yr1000QTL_trad2_intWithin0.2_intAcross0.5_truth0_vdp10x50.RData'), newfilename = NULL)

defArgs <- list(filenames = NULL, newfilename = NULL)
defArgs <- getComArgs(defArgs)

if(is.null(defArgs[["filenames"]])) stop("please provide a filenames!")

attach(defArgs)

experiment <- gsub("/.*", "", gsub("results/", "", filenames))
vdp <- gsub(".*vdp|\\.RData", "", filenames)

if(length(filenames) > 1 & is.null(newfilename)) stop("multiple files, please provide an aregument to newfilename!")
# projDir <- sapply(strsplit(filenames, "/"), "[[", 2)	
# if(is.null(newfilename)) newfilename <- paste0(parDir, "/resultSummary/", projDir, "/", gsub(".*/", "",filenames))
if(is.null(newfilename)) newfilename <- paste0(parDir, "/", gsub("results/", "resultSummary/", filenames))

popList <- list()
for(i in 1:length(filenames)){
	# load(paste0(parDir, "/results/", projDir, simNameij, "/", simNameij, ".RData"))
	load(filenames[i])
	cat(filenames[i], "has", length(simrun), "founder pops with", length(simrun[[1]]), "reps...\n")
	popList[[i]] <- simrun 
	nullsim <- sapply(simrun, is.null)
	if(any(nullsim)) {
		cat(simNameij, "has", sum(nullsim), "missing replicates!!!\n")
		system(paste0("cat ", simNameij, "has", sum(nullsim), "missing replicates!!! >> summaryErrors.txt"))
		q()
	}
}
names(popList) <- gsub("results/", "", filenames)

popList <- rlapply(popList, level = 2, f = unlist, recursive = FALSE)
statList <- formatPop(popList, depth = 2)
save(statList, file = newfilename)
