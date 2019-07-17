# popL = popList; depth = 4; meanVariety = TRUE; removeErr = TRUE

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
		classx <- lapply(x, class)
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

	# class(popL$`1000QTL`$quadprog_fin0.005_fout0.2$truth1$vdp4x25$rgsc0.2[[1]])
	# popL$`1000QTL`$quadprog_fin0.005_fout0.2$truth1$vdp4x25$rgsc0.2[[1]]
	# popL$`1000QTL`$quadprog_fin0.005_fout0.2$truth1$vdp4x25$rgsc0.2[[50]]
	# x <- popL$`1000QTL`$quadprog_fin0.005_fout0.2$truth1$vdp4x25$rgsc0.2 
# names(popL[[1]][[1]][[1]][[1]][[1]])


	# parm <- rlapply(popList, level = 6, f= function(x) {x[["paramL"]][["nFam"]]}, combine = list)
	# classes <- rlapply(popL, level = depth, f= class, combine = c)

	if(removeErr) errs <- rlapply(popL, level = depth - 1, f= getErr)
	if(removeErr) errsIndex <- rlapply(popL, level = depth - 1, f= whichErr)
	if(removeErr) popL <- rlapply(popL, level = depth - 1, f= rmErrList)
	# popL$`1000QTL`$quadprog_fin0.01_fout0.3_N1_pull0$truth2$vdp10x50$rgsc0.2[c(8, 30, 43)]


	# parm <- rlapply(popL, level = 6, f= function(x) {is.null(x[["paramL"]])}, combine = c)
	# popL[["1000QTL"]]$trunc$truth0$vdp40x200$rgsc0.2[parm[["1000QTL"]]$trunc$truth0$vdp40x200$rgsc0.2]
 
	# problems <- which(unlist(rlapply(parm, level = depth -1, f = any, combine = c)))


    nVar <- unique(round(unlist(rlapply(popL, level = depth, f = function(x) {tail(with(x[["paramL"]], nFam * famSize * cumprod(selectTrials)), 1)}, combine = c))))
    cyclePerYr <- unique(round(unlist(rlapply(popL, level = depth, f = function(x) {x[["paramL"]][["cyclePerYr"]]}, combine = c))))

    simStats <- rlapply(popL, function(x) {x[!names(x) %in% c("SP", "paramL", "VgVDP", "gvVDP", "VDPacc")]}, level = depth)
   
    # invert list is failing here. dont know why....
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
	RGSCacc <- rlapply(simAvg, "[[", level = depth - 1, i = "RGSCacc", rbind)
	VDPinAcc <- rlapply(simAvg, "[[", level = depth - 1, i = "VDPinAcc", rbind)
	RGSCoutAcc <- rlapply(simAvg, "[[", level = depth - 1, i = "RGSCoutAcc", rbind)
	# varL <- rlapply(simAvg, "[[", level = depth - 1, i = "predAcc", rbind)
	gsL <- rlapply(simAvg, "[[", level = depth - 1, i = "gvRGSC", rbind)
	Vg <- rlapply(simAvg, "[[", level = depth - 1, i = "VgRGSC", rbind)
	SL <- rlapply(simAvg, "[[", level = depth - 1, i = "sVDP", rbind)
	iL <- rlapply(simAvg, "[[", level = depth - 1, i = "iVDP", rbind)

	# lapply(iL, rowMeans)
	# getIntensity
	# factorNames <- lapply(0:(depth - 2), get1names, l = simAvg)
	# names(factorNames) <- names(factors)

	list(nVar = nVar, cyclePerYr = cyclePerYr, RGSCyr = RGSCyr, RGSCgen = RGSCgen, yr = yr, xlims = xlims, ylims = ylims, VDPinAcc = VDPinAcc, RGSCoutAcc = RGSCoutAcc,
		 varL = varL, Vg = Vg, gsL = gsL, RGSCacc = RGSCacc, SL = SL, iL = iL, errors = errs)
}

# popL <- popList[[1]]; depth = 6; fLabs = NULL; plotReps = FALSE; plotVg = TRUE; plotSelInt = TRUE



# truncVSquadprog20yr1000QTL_quadprog_fin0.005_fout0.2_truth2_rgsc0.2_vdp40x200  
# truncVSquadprog20yr1000QTL_trunc_truth1_rgsc0.2_vdp10x50
# truncVSquadprog20yr1000QTL
# i = "1000QTL"
# j = "quadprog_fin0.005_fout0.2"
# k = "truth2"
# l = "rgsc0.2"
# m = "vdp40x200"

parDir <- getwd()
source(paste0(parDir, "/alphaTools.R"))

# defArgs <- list(QTL = c("1000"),#, "100QTL", "10QTL"),
# 				select = c("quadprog"),
# 				fin = c(0.005, 0.01, 0.05),
# 				fout = c(0.1),
# 				N = c(1),
# 				pull = c(0),
# 				truth = c("0", "1", "2"),
# 				crossSel = c("rgsc0.2"), 
# 				VDP = c( '10x50','20x100'),
# 				simName = "truncVSquadprog20yr", 
# 				projDir = "truncVSquadprog/", 
# 				figDir = "figures/temp"
# 				)

defArgs <- list(QTL = c("1000"),#, "100QTL", "10QTL"),
				select = c("quadprog"),
				fin = c("0.01"),
				fout = c("0.1"),
				N = c("1"),
				pull = c(1),
				truth = c("0", "1", "2"),
				crossSel = c("rgsc0.2"), 
				VDP = c("10x50"),
				simName = "truncVSquadprog20yr", 
				projDir = "truncVSquadprog/", 
				figDir = "figures/temp"
				)

byVDP <- TRUE
defArgs <- getComArgs(defArgs)


defArgs$QTL <- paste0(defArgs$QTL, "QTL")
defArgs$fin[defArgs$fin != "" & !is.na(defArgs$fin)] <- paste0("fin", defArgs$fin[defArgs$fin != "" & !is.na(defArgs$fin)])
defArgs$fout[defArgs$fout != "" & !is.na(defArgs$fout)] <- paste0("fout", defArgs$fout[defArgs$fout != "" & !is.na(defArgs$fout)])
defArgs$N[defArgs$N != "" & !is.na(defArgs$N)] <- paste0("N", defArgs$N[defArgs$N != "" & !is.na(defArgs$N)])
defArgs$pull[defArgs$pull != "" & !is.na(defArgs$pull)] <- paste0("pull", defArgs$pull[defArgs$pull != "" & !is.na(defArgs$pull)])
defArgs$truth <- paste0("truth", defArgs$truth)
defArgs$VDP <- paste0("vdp", defArgs$VDP)

attach(defArgs)
# argList <- as.list(commandArgs(TRUE))
# argList <- c("truncSel_truncCross" "truncSel_expDistPairs" "truncSel_simDHdistPairs")

factors <- defArgs[c("QTL", "select", "fin", "fout", "N", "pull", "truth", "crossSel", "VDP")]
factl <- sapply(factors, length)
vary <- which(factl > 1)
vary <- vary[!names(vary) %in% c("truth", "VDP")]

if(length(vary) == 0) vary <- which(names(factors) == "select") 

popList <- list()
for(i in factors[[vary]]){
	facti <- factors
	facti[[vary]] <- i
	for(j in truth) {
		factij <- facti
		factij[["truth"]] <- j
		if(factij[["pull"]][[1]] == "pull3") {factij[["fout"]] <- factij[["N"]] <- ""}
		for(k in VDP){	
			factijk <- factij
			factijk[["VDP"]] <- k
			simNameij <- gsub("[_]+", "_", paste0(c(simName, unlist(factijk)), c("", "_", "_", "_", "_", "_", "_", "_", "_", ""), collapse = ""))
			load(paste0(parDir, "/results/", projDir, simNameij, "/", simNameij, ".RData"))
			popList[[i]][[j]][[k]] <- simrun 
			nullsim <- sapply(simrun, is.null)
			if(any(nullsim)) cat(simNameij, "has", sum(nullsim), "missing replicates!!!\n")
		}
	}
}

facti <- factors
for(j in truth) {
	factij <- facti
	factij[["truth"]] <- j

	factij[["select"]] <- "trunc"
	for(k in VDP){	
		factijk <- factij
		factijk[["VDP"]] <- k

		simNameTruncj <- gsub("[_]+", "_", paste0(c(simName, unlist(factijk[c("QTL", "select", "truth", "crossSel", "VDP")])), c("",  "_", "_", "_", "_", ""), collapse = ""))

		load(paste0(parDir, "/results/", projDir, simNameTruncj, "/", simNameTruncj, ".RData"))
		popList[["trunc"]][[j]][[k]] <- simrun 
		nullsim <- sapply(simrun, is.null)
		if(any(nullsim)) cat(simNameTruncj, "has", sum(nullsim), "missing replicates!!!\n")
	}
}


# popList <- list()
# for(i in QTL){
# 	for(j in select[!select %in% "trunc"]){
# 		for(k in fin[!(fin %in% "" | is.na(fin))]){
# 			for(l in fout[!(fout %in% "" | is.na(fout))]){
# 				for(m in N[!(N %in% "" | is.na(N))]){
# 					for(n in pull[!(pull %in% "" | is.na(pull))]){
# 						for(o in truth){
# 							for(p in crossSel){
# 								for(q in VDP){
# 									if(n == "pull3") {l <- m <- ""}
# 									# simNameijklm <- paste0(c(simName, c(i, j, k, l, m)), c("", "", "_", "", "_", ""), collapse = "")
# 									simNameijklmnopq <- gsub("[_]+", "_", paste0(c(simName, c(i, j, k, l, m, n, o , p, q)), c("", "_", "_", "_", "_", "_", "_", "_", "_", ""), collapse = ""))
# 									load(paste0(parDir, "/results/", projDir, simNameijklmnopq, "/", simNameijklmnopq, ".RData"))
# 									if(byVDP) popList[[i]][[j]][[k]][[m]][[l]][[n]][[o]][[p]][[q]] <- simrun 
# 									nullsim <- sapply(simrun, is.null)
# 									if(any(nullsim)) cat(simNameijklmnopq, "has", sum(nullsim), "missing replicates!!!\n")
# 								}
# 							}
# 						}
# 					}
# 				}
# 			}
# 		}
# 	}
# }

# # load trunc
# for(i in QTL){
# 	for(j in select[select %in% "trunc"]){
# 		for(n in truth){
# 			for(o in crossSel){
# 				for(q in VDP){
# 					simNameijnop <- gsub("[_]+", "_", paste0(c(simName, c(i, j, n, o, q)), c("", "_", "_", "_", "_", ""), collapse = ""))
# 					if(j == "trunc") {k <- l <- m <- "trunc"}
# 					load(paste0(parDir, "/results/", projDir, simNameijnop, "/", simNameijnop, ".RData"))
# 					if(byVDP) popList[[i]][[j]][[k]][[m]][[l]][[n]][[o]][[p]][[q]] <- simrun 
# 					nullsim <- sapply(simrun, is.null)
# 					if(any(nullsim)) cat(simNameijnop, "has", sum(nullsim), "missing replicates!!!\n")
# 				}
# 			}
# 		}
# 	}
# }


# cols <- c("#006600", "#00008C", "#660000", "#8C8C00", "#660066", "#006666")
# cols <- cols[1:length(popList[[1]])]


# popListHOLD <- popList
# popList <- popListHOLD

popList <- rlapply(popList, level = 3, f = unlist, recursive = FALSE)
statList <- formatPop(popList, depth = 4)

getYranges <- function(x, indexL, j, k) {
	ll <- list()
	for(i in indexL) ll <- c(ll, x[[i]][[j]][k,])
	range(unlist(ll)) * c(0.9, 1.1)
}

figDir <- paste0(parDir, "/", figDir)
system(paste0("simdir=", figDir, "
if [ ! -d $simdir ]; then
	mkdir -p $simdir
else 
	echo 'Directory already exists! Any results therein will be overwritten without further notice...'
fi"))

# select <- names(statList[["varL"]])
select = c("trunc", names(statList[["varL"]])[names(statList[["varL"]]) != "trunc"])
# select = c("quadprog_fin0.01_pull3", "quadprog_fin0.01_fout0.1_N1_pull0", "quadprog_fin0.01_fout0.2_N1_pull0", "quadprog_fin0.01_fout0.3_N1_pull0")

# for(i in names(statList$varL[[1]])){ # "trunc" "quadprog_fin0.005_fout0.2" "quadprog_fin0.01_fout0.2" 
for(j in names(statList$varL[[1]])){ # truth0  truth1 truth2
	for(k in rownames(statList$varL[[1]][[1]])){ # VDP
		
		pdfName <- paste0(figDir, "/", paste(c("truncVSquadprog", names(factors)[vary], j, k), collapse = "_"), "_", paste(factors[sapply(factors, length) == 1], collapse = "_"), ".pdf") 
		pdf(pdfName, width = 12, height = 7)

		ylims <- getYranges(statList$varL, select, j, k)
		plot(statList$RGSCyr, statList$varL[[select[1]]][[j]][k,], type = "l", ylim = ylims, xlab = "year", ylab = "Variety Mean", main = "Variety Means", xaxt = "n")
	    axis(1, at = c(0, statList$RGSCyr), labels = c(0, statList$yr))
		for(s in 2:length(select)) lines(statList$RGSCyr, statList$varL[[select[s]]][[j]][k,], lty = s)
		legend("topleft", legend = select, lty = 1:length(select))

		ylims <- getYranges(statList$gsL, select, j, k)
		plot(statList$RGSCgen, statList$gsL[[select[1]]][[j]][k,], type = "l", ylim = ylims, xlab = "year", ylab = "Recurrent Population Mean", main = "Recurrent Population Mean", xaxt = "n")
	    axis(1, at = c(0, statList$RGSCyr), labels = c(0, statList$yr))
		for(s in 2:length(select)) lines(statList$RGSCgen, statList$gsL[[select[s]]][[j]][k,], lty = s)
		legend("topleft", legend = select, lty = 1:length(select))

		ylims <- getYranges(statList$Vg, select, j, k)
		plot(statList$RGSCgen, statList$Vg[[select[1]]][[j]][k,], type = "l", ylim = ylims, xlab = "year", main = "Recurrent Population Variance", ylab = "Variance", xaxt = "n")
	    axis(1, at = c(0, statList$RGSCyr), labels = c(0, statList$yr))
		for(s in 2:length(select)) lines(statList$RGSCgen, statList$Vg[[select[s]]][[j]][k,], lty = s)
		legend("topright", legend = select, lty = 1:length(select))

		ylims <- c(0, 1)
		plot(statList$RGSCgen[-length(statList$RGSCgen)], statList$RGSCacc[[select[1]]][[j]][k,], type = "l", ylim = ylims, xlab = "year",  main = "Recurrent Population Prediction Accuracy", ylab = "Prediction Accuracy", xaxt = "n")
	    axis(1, at = c(0, statList$RGSCyr), labels = c(0, statList$yr))
		for(s in 2:length(select)) lines(statList$RGSCgen[-length(statList$RGSCgen)], statList$RGSCacc[[select[s]]][[j]][k,], lty = s)
		legend("topright", legend = select, lty = 1:length(select))

		if(k %in% rownames(statList$RGSCoutAcc[[select[1]]][[j]])){
			ylims <- c(0, 1)
			plot(statList$RGSCyr, statList$RGSCoutAcc[[select[1]]][[j]][k,], type = "l", ylim = ylims, xlab = "year",  main = "Prediction Accuracy of material out of Recurrent Population", ylab = "Prediction Accuracy", xaxt = "n")
		    axis(1, at = c(0, statList$RGSCyr), labels = c(0, statList$yr))
			for(s in 2:length(select)) lines(statList$RGSCyr, statList$RGSCoutAcc[[select[s]]][[j]][k,], lty = s)
			legend("topright", legend = select, lty = 1:length(select))
		}
		if(k %in% rownames(statList$VDPinAcc[[select[1]]][[j]])){
			ylims <- c(0, 1)
			plot(statList$RGSCyr, statList$VDPinAcc[[select[1]]][[j]][k,], type = "l", ylim = ylims, xlab = "year",  main = "Prediction Accuracy of material into VDP", ylab = "Prediction Accuracy", xaxt = "n")
		    axis(1, at = c(0, statList$RGSCyr), labels = c(0, statList$yr))
			for(s in 2:length(select)) lines(statList$RGSCyr, statList$VDPinAcc[[select[s]]][[j]][k,], lty = s)
			legend("topright", legend = select, lty = 1:length(select))

		}
		dev.off()	
	}
}


