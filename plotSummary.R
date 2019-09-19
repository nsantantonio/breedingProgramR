parDir <- getwd()
source(paste0(parDir, "/alphaTools.R"))

f <- list.files("results/phenoRGSC")
qp <- f[grepl("quadprog", f) & grepl("phenoRGSC", f)]
p0 <- gsub(".*quadprog_", "", qp[grep("pull0", qp)])
p0df <- data.frame(fin = gsub("fin|_fout.*", "", p0), fout = gsub(".*fout|_N1.*", "", p0), phRS = gsub(".*phRS|_sepTrn.*", "", p0), vdp = gsub(".*vdp|\\.RData", "", p0))
levs <- sapply(p0df, unique)

v <- "varL"
panelYrs <- c(5, 10, 20, 30)
vdp = c("10x50", "20x75", "40x100")
models <- list(phenoRGSC = list(fin = c(0.001, 0.005, 0.01), fout = c(0.01, 0.05, 0.1,  0.2), pull = 0, phRS = c(0.2, 0.6), vdp = vdp),
	quadprog = NULL, 
	traditional = NULL, 
	traditionalBest = list(trad = 2, intWithin = 0.2, intAcross = 0.5, vdp = vdp),
	traditionalLimitI = NULL, 
	truncation = NULL)

if(any(unlist(rlapply(models, length, level = 2)) > 1)) 
combos <- lapply(models[sapply(models, length)> 1], expand.grid)



phenoRGSC <- paste0("phenoRGSC30yr1000QTL_quadprog_fin", combos$phenoRGSC$fin, "_fout", combos$phenoRGSC$fout, "_N1_pull", combos$phenoRGSC$pull, "_phRS", combos$phenoRGSC$phRS, "_sepTrn0_truth0_rgsc0.2_vdp", combos$phenoRGSC$vdp, ".RData")
traditionalBest <- paste0("traditionalBest30yr1000QTL_trad", combos$traditionalBest$trad, "_intWithin", combos$traditionalBest$intWithin, "_intAcross", combos$traditionalBest$intAcross, "_truth0_vdp", combos$phenoRGSC$vdp, ".RData")

popList <- list()
sumDirs <- list.files("resultSummary/")
for(i in sumDirs) {
	files <- list.files(paste0("resultSummary/", i))
	for(j in files){
		load(paste0("resultSummary/", i, "/", j))
		popList[[i]][[j]] <- statList
	}
}


if(any(!phenoRGSC %in% names(popList[["phenoRGSC"]]))) {
	print(phenoRGSC[!phenoRGSC %in% names(popList[["phenoRGSC"]])])
	stop("These files are missing!")
}

Y <- do.call(rbind, unlist(lapply(popList[["phenoRGSC"]][phenoRGSC], "[[", v), recursive = FALSE))
# yPanel<- yL[, panelYrs]

head(combos$phenoRGSC)

ylims <- range(Y)
phRS <- c(0.2, 0.4, 0.6)
phRS <- c(0.6)
for(p in phRS){
	pdf(paste0("figures/optContBranchPhRS/finFoutInteraction", p, ".pdf"))
	# par(mfrow = c(length(panelYrs), length(vdp)))
	par(mfrow = c(length(panelYrs), length(vdp)), mar = c(0, 0, 0, 0), oma = c(4, 5, 0.5, 0.5), tcl = -0.25, mgp = c(2, 0.6, 0))

	for(i in panelYrs){
		for(k in vdp){		
			whichOnes <- combos$phenoRGSC$vdp == k & combos$phenoRGSC$phRS == p
			x <- combos$phenoRGSC[whichOnes,]
			y <- Y[whichOnes, i]
			# plot(NA, xlim = range(x$fin), ylim = range(y), ylab = "", )
			pry <- if(k == vdp[1]) 's' else 'n'
			# prx <- if(i == panelYrs[length(panelYrs)]) 's' else 'n'
			plot(NA, xlim = range(x$fin), ylim = ylims, ylab = "", xaxt = 'n', yaxt = pry)
			if(i == panelYrs[length(panelYrs)]) axis(1, models$phenoRGSC$fin, models$phenoRGSC$fin)
			for(j in models$phenoRGSC$fout) lines(x$fin[x$fout == j], y[x$fout == j], lty = which(j == models$phenoRGSC$fout))
			if (k == vdp[1]) legend("topleft", legend = paste0("year ", i), bty = 'n', cex = 1.5)
			if (i == panelYrs[1]) legend("topright", legend = k, bty = 'n', cex = 1.5)

			# if(k == vdp[length(vdp)] & i == panelYrs[1]) legend("topleft", legend = bquote(Delta[fb] == .(models$phenoRGSC$fout)), lty = 1:length(models$phenoRGSC$fout))
			if(k == vdp[length(vdp)] & i == panelYrs[1]) legend("topleft", legend = do.call( 'expression', 
                                list(bquote(Delta[fb] == .(models$phenoRGSC$fout[[1]])),
                                bquote(Delta[fb] == .(models$phenoRGSC$fout[[2]])),
                                bquote(Delta[fb] == .(models$phenoRGSC$fout[[3]])),
                                bquote(Delta[fb] == .(models$phenoRGSC$fout[[4]])))), lty = 1:length(models$phenoRGSC$fout))
     # bquote(a == .(a))
		}
	}
	mtext(expression(Delta[f]), side = 1, outer = TRUE, padj = 1.2, cex = 2)
	mtext("Variety Mean", side = 2, outer = TRUE, padj = -1.5, cex = 2)
	dev.off()
}

