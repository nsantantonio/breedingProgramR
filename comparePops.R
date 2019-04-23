parDir <- getwd()
source(paste0(parDir, "/alphaTools.R"))

argList <- as.list(commandArgs(TRUE))

# argList <- getComArgs()
# print(argList)
# argList <- list("selectMeanDefault", "selectMeanSingleFam") #, "selectVarDefault", "selectVarSingleFam" )
# argList <- list("RGSCintesity0.1", "RGSCintesity0.5", "RGSCintesity0.9")
# argList <- list("RGSCintesity0.1", "RGSCintesity0.3", "RGSCintesity0.5", "RGSCintesity0.7", "RGSCintesity0.9")
# argList <- list("truncSel_expDistPairs", "truncSel_truncCross", "truncSel_simDHdistPairs")
# argList <- list("truncSel_expDistPairs", "truncSel_truncCross", "truncSel_simDHdistPairs")
# argList <- list("expDist_truncCross", "truncSel_truncCross", "simDHdist_truncCross")
# argList <- list("expDist_expDistPairs", "truncSel_truncCross", "simDHdist_simDHdistPairs")
# argList <- c("RGSCintensity0.05_nNucl1000", "RGSCintensity0.01_nNucl1000", "RGSCintensity0.1_nNucl1000", "RGSCintensity0.2_nNucl1000", "RGSCintensity0.3_nNucl1000", "RGSCintensity0.5_nNucl1000", "RGSCintensity0.7_nNucl1000", "RGSCintensity0.9_nNucl1000")
# argList <- c("RGSCintensity0.05_nNucl1000", "RGSCintensity0.01_nNucl1000", "RGSCintensity0.1_nNucl1000", "RGSCintensity0.2_nNucl1000", "RGSCintensity0.3_nNucl1000", "RGSCintensity0.5_nNucl1000")


# argList <- c("returnVDPtoRGSC_trial1_0.1", "returnVDPtoRGSC_trial2_0.1", "returnVDPtoRGSC_trial3_0.2", "returnVDPtoRGSC_trial4_0.5", "returnVDPtoRGSC_trial5_1", "returnVDPtoRGSC_variety_1")

# argList <- list("RGSCintensity0.05", "RGSCintensity0.1", "RGSCintensity0.2", "RGSCintensity0.3", "RGSCintensity0.5")

# cols <- c("#006600", "#00008C", "#660000", "#000000")
cols <- c("#006600", "#00008C", "#660000", "#8C8C00", "#660066", "#006666")
# pie(rep(1, length(cols)), col = cols)
cols <- if(length(argList) == 1) "#000000" else if(length(argList) <= 6) cols[1:length(argList)] else stop("I can only do up to 4 sims!")

popList <- list()
for(i in argList){
	simName <- i#argList[[i]]
	load(paste0(parDir, "/results/", simName, "/", simName, ".RData"))
	popList[[i]] <- simrun
}

params <- lapply(popList, function(x) x[[1]][["paramL"]])
llunion <- function(ll){
	lunion <- function(l1, l2) l1[union(names(l1)[l1 %in% l2], names(l2)[l2 %in% l1])]
	inCom <- Reduce(lunion, ll)
	notInCom <- lapply(ll, function(x) x[!x %in% inCom])
	unique(unlist(lapply(notInCom, names)))
}

lapply(params, "[[", "useTrue")
lapply(params, "[[", "nNuclear")

# I dont understand why this doesnt work properly...
# ll <- lapply(params, "[", c("simName", "weight", "useTrue"))
# ll[[2]] <- ll[[2]][-2]

# llunion(ll)

ldiff <- llunion(params)
# lapply(params, "[", ldiff)
print(ldiff)

pdfName <- paste0(parDir, "/figures/", paste(argList, collapse = "_"), ".pdf") 
pdf(pdfName, width = 12, height = 7)
simPlot(popList, cols, varLine = "poly")
dev.off()


pdfName <- paste0(parDir, "/figures/", paste(argList, collapse = "_"), ".pdf") 
pdf(pdfName, width = 12, height = 7)
simPlot(popList, cols, varLine = "poly", meanVariety = TRUE)
dev.off()


if(system("hostname", intern = TRUE) == "cbsurobbins.biohpc.cornell.edu") system(paste0("scp ", pdfName, " Bender:~/Dropbox/optibreedSim/figures/"))