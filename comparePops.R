parDir <- getwd()
source(paste0(parDir, "/alphaTools.R"))

argList <- as.list(commandArgs(TRUE))
# argList <- c("testACquant", "testDefault")

cols <- c("#006600", "#00008C", "#660000", "#8C8C00", "#660066", "#006666")
cols <- if(length(argList) == 1) "#000000" else if(length(argList) <= 6) cols[1:length(argList)] else stop("I can only do up to 4 sims!")

popList <- list()
for(i in argList){
	simName <- i
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

# lapply(params, "[[", "useTrue")

ldiff <- llunion(params)
print(ldiff)

pdfName <- paste0(parDir, "/figures/", paste(argList, collapse = "_"), ".pdf") 
pdf(pdfName, width = 12, height = 7)
simPlot(popList, cols, varLine = "poly")
dev.off()

if(system("hostname", intern = TRUE) == "cbsurobbins.biohpc.cornell.edu") system(paste0("scp ", pdfName, " Bender:~/Dropbox/optibreedSim/figures/"))