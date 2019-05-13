# set.seed(12345)
 
# pop = RGSC[[lastRGSCgen]]; GSfit = GSmodel[[lastGSmodel]];
# save(pop, GSfit, SP, file = "tempPopGSfitGen6.RData")
# initFunc <- rank

# x = ebv(pop); n = 20; targetFunc = popQuant; pherFunc = pherFunc; evapRate = 0.05; nAnts = 500; maxIter = 300; countThresh = 1000; showResult = TRUE; pherPower = 1.5; dumbAnt = FALSE; plateau = 100; returnStats = FALSE
acOpt <- function(x, n, targetFunc, pherFunc, xAt0 = FALSE, evapRate = 0.05, nAnts = 500, maxIter = 300, countThresh = 1000, showResult = FALSE, pherPower = 1, dumbAnt = FALSE, plateau = 100, returnStats = FALSE){
	evap <- function(oldP, newP, evapRate) {(1 - evapRate) * oldP + newP }
	# xindex <- 1:length(x)
	# xord <- order(x)
	# if(!identical(xord, xindex) & !identical(xord, rev(xindex))) stop("x is not sorted!")
	if(is.matrix(x)) {if(ncol(x) == 1) x <- c(x) else stop("I cant handle multiple traits yet!")}
	
	if(xAt0){
		minx <- min(x)
		x <- x - minx
	}
	N <- length(x)
	pheromone <- initFunc(x)

	bestAnt <- NULL
	noP = rep(0, N)
	shortestPath = 0
	pathCounter = 0 
	iter = 0
	if(showResult) {
		plotBest <- NULL 
		plotMean <- NULL 
	}

	while(pathCounter <= countThresh & iter < maxIter){
	# while(!identical(optSol, sort(bestAnt)) & iter < maxIter){
		iter = iter + 1
		ant <- list()
		path <- list()
		antPheromone <- list()
		for(i in 1:nAnts){
			ant[[i]] <- sample(1:N, n, prob = pheromone / sum(pheromone))
			path[[i]] <- targetFunc(sel = ant[[i]], ebvs = x)
			# path[[i]] <- quant(x[ant[[i]]], sigma = 2, w = 0.5)
			Pi <- noP
			Pi[ant[[i]]] <- path[[i]]
			antPheromone[[i]] <- Pi
		}
		bestPath <- which.max(path)
		if(path[[bestPath]] > shortestPath) {
			shortestPath <- path[[bestPath]] 
			bestAnt <- ant[[bestPath]]
			pathCounter <- 0
		} else {
			pathCounter <- pathCounter + 1
		}
		newPheromone <- pherFunc(Reduce("+", antPheromone), p = pherPower)
		
		if(!dumbAnt) pheromone <- evap(pheromone, newPheromone, evapRate = evapRate)
		
		if(showResult) {
			plotBest <- c(plotBest, path[[bestPath]])
			plotMean <- c(plotMean, mean(unlist(path)))
			ylims <- range(c(plotBest, plotMean))
			cols <- rep("gray80", N)
			cols[ant[[bestPath]]] <- "firebrick"
			cols[bestAnt] <- "black"
			par(mfrow = c(2, 1))
			xord <- order(x)
			barplot(pheromone[xord], col = cols[xord])
			plot(1:iter, plotBest, col = "firebrick", type = "b", pch = 16, ylim = ylims)
			lines(1:iter, plotMean, col = "black", type = "b", pch = 16)
			if(iter > plateau + 3){
				abline(h = mean(plotMean[plateau:iter]), col = "black", lty = 2)
				abline(lm({plotMean[plateau:iter]} ~ {plateau:iter}), col = "black")
			}
			# print(iter)
			# Sys.sleep(0.05)
		}
	}
	# sort(bestAnt)
	if(returnStats){
		return(list(iter = 1:iter, best = plotBest, mean = plotMean, pheromone = pheromone, bestAnt = sort(bestAnt), ylims = ylims))
	} else {
		return(bestAnt)
	}
}



load("tempPopGSfitGen6.RData")
library(AlphaSimR)
source("alphaTools.R")
# x <- sort(rnorm(100), decreasing = TRUE)
# ebvs <- ebv(pop)

initFunc <- function(x) rep(1, length(x))

quant <- function(x, sigma, w) {2 * (w * mean(x) + (1 - w) * sigma * sd(x))}
# quant <- function(x, sigma) {mean(x) + sigma * sd(x)}
popQuant <- function(sel, ebvs, sigma = 2, w = 0.5) {quant(x = ebvs[sel], sigma = sigma, w = 0.5)}

pherFuncID <- function(x, p) {x}
pherFunc <- function(x, p) {(x / max(x))^p}

gsQuant <- function(sel, ebvs, sigma = 2, w = 0.5, nCrosses = 100, nProgenyPerCross = 1) {
	if(!all(ebv(pop[sel]) == ebvs[sel])) stop("ebvs of pop dont match those of ant!")
	simCrosses <- randomCross(pop[sel], nFam = nCrosses, nProgeny = nProgenyPerCross)
	simPop <- makeCross(pop[sel], simCrosses)
	simPop <- setEBV(simPop, GSfit)
	quant(ebv(simPop), sigma = sigma, w = w)
}

truncPop <- function(x, n, returnSel = TRUE) {
	if(is.matrix(x)) x <- c(x)
	xindex <- 1:length(x)
	xord <- order(x, decreasing = TRUE)
	if(returnSel) x[xord][1:n] else xindex[xord][1:n]
} 

getI <- function(sel, pop){
	(mean(sel) - mean(ebv(pop)) )/ varA(pop)
}

# have target function make the cross then eval pop mean / or best ?
# plot.new()

ebvs <- ebv(pop)

selTrunc <- truncPop(ebvs, 20, FALSE)
mean(ebvs[selTrunc])
popQuant(selTrunc, ebvs, sigma = 2, w = 0.5)
getI(ebvs[selTrunc], pop)
varA(pop[selTrunc]) - varA(pop)

set.seed(1234)
mu <- list()
sigmasq <- list()
pq <- list()
int <- list()
for(i in 1:10){
	sol1 <- acOpt(ebv(pop), n = 20, xAt0 = TRUE, targetFunc = popQuant, pherFunc = pherFunc, showResult = TRUE, evapRate = 0.05, nAnts = 500, pherPower = 1.5)
	selsol1 <- ebvs[sol1]
	mu[[i]] <- mean(ebvs[sol1])
	pq[[i]] <- popQuant(sol1, ebv(pop))
	int[[i]] <- getI(ebvs[sol1], pop)
	sigmasq[[i]] <- varA(pop[sol1])
}
# mean(unlist(lapply(sigmasq, function(x) x- varA(pop))))



> varA(pop[sol1]) - varA(pop)
             [,1]                                               
[1,] 0.0001187818


set.seed(1234)
sol2 <- acOpt(ebv(pop), n = 20, xAt0 = FALSE, targetFunc = popQuant, pherFunc = pherFunc, showResult = TRUE, evapRate = 0.05, nAnts = 500, pherPower = 1.5)
popQuant(x[sol2])

set.seed(1234)
sol <- acOpt(ebv(pop), n = 20, targetFunc = gsQuant, pherFunc = pherFunc, showResult = TRUE, evapRate = 0.05, nAnts = 200, pherPower = 1.25, maxIter = 1000)

# sol <- acOpt(sort(ebv(pop)), n = 20, targetFunc = popQuant, pherFunc = pherFunc, showResult = TRUE, evapRate = 0.05, nAnts = 500, pherPower = 1.5)


# N <- 100
# n <- 20
# load(founderRData)
# founderRData <- "founderPop/testAlphaSimR1000SegSite.RData"
# founderPop <- founderPop[sample(1:nInd(founderPop), N)]
# pop <- founderPop
# GS <- RRBLUP(pop)
# mu <- list()
# sigmasq <- list()
# pq <- list()
# int <- list()
# for(i in nGen){

# 	sol1 <- acOpt(ebv(pop), n = 20, xAt0 = TRUE, targetFunc = popQuant, pherFunc = pherFunc, showResult = TRUE, evapRate = 0.05, nAnts = 500, pherPower = 1.5)
# 	selsol1 <- ebvs[sol1]
# 	mu[[i]] <- mean(ebvs[sol1])
# 	pq[[i]] <- popQuant(sol1, ebv(pop))
# 	int[[i]] <- getI(ebvs[sol1], pop)
# 	sigmasq[[i]] <- varA(pop[sol1])

# 	mergePop(pheno, pop[sel])
# }




acSol <- list()
for(i in c(0.01, 0.05, 0.1)){
		print(paste0("evap: ", i))
		acSol[["id"]][[paste0("evap", i)]] <- acOpt(ebvs, n = 20, 
													targetFunc = popQuant, 
													pherFunc = pherFuncID, 
													showResult = TRUE, 
													evapRate = i, 
													nAnts = 500, 
													returnStats = TRUE)
	for(j in c(1, 2)){
		print(paste0("    pherP: ", j))
		acSol[["percMax"]][[paste0("evap", i)]][[paste0("pherP", j)]] <- acOpt(x, n = 20,
																		  targetFunc = popQuant, 
																		  pherFunc = pherFunc, 
																		  showResult = TRUE, 
																		  evapRate = i, 
																		  nAnts = 500, 
																		  returnStats = TRUE,
																		  pherPower = j)
	}
}


# acSol2 <- list()
for(i in c(0.005, 0.01, 0.025, 0.05, 0.1)){
	print(paste0("evap: ", i))
	for(j in c(1, 1.5, 2, 3)){
		print(paste0("    pherP: ", j))
		acSol2[["percMax"]][[paste0("evap", i)]][[paste0("pherP", j)]] <- acOpt(x, n = 20,
																		  targetFunc = popQuant, 
																		  pherFunc = pherFunc, 
																		  showResult = TRUE, 
																		  evapRate = i, 
																		  nAnts = 1000, 
																		  returnStats = TRUE,
																		  pherPower = j, 
																		  maxIter = 500)
	}
}

names(acSol2)

acSol[["percMax"]][[1]][[1]][["best"]]

pheromoneList <- rlapply(acSol[[1]], "[[", level = 2, i = "pheromone") 
bestAntList <- rlapply(acSol[[1]], "[[", level = 2, i = "bestAnt") 
iterlist <- rlapply(acSol[[1]], "[[", level = 2, i = "iter")
bestlist <- rlapply(acSol[[1]], "[[", level = 2, i = "best")
meanlist <- rlapply(acSol[[1]], "[[", level = 2, i = "mean")

ylims <- range(c(bestlist, meanlist))
par(mfrow = c(4, 5), mar = c(0, 0.5, 0.5, 0), oma = c(4, 5, 0.5, 0.5), tcl = -0.25, mgp = c(2, 0.6, 0))
for(j in c(1, 1.5, 2, 3)){
	for(i in c(0.005, 0.01, 0.025, 0.05, 0.1)){
			cols <- rep("gray80", length(ebvs))
			cols[bestAntList[[i]][[j]]] <- "black"
			barplot(pheromoneList[[i]][[j]], col = cols)
			# plot(iterlist[[i]][[j]], bestlist[[i]][[j]], col = "firebrick", type = "b", pch = 16, ylim = ylims)
			# lines(iterlist[[i]][[j]], meanlist[[i]][[j]], col = "black", type = "b", pch = 16)
	}
}

dev.new()
ylims <- range(c(bestlist, meanlist))
par(mfrow = c(4, 5), mar = c(0, 0.5, 0.5, 0), oma = c(4, 5, 0.5, 0.5), tcl = -0.25, mgp = c(2, 0.6, 0))
for(j in c(1, 1.5, 2, 3)){
	for(i in c(0.005, 0.01, 0.025, 0.05, 0.1)){
			# cols <- rep("gray80", length(ebvs))
			# cols[bestAntList[[i]][[j]]] <- "black"
			# barplot(pheromoneList[[i]][[j]], col = cols)
			plot(iterlist[[i]][[j]], bestlist[[i]][[j]], col = "firebrick", type = "b", pch = 16, ylim = ylims)
			lines(iterlist[[i]][[j]], meanlist[[i]][[j]], col = "black", type = "b", pch = 16)
	}
}


rlapply(acSol, level = 3, length, rbind)

# save(acSol, file = "tempACsol.RData")


# acSol <- acOpt(x, n = 20, targetFunc = popQuant, pherFunc = pherFuncID, showResult = TRUE, evapRate = 0.05, nAnts = 500, returnStats = TRUE)
# acSol <- acOpt(x, n = 20, targetFunc = popQuant, pherFunc = pherFuncID, showResult = TRUE, evapRate = 0.5, nAnts = 500, returnStats = TRUE)
# acSol <- acOpt(x, n = 20, targetFunc = popQuant, pherFunc = pherFuncID, showResult = TRUE, evapRate = 0.9, nAnts = 500, returnStats = TRUE)


# acSol <- acOpt(x, n = 20, targetFunc = popQuant, pherFunc = pherFunc, showResult = TRUE, evapRate = 0.05, nAnts = 500, pherPower = 1)
# acSol <- acOpt(x, n = 20, targetFunc = popQuant, pherFunc = pherFunc, showResult = TRUE, evapRate = 0.5, nAnts = 500, pherPower = 1)
# acSol <- acOpt(x, n = 20, targetFunc = popQuant, pherFunc = pherFunc, showResult = TRUE, evapRate = 0.9, nAnts = 500, pherPower = 1)

# acSol <- acOpt(x, n = 20, targetFunc = popQuant, pherFunc = pherFunc, showResult = TRUE, evapRate = 0.05, nAnts = 500, pherPower = 2)
# acSol <- acOpt(x, n = 20, targetFunc = popQuant, pherFunc = pherFunc, showResult = TRUE, evapRate = 0.5, nAnts = 500, pherPower = 2)
# acSol <- acOpt(x, n = 20, targetFunc = popQuant, pherFunc = pherFunc, showResult = TRUE, evapRate = 0.9, nAnts = 500, pherPower = 2)



solID <- acOpt(x, n = 20, targetFunc = quant, pherFunc = pherFuncID, showResult = TRUE, evapRate = 0.05, nAnts = 500, returnStats = TRUE)
solID2 <- acOpt(x, n = 20, targetFunc = quant, pherFunc = pherFuncID, showResult = TRUE, evapRate = 0.9, nAnts = 500, returnStats = TRUE)



sol1 <- acOpt(x, n = 20, targetFunc = quant, pherFunc = pherFunc, showResult = TRUE, evapRate = 0.05, nAnts = 500, pherPower = 1, returnStats = TRUE)
sol1.1 <- acOpt(x, n = 20, targetFunc = quant, pherFunc = pherFunc, showResult = TRUE, evapRate = 0.05, nAnts = 500, pherPower = 1.25, returnStats = TRUE, maxIter = 1000)
max(sol1.25$best)


sol1.25 <- acOpt(x, n = 20, targetFunc = quant, pherFunc = pherFunc, showResult = TRUE, evapRate = 0.05, nAnts = 500, pherPower = 1.25, returnStats = TRUE, maxIter = 500)
max(sol1.25$best)
sol1.5 <- acOpt(x, n = 20, targetFunc = quant, pherFunc = pherFunc, showResult = TRUE, evapRate = 0.05, nAnts = 500, pherPower = 1.5, returnStats = TRUE)
max(sol1.5$best)

sol1.5 <- acOpt(x, n = 20, targetFunc = quant, pherFunc = pherFunc, showResult = TRUE, evapRate = 0.05, nAnts = 500, pherPower = 2, returnStats = TRUE)

soldumb <- acOpt(x, n = 20, targetFunc = quant, pherFunc = pherFunc, showResult = TRUE, evapRate = 0.05, nAnts = 500, pherPower = 2, returnStats = TRUE, dumbAnt = TRUE)


sol <- acOpt(x, n = 20, targetFunc = quant, pherFunc = pherFunc, showResult = TRUE, evapRate = 0.05, nAnts = 500, pherPower = 2)
acOpt(x, n = 20, targetFunc = popQuant, pherFunc = pherFunc, showResult = TRUE, evapRate = 0.5, nAnts = 500, pherPower = 2)
acOpt(x, n = 20, targetFunc = popQuant, pherFunc = pherFunc, showResult = TRUE, evapRate = 0.9, nAnts = 500, pherPower = 2)





randSol	<- acOpt(x, n = 20, targetFunc = quant, dumbAnt = TRUE, showResult = TRUE)


reps <- 10
result <- list()
for(i in 1:reps){
	print(i)
	x <- sort(rnorm(100), decreasing = TRUE)
	result[[i]] <- c(acSol, randSol)
}
rmat <- do.call(rbind, result)
rmat[, 1] / rmat[, 2] - 1
acOpt(x = ebv(pop), n = 20, showResult = TRUE)