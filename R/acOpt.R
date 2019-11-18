#' acOpt function
#'
#' function to (do something)
#'
#' @param x [value]
#' @param n [value]
#' @param targetFunc [value]
#' @param pherFunc [value]
#' @param xAt0 [value]. Default is FALSE
#' @param evapRate [value]. Default is 0.05
#' @param nAnts [value]. Default is 500
#' @param maxIter [value]. Default is 300
#' @param countThresh [value]. Default is 1000
#' @param showResult [value]. Default is FALSE
#' @param pherPower [value]. Default is 1
#' @param dumbAnt [value]. Default is FALSE
#' @param plateau [value]. Default is 100
#' @param returnStats [value]. Default is FALSE
#' @return [value]
#' @details [fill in details here]
#' @examples none
#' @export
acOpt <- function(x, n, targetFunc, pherFunc, xAt0 = FALSE, evapRate = 0.05, nAnts = 500, maxIter = 300, countThresh = 1000, showResult = FALSE, pherPower = 1, dumbAnt = FALSE, plateau = 100, returnStats = FALSE){
	evap <- function(oldP, newP, evapRate) {(1 - evapRate) * oldP + newP }

	if(is.matrix(x)) {if(ncol(x) == 1) x <- c(x) else stop("I cant handle multiple traits yet! dim = ", dim(x))}
	
	if(xAt0){
		minx <- min(x)
		x <- x - minx
	}
	N <- length(x)
	pheromone <- initFunc(x)

	bestAnt <- NULL
	noP = rep(0, N)
	bestestPath = 0
	lastBestPath <- bestestPath
	pathCounter = 0 
	iter = 0
	if(showResult) {
		plotBest <- NULL 
		plotMean <- NULL 
	}

	while(pathCounter <= countThresh & iter < maxIter){
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
		if (path[[bestPath]] == lastBestPath) {
			pathCounter <- pathCounter + 1
		} else {
			lastBestPath <- path[[bestPath]]
			pathCounter <- 0
		}
		if(path[[bestPath]] >= bestestPath) {
			bestestPath <- path[[bestPath]] 
			bestAnt <- ant[[bestPath]]
		}
		newPheromone <- pherFunc(Reduce("+", antPheromone), p = pherPower)	
		pheromone <- evap(pheromone, newPheromone, evapRate = evapRate)
	}
	return(bestAnt)
}
