test <- function(i, N, n, showResult = FALSE){
	target <- function(x) mean(x) + 2*sd(x)
	pherFunc <- function(x, p) (x / max(x))^p

	ebvs <- sort(rnorm(N), decreasing = TRUE) 
	sel <- sort(ebvs, decreasing = TRUE)[1:n]

	nCombo <- choose(N, n)
	if(nCombo < 1e6) allCombo <- combn(1:N, n, simplify = FALSE) else stop("I cant do it!!!!!")
	allSol <- sapply(allCombo, function(x) target(ebvs[x]))  

	txtdensity(allSol)
	# range(allSol) 
	target(sel)
	optsol <- allCombo[[which.max(allSol)]]
	print(optsol)

	pheromone <- rep(1, N)
	evap <- function(oldP, newP, evapRate) evapRate * oldP + newP 

	evapRate <- 0.8
	nAnts <- 1000
	noP <- rep(0, N)
	shortestPath <- 0
	pathCounter <- 0
	maxIter <- 1000
	countThresh <- 50
	iter <- 0
	while(pathCounter <= countThresh & iter <= maxIter){
		iter <- iter + 1
		ant <- list()
		path <- list()
		antPheromone <- list()
		for(i in 1:nAnts){
			ant[[i]] <- sample(1:N, n, prob = pheromone / sum(pheromone))
			path[[i]] <- target(ebvs[ant[[i]]])
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
		newPheromone <- pherFunc(Reduce("+", antPheromone), p = power)
		pheromone <- evap(pheromone, newPheromone, evapRate)
		
		if(showResult) {
			cols <- rep("gray80", N)
			cols[bestAnt] <- "black"
			barplot(pheromone, col = cols)
			print(iter)
			Sys.sleep(0.1)
		}
	}
	identical(sort(bestAnt), optsol)
}

set.seed(12345)
N <- 30
n <- 6
reps <- 20
power = 2
test(i, N, n, TRUE)

library(doMC)
if(system("hostname", intern = TRUE) == "Bender") {
	setMKLthreads(1)
	registerDoMC(min(20, reps))
} else {
	registerDoMC(min(20, reps))
}

same <- foreach(i = 1:reps, .combine = c) %dopar% test(i, N, n)
sum(same)