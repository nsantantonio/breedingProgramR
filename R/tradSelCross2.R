#' tradSelCross2 function
#'
#' function to (do something)
#'
#' @param pop [value]
#' @param elitepop [value]
#' @param nCrosses [value]
#' @param nFam [value]
#' @param familySize [value]
#' @param families [value]
#' @param use [value]
#' @param best [value]. Default is TRUE
#' @param nProgeny [value]. Default is 1
#' @param intWithin [value]. Default is 0.2
#' @param intAcross [value]. Default is 1
#' @param equalWeight [value]. Default is FALSE
#' @param useFamPrior [value]. Default is FALSE
#' @param verbose [value]. Default is FALSE
#' @return [value]
#' @details [fill in details here]
#' @examples none
#' @export
tradSelCross2 <- function(pop, elitepop, nCrosses, nFam, familySize, families, use, simParam = NULL, best = TRUE, nProgeny = 1, intWithin = 0.2, intAcross = 1, equalWeight = FALSE, useFamPrior = FALSE, verbose = FALSE){
	if(is.character(use)) use <- match.fun(use)

	if(any(table(elitepop@id) > 1)) msg(2, "WARNING: some elites repeated!")
	if(!all(pop@id %in% unlist(families))) stop("wrong family information... fix me!")
	famL <- lapply(families, function(x) x[x %in% pop@id])
	repFam <- sapply(famL, length) > 0
	famNum <- sum(repFam)
		
	if(famNum == 1) msg(2, "NOTE: Only one family represented!")
	# if(nInd(elitepop) > nFam) elitepop <- elitepop[sample(nInd(elitepop), nFam)]

	nFamSel <- ceiling(nFam * intAcross)
	if(famNum < nFamSel) msg(2, "NOTE: insufficient families represented to use specified intAcross. Using individuals from", famNum ,"represented families")
	msg(2, "Number of families available for crossing", min(famNum, nFamSel))

	famMeans <- sapply(famL[repFam], function(x) mean(use(pop[x])))
	famSel <- names(sort(famMeans, decreasing = TRUE))[1:min(famNum, nFamSel)]

	candFamL <- famL[famSel]
	inFamNum <- round(familySize * intWithin)

	parentL <- list()
	for(i in famSel) {
		parentL[[i]] <- selectInd(pop[candFamL[[i]]], nInd = min(nInd(pop[candFamL[[i]]]), inFamNum), simParam = simParam)
		parentL[[i]] <- parentL[[i]][order(use(parentL[[i]]), decreasing = TRUE)]
	}

	candL <- lapply(parentL, function(x) x@id)
	selFamMeans <- sapply(parentL, function(x) mean(use(x)))
	
	parents <- mergePopsRec(parentL)
	parents <- parents[order(use(parents), decreasing = TRUE)]
	parNames <- parents@id
	elNames <- elitepop@id
	
	needToAvoidInbreeding <- any(elNames %in% parNames)
	if(needToAvoidInbreeding) {
		msg(2, "Selection candidates already in elite pop! Avoiding within family crosses...")
	}

	if(nInd(parents) < nCrosses) msg(2, "NOTE: insufficient individuals to make desired number of unique crosses! Some parents will be reused...")

	if(!best){
		fm <- if (useFamPrior) famMeans[famSel] else selFamMeans 
		inFamWeights <- lapply(parentL, function(x) use(x)^2 / sum(use(x)^2))
		w <- if (equalWeight) rep(1/nFamSel, nFamSel) else fm^2 / sum(fm^2)
		whichFams <- sample(famSel, nCrosses, replace = TRUE, prob = w)
	} else {
		whichFams <- sapply(parNames, function(x) names(candL)[sapply(candL, function(xx) x %in% xx)])
	}

	selection <- list()
	for(i in 1:nCrosses){
		p1 <- if(best) parNames[i] else sample(candL[[whichFams[i]]], 1, prob = c(inFamWeights[[whichFams[i]]]))
		notRelated <- !elNames %in% famL[[whichFams[i]]]
		if(sum(notRelated) == 0) {
			nR <- !elitepop@id %in% famL[[whichFams[i]]]
			if(sum(nR) == 0){
				msg(1, "No unrelated selection candidates! Using lines from another family in same generation... ")
				otherFam <- sample(famSel[!famSel %in% whichFams[i]], 1)
				if(best) {
					p2 <- candL[[otherFam]][1]
				} else {
					p2 <- sample(candL[[otherFam]], 1, prob = c(inFamWeights[[otherFam]]))
				}
			} else {
				p2 <- sample(elitepop@id[nR], 1)
			}
		} else {		
			p2cand <- elNames[notRelated]
			elPar <- which(elNames == sample(p2cand, 1))
			p2 <- elNames[elPar]
			elNames <- elNames[-elPar]
		}
		selection[[i]] <- c(p1, p2)
		if(length(p1) > 1 | length(p2) > 1) break
	}

	selection <- do.call(rbind, selection)
	crossFamilies <- lapply(candL, function(x) x[x %in% selection[, 1]])
	if(verbose) msg(2, "Number of lines selected from each family for crosses:", sapply(crossFamilies, length))

	if(nProgeny > 1) selection <- selection[rep(1:nrow(selection), each = nProgeny), ] 
	newpop <- makeCross(mergePopsRec(list(pop, elitepop)), crossPlan = selection, simParam = simParam) 

	if(verbose) msg(2, "Selection mean:", round(mean(gv(newpop)), 6), "from pop mean:", round(mean(gv(pop)), 6))
	newpop
}
