#' tradSelCross function
#'
#' function to (do something)
#'
#' @param pop [value]
#' @param nCrosses [value]
#' @param nFam [value]
#' @param famSize [value]
#' @param families [value]
#' @param use [value]
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
tradSelCross <- function(pop, nCrosses, nFam, famSize, families, use, nProgeny = 1, intWithin = 0.2, intAcross = 1, equalWeight = FALSE, useFamPrior = FALSE, verbose = FALSE){
# RCRS[[gen(j)]] <- tradSel(pop = selPop, nInd = min(selectRCRSi, nInd(selPop)), use = useIn,  trait = 1, simParam = simParam, nCrosses = nNuclear, nProgeny = nProgenyPerCrossIn)
	# lpop <- function(x, l, whichElem = NULL) {
	# 	whichElem <- if(is.null(whichElem)) which(sapply(l, function(elem) x %in% elem))
	# 	for(i in whichElem) l[[i]] <- l[[i]][!l[[i]] %in% x]
	# 	l
	# }
	if(is.character(use)) use <- match.fun(use)

	if(!all(pop@id %in% unlist(families))) stop("wrong family information... fix me!")
	famL <- lapply(families, function(x) x[x %in% pop@id])
	# famL[10] <- NULL
	repFam <- sapply(famL, length) > 0
	famNum <- sum(repFam)
	if(famNum == 1){
		msg(2, "NOTE: Only one family represented! Continuing by mating within familiy...")
		newpop <- selectCross(pop, nInd = round(nInd(pop) * intWithin), nCrosses = nCrosses)
	} else {
		nFamSel <- ceiling(nFam * intAcross)
		if(famNum < nFamSel) msg(2, "NOTE: insufficient families represented to use specified intAcross. Using individuals from all", famNum ,"represented families")
		msg(2, "Number of families used for crossing", min(famNum, nFamSel))

		famMeans <- sapply(famL[repFam], function(x) mean(use(pop[x])))
		famSel <- names(sort(famMeans, decreasing = TRUE))[1:min(famNum, nFamSel)]

		candFamL <- famL[famSel]
		inFamNum <- famSize * intWithin

		parentL <- list()
		for(i in famSel) {
			parentL[[i]] <- selectInd(pop[candFamL[[i]]], nInd = min(nInd(pop[candFamL[[i]]]), inFamNum))
			parentL[[i]] <- parentL[[i]][order(use(parentL[[i]]), decreasing = TRUE)]
		}

		selFamMeans <- sapply(parentL, function(x) mean(use(x)))
		
		parents <- mergePopsRec(parentL)
		parents <- parents[order(use(parents), decreasing = TRUE)]
		parNames <- parents@id

		if(choose(nInd(parents), 2) < nCrosses) {
			msg(2, "NOTE: insufficient individuals to make desired number of unique crosses! Some mate pairs will be repeated...")
			nTimes <- ceiling(nCrosses / choose(nInd(parents), 2))
			selL <- list(t(combn(parents@id, 2)))
			selection <- do.call(rbind, rep(selL, nTimes - 1))
			selection <- rbind(selection, selL[[1]][sample(1:nrow(selL[[1]])), ])[1:nCrosses, ]
		} else {
			selection <- list()
			if(intWithin == 1) {
				wp <- 1
				for(i in 1:nCrosses) {
					p1 <- parNames[1]
					diffFam <- !sapply(candFamL, function(f) p1 %in% f)
					cand <- parNames[parNames %in% unlist(candFamL[diffFam])]
					if(length(parNames) <= 1) {
						parNames <- parents@id
						wp <- wp + 1
					}
					p2 <- parNames[parNames %in% unlist(candFamL[diffFam])][[wp]]
					selection[[i]] <- c(p1, p2)
					parNames <- parNames[!parNames %in% c(p1, p2)]
				}
			} else {
				nFamPairs <- choose(nFamSel, 2)
				famPairs <- combn(famSel, 2)
				fm <- if (useFamPrior) famMeans[famSel] else selFamMeans 
				inFamWeights <- lapply(parentL, function(x) use(x)^2 / sum(use(x)^2))
				pbar <- combn(fm, 2, FUN = mean)
				w <- if (equalWeight) rep(1/nFamPairs, nFamPairs) else pbar^2 / sum(pbar^2)
	 			withRep <- if(nCrosses > nFamPairs) TRUE else FALSE
				whichPairs <- sample(1:ncol(famPairs), nCrosses, replace = withRep, prob = w)
				for(i in 1:nCrosses){
					pairi <- famPairs[, whichPairs[i]]
					p1 <- sample(parentL[[pairi[1]]]@id, 1, prob = c(inFamWeights[[pairi[1]]]))
					p2 <- sample(parentL[[pairi[2]]]@id, 1, prob = c(inFamWeights[[pairi[2]]]))
					selection[[i]] <- c(p1, p2)
				}
			}
			if(any(sapply(selection, function(x) length(unique(x))) != 2)) stop("something wrong happened!")
			selection <- do.call(rbind, selection)
		}

		if(nProgeny > 1) selection <- selection[rep(1:nrow(selection), each = nProgeny), ] 
		newpop <- makeCross(pop, crossPlan = selection) 
	}
	if(verbose) msg(2, "Selection mean:", round(mean(gv(newpop)), 6), "from pop mean:", round(mean(gv(pop)), 6))
	newpop
}
