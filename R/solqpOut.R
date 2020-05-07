#' solqpOut function
#'
#' function to (do something)
#'
#' @param pop [value]
#' @param GSfit [value]
#' @param use [value]
#' @param nSel [value]
#' @param nProgeny [value]
#' @param nGenOut [value]
#' @param nGenThisYr [value]
#' @param simParam [value]
#' @param limitN [value]. Default is 0
#' @param lambdaOut [value]. Default is NULL
#' @param fthreshOut [value]. Default is NULL
#' @param gainOut [value]. Default is NULL
#' @param truncqpOut [value]. Default is NULL
#' @param varFunc [value]. Default is varA
#' @param verbose [value]. Default is FALSE
#' @return [value]
#' @details [fill in details here]
#' @examples none
#' @export
solqpOut <- function(pop, GSfit, use, nSel, nProgeny, nGenOut, nGenThisYr, simParam = NULL, limitN = 1, lambdaOut = NULL, fthreshOut = NULL, gainOut = NULL, truncqpOut = NULL, varFunc = varA, verbose = FALSE, ...){
	params <- list(lambdaOut = lambdaOut, fthreshOut = fthreshOut, gainOut = gainOut, truncqpOut = truncqpOut)
	for(i in names(params)){
		if (length(params[[i]]) > 1 & !is.list(params[[i]])) stop(paste0(i, " must be length 1 or a list."))
		if(length(params[[i]]) <= 1) params[[i]] <- rep(list(params[[i]]), nGenOut) else if(length(params[[i]]) != nGenOut) stop(paste0(i, " is the wrong length, must be 1 or cyclePerYr - pullCycle"))
	}
	pop <- setEBV(pop, GSfit, simParam = simParam)
	i <- 0
	popParams <- list()
	if(nGenOut == 0){
		pop <- selectInd(pop, nInd = nSel, use = use, simParam = simParam)
	} else {
		N <- if(limitN > 1) rep(nSel, nGenOut) else if(limitN > 0) c(rep(nInd(pop), nGenOut - 1), nSel) else rep(nInd(pop), nGenOut)
		if(verbose) msg(2, "solqpOut Initial Pop Mean:", round(mean(pop@gv), 6), "PopVar", round(varFunc(pop, simParam = simParam), 6))

		while(i < nGenOut){
			i <- i + 1
			popParams[["mu"]][i] <- mean(pop@gv)
			popParams[["sigmasq"]][i] <- varFunc(pop, simParam = simParam)
			sol <- do.call(solqp, getArgs(solqp, pop = pop, GSfit = GSfit, use = use, nCrosses = N[i], simParam = simParam,
						   nProgeny = nProgeny, lambda = params$lambdaOut[[i]], fthresh = params$fthreshOut[[i]], gain = params$gainOut[[i]], truncqp = params$truncqpOut[[i]], ...))
			if(is.list(sol)){
					solClass <- sapply(sol, class)
					popParams[["lambda"]] <- c(popParams[["lambda"]], sol[["lambda"]])
					pop <- sol[[which(solClass == "Pop")]]
			} else {
				pop <- sol
			}
			pop <- setEBV(pop, GSfit, simParam = simParam)

			if(verbose) msg(2, "solqpOut Pop Mean:", round(mean(pop@gv), 6), "PopVar", round(varFunc(pop, simParam = simParam), 6))
			if(verbose) msg(2, "solqpOut pop size:", nInd(pop))
		}
	}
	popParams[["mu"]][i+1] <- mean(pop@gv)
	popParams[["sigmasq"]][i+1] <- varFunc(pop, simParam = simParam)

	if(limitN < 1) {
		pop <- selectInd(pop, nInd = nSel, use = use, simParam = simParam)
		popParams[["mu"]][i+2] <- mean(pop@gv)
		popParams[["sigmasq"]][i+2] <- varFunc(pop, simParam = simParam)
	}
	if(verbose) msg(2, "solqpOut Final Pop Mean:", round(mean(pop@gv), 6), "PopVar", round(varFunc(pop, simParam = simParam), 6))
	return(c(list(pop = pop), popParams))
}
