library(breedingProgramR)
library(AlphaSimR)

####################################################
# make population and set its genetic architecture #
####################################################

# some population parameters to make founder population
nChrom <- 10 # number of chromosomes
nLoci <- 1000 # number of loci per chromosome
nM <- 100 # number of markers
nQTL <- 100 # number of QTL per chromosome
nFounder <- 100 # size of founder population

Vg <- 1 # genetic variance 
founderh2 <- 0.3 # heritability
Ve <- (1 - founderh2) / founderh2 * Vg # error variance

# create simple founder pop, either simple or from loaded MaCS sim
founderPop <- quickHaplo(nFounder, nChrom, nLoci, genLen = 1, ploidy = 2L, inbred = FALSE)

# Setting Simulation Parameters
SP = SimParam$new(founderPop)

# add trait 
SP$addTraitA(nQtlPerChr=nQTL, mean = 0, var = Vg)

# add error variance as function of heritability
SP$setVarE(h2 = founderh2)

# add genotype platform
SP$addSnpChip(nM)

loci <- pullLoci(SP)
Reduce(intersect, loci)

####################################################
# Define User Arguemnts to be passed to simulation #
####################################################

# set user arguments, as a list that is provided to the paramL argument for simSingleTraitInbred()
 # note: there are a lot more options! see ?simSingleTraitInbred 
userArgs <- list(
# family parameters
	nFounder = 100, # number of population founders
	nRCRS = 100, # number of individuals in RCRS
	nFam = 20, # number of families (i.e. crosses) made each year
	famSize = 75, # number of lines per family
# genetic parameters
	Vg = 1, # genetics variance
	h2 = c(0.3, 0.3, 0.3, 0.5), # heritability within each trial
# program parameters
	nYr = 30, # number of years
	cyclePerYr = 3, # number of years
	# selectTrials = c(0.5, 0.2, 0.5, 0.4), # selection intensity within each trial, note you can also provide a positive integer to indicate the number of lines at each stage
	selectTrials = c(0.5, 0.10, 0.4, 1/3),
	trialLocs = c(1, 2, 3, 5), # number of locations
	trialReps = c(1, 2, 3, 3) # number of reps per trial location
)

# fit a rapid cycle recurrent selection program for 30 years
rapidCycleTruncation <- simSingleTraitInbred(founderPop = founderPop, simParam = SP, paramL = userArgs, returnFunc = getPopStats)

# make a plot of the variety means and recurrent population mean and sd
pdf("RT.pdf")
plot(NA, ylim = c(0, max(rapidCycleTruncation$vy)), xlim = c(0, rapidCycleTruncation$paramL$nYr), ylab = "Genetic Value")
plotPop(rapidCycleTruncation, alphaMean = "22")
legend("bottomright", legend = c("Recurrent population", "Varieties"), pch = c(NA, 1), lty = c(1, NA), col = rep("#000000", 2), bty = "n")
dev.off()


# now lets compare that to a traditional program
userArgsTrad <- userArgs # duplicate userArgs
userArgsTrad$traditional <- 3 # select the best lines out of the third year trials to make crosses to elite germplasm
userArgsTrad$useIn <- "pheno" # select using phenotypes only
userArgsTrad$useOut <- "pheno" # select using phenotypes only
userArgsTrad$intAcross <- 0.5 # select the best 50% of families, then the best lines within each familiy
# note: the traditional program keeps a list of the best nFam elite lines, and mates each canidate line to each elite randomly
traditionalTruncation <- simSingleTraitInbred(founderPop = founderPop, simParam = SP, paramL = userArgsTrad, returnFunc = getPopStats)

# make a plot of the variety means and recurrent population mean and sd
pdf("RTvsTR.pdf")
plot(NA, ylim = c(0, max(c(rapidCycleTruncation$vy, traditionalTruncation$vy))), xlim = c(0, rapidCycleTruncation$paramL$nYr), ylab = "Genetic Value")
plotPop(rapidCycleTruncation)
plotPop(traditionalTruncation, popcol = "#B22222")
legend("bottomright", legend = c("Recurrent population", "Varieties"), pch = c(NA, 1), lty = c(1, NA), col = rep("#000000", 2), bty = "n")
legend("topleft", legend = c("Rapid Cycle Truncation", "Traditional"), pch = 1, lty = 1, col = c("#000000", "#B22222"), bty = "n")
dev.off()


# define an arbitrary function for selection in RGSC, this one mates the best to the best and the worst to the worst, should not change mean, but will change variance 
# selFuncIn
varianceGenerator <- function(pop, GSfit, simParam, bidirIntensity = 0.2, use = "ebv"){
	if(is.character(use)) use <- match.fun(use) # make use a function if provided a character
	pop <- setEBV(pop, GSfit, simParam = simParam) # set ebv with provided model fit, not this doesnt necessarily need to be done
	gebv <- c(use(pop))
	n <- nInd(pop)
	int <- bidirIntensity / 2

	hilow <- quantile(gebv, probs = c(int, (1-int)))
	whichLo <- which(gebv <= hilow[[1]])
	whichHi <- which(gebv >= hilow[[2]])
	
	randomPairs <- function(x, y) {matrix(c(sample(x, 2), sample(y, 2)), 2, 2, byrow = TRUE)}
	selection <- randomPairs(whichLo, whichHi)
	while(nrow(selection) < n) selection <- rbind(selection, randomPairs(whichLo,whichHi))
	makeCross(pop, crossPlan = selection, simParam = simParam) 
}

# this is a function to select the best and worst lines for selection out of the recurrent population and advancement through the VDP
#  selFuncOut & selFuncVDP
divergentSelection <- function(pop, nSel, use = "pheno"){
	if(is.character(use)) use <- match.fun(use) # make use a function if provided a character
	b <- c(use(pop))
	n <- nInd(pop)
	int <- nSel / 2 / n
	hilow <- quantile(b, probs = c(int, (1-int)))
	pop[b <= hilow[[1]] | b >= hilow[[2]]]
}

# modify user args
userArgsDS <- userArgs
userArgsDS$selFuncIn <- varianceGenerator # notice we provide 
userArgsDS$selFuncOut <- divergentSelection
userArgsDS$selFuncVDP <- divergentSelection

# run divergent selection program
diverge <- simSingleTraitInbred(founderPop = founderPop, simParam = SP, paramL = userArgsDS, returnFunc = getPopStats)

# plot results
pdf("RTvsTRvsDS.pdf")
plot(NA, ylim = range(c(rapidCycleTruncation$vy, traditionalTruncation$vy, diverge$vy)), xlim = c(0, rapidCycleTruncation$paramL$nYr), ylab = "Genetic Value")
plotPop(rapidCycleTruncation)
plotPop(traditionalTruncation, popcol = "#B22222")
plotPop(diverge, popcol = "#2222B2")
legend("bottomleft", legend = c("Recurrent population", "Varieties"), pch = c(NA, 1), lty = c(1, NA), col = rep("#000000", 2), bty = "n")
legend("topleft", legend = c("Rapid Cycle Truncation", "Traditional", "Divergent Selection"), pch = 1, lty = 1, col = c("#000000", "#B22222", "#2222B2"), bty = "n")
dev.off()



# Lets try something more complicated. Here we use optimal contribution to select within the recurrent population, and branch out every year to push means. This may take some time to run. 
userArgsOCB <- userArgs
userArgsOCB$selFuncIn <- solqp # function to select based on optimal contribution within RCRS
userArgsOCB$selFuncOut <- solqpOut # function to branch and select based on optimal contribution out of RCRS
userArgsOCB$pullCycle <- 0 # cycle to branch material out at minimum is 0 (i.e. at the beginning of each year) up to the number of cycles per year (which doesnt branch anymore).  
userArgsOCB$phenoRCRS <- 0.5

# note that arguments to selFuncIn and selFuncOut are provided directly to simSingleTraitInbred(), and are passed to the right functions within.
optimalContributionWithBranching <- simSingleTraitInbred(founderPop = founderPop, simParam = SP, paramL = userArgsOCB, returnFunc = getPopStats, fthresh = 0.005, fthreshOut = 0.1)

# Lets see how the branching scheme (blue) compares to the rapid cycle (black) and traditional (red)
pdf("RTvsTRvsOCB.pdf")
plot(NA, ylim = c(0, max(c(rapidCycleTruncation$vy, traditionalTruncation$vy, optimalContributionWithBranching$vy))), xlim = c(0, rapidCycleTruncation$paramL$nYr))
plotPop(rapidCycleTruncation)
plotPop(traditionalTruncation, popcol = "#B22222")
plotPop(optimalContributionWithBranching, popcol = "#22B222")
legend("bottomright", legend = c("Recurrent population", "Varieties"), pch = c(NA, 1), lty = c(1, NA), col = rep("#000000", 2), bty = "n")
legend("topleft", legend = c("Rapid Cycle Truncation", "Traditional", "Optimal Contribution Branching"), pch = 1, lty = 1, col = c("#000000", "#B22222", "#22B222"), bty = "n")
dev.off()