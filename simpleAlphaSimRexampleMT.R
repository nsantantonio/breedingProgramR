library(breedingProgramR)
library(AlphaSimR)

####################################################
# make population and set its genetic architecture #
####################################################

# some population parameters to make founder population
nChrom <- 2 # number of chromosomes
nLoci <- 100 # number of loci per chromosome
nM <- 10 # number of markers
nQTL <- 10 # number of QTL per chromosome
nFounder <- 10 # size of founder population

Vg <- c(2, 1) # genetic variance 
founderh2 <- c(0.5, 0.3) # heritability
Ve <- (1 - founderh2) / founderh2 * Vg # error variance

# set economic weights
econWt = c(0.4, 0.6)
# define functions to get variance estimates from mixed model fit
estVarG <- function(fit){fit@Vu}
estVarP <- function(fit){fit@Vu + fit@Ve}

# create simple founder pop, either simple or from loaded MaCS sim
founderPop <- quickHaplo(nFounder, nChrom, nLoci, genLen = 1, ploidy = 2L, inbred = FALSE)

# Setting Simulation Parameters
SP <- SimParam$new(founderPop)

# add multi-trait (i.e. add/sample qtl)
rho <- matrix(c(1, 0.5, 0.5, 1), 2, 2)
SP$addTraitA(nQtlPerChr=nQTL, mean = c(0, 0), var = Vg, corA = rho)

# add error variance as function of heritability
SP$setVarE(h2 = founderh2)

# add genotype platform (i.e add/sample markers)
SP$addSnpChip(nM)

# run selection for nYears
nYr = 10
intensity = 0.2

# Set phenotypes for founders

founders <- newPop(founderPop) # make pop with phenotypes and qtl

# make larger population from random crosses of founders
nFam <- 10 # number of families
nIndPerFam <- 20 # number of individuals per family
pop0 <- selectCross(founders, nInd = nInd(founders), use = "rand", simParam = SP, nCrosses = nFam, nProgeny = nIndPerFam) 
# pheno(pop0)

# calculate number selected every year
nSelected <- intensity * nFam * nIndPerFam

# This is just to get the indexes for half of each family 
GS <- TRUE
if(GS){
	index <- NULL
	start <- seq(1, nFam * nIndPerFam, nIndPerFam)
	for(i in 1:length(start)) index <- c(index, (start[i]):(start[i] + nIndPerFam / 2 - 1))
} else {
	index <- 1:(nFam * nIndPerFam)
}

# pop0 <- setPheno(pop0, varE = Ve) # not necessary, as they already have phenotypes, but use if you want to change the number of replications, Ve etc

# 
selectCriterion <- "ebv" # this could aslo be "pheno", "bv"
updateModel <- TRUE
popList <- list()
modelList <- list()
predAcc <- list()
gvList <- list()
# initialize first generation
pop <- pop0


# train model using only half of population (half of members from each family)
model <- RRBLUP(pop[index], traits = 1:length(econWt), use = "pheno", snpChip = 1, simParam = SP)
b = smithHazel(econWt, estVarG(model), estVarP(model))

for(i in 1:nYr){

	pop <- setPheno(pop, varE = Ve, reps = 1)

	# train model using only half of population (half of members from each family)
	if(updateModel & i > 1) {
		model <- RRBLUP(pop[index], traits = 1:length(econWt), use = "pheno", snpChip = 1, simParam = SP)
		b = smithHazel(econWt, estVarG(model), estVarP(model))
	}

	# set ebv of all individuals based on model trained with half of individuals
	pop <- setEBV(pop, model, simParam = SP)

	# get true mean genetic value
	gvList[[i]] <- mean(bv(pop))

	# get true prediction accuracy
	predAcc[[i]] <- cor(pop@ebv, bv(pop))

	# save current geration in popList
	popList[[i]] <- pop

	# save current model in modelList
	modelList[[i]] <- model

	# select and make new population (i.e. next generation)
	pop <- selectCross(pop, nInd = nSelected, use = "ebv", trait = selIndex, 
					   b = b, simParam = SP, nCrosses = nFam, nProgeny = nIndPerFam) 
	# check out makeCross(), allows the user to define mate pairs
}

# extract breeding values, and plot by year
popGVs <- lapply(popList, gv)
meanGVs <- sapply(popGVs, colMeans)
plot(1:nYr, meanGVs[1, ], pch = 16, ylim = c(0, max(meanGVs)))
points(1:nYr, meanGVs[2, ], pch = 17, col = "firebrick")

# extract Vg and plot by year
popVar <- sapply(popList, function(x) diag(varA(x)))
lines(1:nYr, popVar[1, ])
lines(1:nYr, popVar[2, ], col = "firebrick")

# add line for prediction accuracy
popPA <- sapply(predAcc, diag)
lines(1:nYr, popPA[1, ], lty = 2)
lines(1:nYr, popPA[2, ], lty = 2, col = "firebrick")

