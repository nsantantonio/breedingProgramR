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

Vg <- 1 # genetic variance 
Vgxe <- c(1, 2) # genetic variance 
founderh2 <- c(0.5, 0.3) # heritability
Ve <- (1 - founderh2) / founderh2 * Vg # error variance

# create simple founder pop, either simple or from loaded MaCS sim
founderPop <- quickHaplo(nFounder, nChrom, nLoci, genLen = 1, ploidy = 2L, inbred = FALSE)

# Setting Simulation Parameters
SP <- SimParam$new(founderPop)

# add trait (i.e. add/sample qtl)
rho <- matrix(c(1, 0.5, 0.5, 1), 2, 2)
SP$addTraitAG(nQtlPerChr=nQTL, mean = c(0, 0), var = c(Vg, Vg), varGxE = Vgxe, corGxE = rho)

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
model <- RRBLUP(pop[index], traits = 1, use = "pheno", snpChip = 1, simParam = SP)

for(i in 1:nYr){

	pop <- setPheno(pop, varE = Ve, reps = 1)

	# train model using only half of population (half of members from each family)
	if(updateModel & i > 1) model <- RRBLUP(pop[index], traits = 1, use = "pheno", snpChip = 1, simParam = SP)

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
	pop <- selectCross(pop, nInd = nSelected, use = selectCriterion, simParam = SP, nCrosses = nFam, nProgeny = nIndPerFam) 
	# check out makeCross(), allows the user to define mate pairs
}

# extract breeding values, and plot by year
popGVs <- lapply(popList, gv)
meanGVs <- sapply(popGVs, mean)
plot(1:nYr, meanGVs, pch = 16, ylim = c(0, max(meanGVs)))

# extract Vg and plot by year
popVar <- sapply(popList, varA)
lines(1:nYr, popVar)

# add line for prediction accuracy
lines(1:nYr, unlist(predAcc), lty = 2)


### marker effects change through time
getMeff <- function(x) x@markerEff
mEffs <- sapply(modelList, getMeff)

plot(1:nYr, mEffs[1,], type = "l", ylim = range(mEffs))
for(i in 2:nrow(mEffs)) lines(1:nYr, mEffs[i,])