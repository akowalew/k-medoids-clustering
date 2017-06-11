# 	data-loading.R
# Contains functions to generate and save/load GeneticAlgorithm 
# populations. These populations are generated using `GA` package 
# by minimizing fitnessFunctions from `cec2013` package.
#
# Author: akowalew

library("GA")
library("cec2013")

left <- -100; right <- 100; popSize <- 100

# Utility function, used by `createPopulation` and `getPopulation`
getFileName <- function(D, fnumber) {
	fileName <- paste("cec2013", fnumber, D, popSize, left, right, sep="_")
}

createPopulation <- function(D, fnumber) {
	fileName <- getFileName(D, fnumber)
	print(paste("File with name:", fileName, "doesn't exists. Generating population..."))

	fitness <- function(x) -cec2013(fnumber, x)

	GA <- ga(type="real-valued",
		fitness=fitness,
		min=rep(left, D),
		max=rep(right, D),
		maxiter=999999,
		maxFitness=Inf,
		popSize=popSize,
		parallel=FALSE,
		optim=TRUE)

	save(GA, file=fileName);
	print("Saved generated population!")

	return(GA@population)
}

getPopulation <- function(D, fnumber) {
	fileName <- getFileName(D, fnumber)
	if(file.exists(fileName)) {
		print(paste("File with name:", fileName, "exists. Loading it..."))
		load(fileName)
		return(GA@population)
	} 
	else {
		# return(createPopulation(D, fnumber))
		return(NULL)
	}
}