rm(list=ls())
# graphics.off()

library("GA")
library("cec2013")
library("parallel")
source("k-medoids-clustering.R")

functionNumber <- 7
targetFunction <- function(x) {
	cec2013(functionNumber, x)
}

fitnessFunction <- function(x) {
	-targetFunction(x)
}

runGA <- function(D, popSize, left, right, fitnessFunction) {
	fileName <- paste("cec2013", functionNumber, D, popSize, left, right, sep="_")
	if(file.exists(fileName)) {
		print(paste("File with name:", fileName, "exists. Loading it..."))
		load(fileName)
	} else {
		print(paste("File with name:", fileName, "doesn't exists. Generating population..."))
		GA <- ga(type="real-valued",
			fitness=fitnessFunction,
			min=rep(left, D),
			max=rep(right, D),
			maxiter=999999,
			maxFitness=799.9,
			popSize=popSize,
			# pcrossover=0.8,
			# pmutation=0.1,
			parallel=FALSE,
			optim=TRUE
			#seed=12345
			# elitism=base::max(1, round(popSize*0.05))
			)
		save(GA, file=fileName);
		print("Saved generated population!")
	}
	return(GA)
}
