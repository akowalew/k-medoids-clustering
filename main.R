# clear all
rm(list = ls())
graphics.off()

library("GA")
library("cec2013")

generatePopulation <- function(popSize) {
	functionNumber <- 5
	fun <- function(x) {
		cec2013(functionNumber, x)
	}

	D <- 2
	left <- -100
	right <- 100
	# GA <- ga(type = "real-valued",
	# 	fitness = function(x) -fun(x),
	# 	min = rep(left, D),
	# 	max = rep(right, D),
	# 	maxiter = 500,
	# 	maxFitness = 999.99,
	# 	popSize = 50,
	# 	pcrossover = 0.8,
	# 	pmutation = 0.1,
	# 	elitism = base::max(1, round(popSize*0.05)),
	# 	seed = 12349876)
	# population <- GA@population
	load("population")
	return(population)
}

popSize <- 50
population <- generatePopulation(popSize)

clustersSize <- 5
distMethod <- "euclidean"

# get initial random medoids
# and remove them from population
medoids <- vector("list", clustersSize)
idxs <- sample(1:popSize, clustersSize)
for(i in 1:clustersSize) {
	idx <- idxs[i]
	medoid <- population[idx,]
	medoids[[i]] <- medoid
}
population <- population[-idxs,]

# init clusters
restPopSize <- popSize - clustersSize
clusters <- vector("list", clustersSize)
for(i in 1:clustersSize) {
	cluster <- list()
	clusters[[i]] <- cluster
}

# assign other elements of population
# to clusters represented by medoids
for(idx in 1:(popSize - clustersSize)) {
	elem <- population[idx,]

	# find the nearest medoid
	minDistance <- Inf
	minDistanceIdx <- -1
	for(i in 1:clustersSize) {
		medoid <- medoids[[i]]
		distance <- dist(rbind(elem, medoid), method=distMethod)
		if(distance < minDistance) {
			minDistance <- distance
			minDistanceIdx <- i
		}
	}

	# push elem to nearest cluster
	cluster <- clusters[[minDistanceIdx]]
	clusters[[minDistanceIdx]] <- c(cluster, list(elem))
}

# sum of all distances to medoid 
# of elements inside cluster
calculateCost <- function(medoid, cluster) {
	totalDistance <- 0
	for(elem in cluster) {
		distance <- dist(rbind(elem, medoid), method=distMethod)
		totalDistance <- totalDistance + distance
	}
	return(totalDistance)
}

# for each cluster find better medoid
wasChanged <- FALSE
for(i in 1:clustersSize) {
	cluster <- clusters[[i]]
	medoid <- medoids[[i]]
	currentCost <- calculateCost(medoid, cluster)

	# by replacing each cluster elem with medoid
	# find the lowest cost from all configurations
	minCost <- Inf
	minCostIdx <- -1
	for(j in 1:length(cluster)) {
		# swap medoid with j-th cluster's elem
		newMedoid <- cluster[[j]]
		cluster[[j]] <- medoid
		medoid <- newMedoid

		# check, if the cost was reduced
		cost <- calculateCost(medoid, cluster)
		if(cost < minCost) {
			minCost <- cost
			minCostIdx <- j
		}

		# undo the swap
		oldMedoid <- cluster[[j]]
		cluster[[j]] <- medoid
		medoid <- oldMedoid
	}

	# if best found configuration is better than current
	if(minCost < currentCost) {
		# swap best elem with current medoid
		bestElem <- cluster[[minCostIdx]]
		cluster[[minCostIdx]] <- medoids[[i]]
		medoids[[i]] <- bestElem
		clusters[[i]] <- cluster
		wasChanged <- TRUE
	}
}