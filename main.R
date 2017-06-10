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

showResults <- function(gaPopulation, result, D) {
	plot(gaPopulation[,2], gaPopulation[,3])

	clusters <- result$clusters
	medoids <- result$medoids

	colors <- rainbow(length(clusters))
	# plot clusters members
	for(idx in 1:length(clusters)) {
		cluster <- clusters[[idx]]
		if(length(cluster) == 0) {
			next
		}

		clusterPoints <- matrix(unlist(cluster), ncol=D, byrow=TRUE)
		points(clusterPoints[,2], clusterPoints[,3], "col"=colors[idx])
	}

	# plot clusters medoids
	for(idx in 1:length(clusters)) {
		medoid <- medoids[[idx]]
		points(medoid[2], medoid[3], "col"=colors[idx], "pch"=23, lwd=4)
	}
}

# calculate sum of costs for each element in cluster
calculateCostsMean <- function(clusters, medoids, distMethod) {
	totalCost <- 0
	for(i in 1:length(clusters)) {
		cluster <- clusters[[i]]
		medoid <- medoids[[i]]
		cost <- calculateCost(medoid, cluster, distMethod)
		totalCost <- totalCost + cost
	}
	costsMean <- totalCost / length(clusters)
	return(costsMean)
}

calculatePopValues <- function(gaPopulation) {
	popSize <- dim(gaPopulation)[1];

	values <- matrix(data=NA, ncol=1, nrow=popSize)
	for(j in 1:popSize) {
		elem <- gaPopulation[j,]
		elemValue <- targetFunction(elem)
		values[j,] <- elemValue
	}

	return(values)
}

# calculate mean of output values for whole population
calculatePopMean <- function(gaPopulation) {
	popSize <- dim(gaPopulation)[1]
	values <- calculatePopValues(gaPopulation)
	return(mean(values))
}

# calculate variance of output values for whole population
calculatePopVariance <- function(gaPopulation) {
	popDim <- dim(gaPopulation)
	popSize <- popDim[1]; D <- popDim[2]

	popMean <- calculatePopMean(gaPopulation)
	variance <- 0
	for(j in 1:popSize) {
		elem <- gaPopulation[j,]
		elemValue <- targetFunction(elem)
		variance <- variance + (elemValue - popMean)^2
	}	
	variance <- variance / popSize

	return(variance)
}

# calculate variance of output values for specified cluster according to its medoid
calculateClusterVariance <- function(cluster, medoid) {
	if(length(cluster) == 0) {
		return(0)
	}

	medoidValue <- targetFunction(medoid)
	variance <- 0
	for(elem in cluster) {
		elemValue <- targetFunction(elem)
		variance <- variance + (elemValue - medoidValue)^2
	}
	variance <- variance / length(cluster)

	return(variance)
}

# calculate total dispersion between clusters
calculateClusteringDispersion <- function(gaPopulation, clusters, medoids) {
	clustersVariancesSum <- 0

	for(i in 1:length(clusters)) {
		cluster <- clusters[[i]]
		medoid <- medoids[[i]]
		clusterVariance <- calculateClusterVariance(cluster, medoid)
		clustersVariancesSum <- clustersVariancesSum + clusterVariance
	}

	popVariance <- calculatePopVariance(gaPopulation)
	dispersion <- clustersVariancesSum / popVariance
	return(dispersion)
}

test1 <- function(gaPopulation, nclusters) {
	# TEST#1
	# test every method's dispersions and costs change, modifying number of clusters
	distMethods <- c("euclidean", "manhattan")
	methodsDispersions <- vector("list", length(distMethods))
	methodsCostsMeans <- vector("list", length(distMethods))

	nclustersRange <- 1:30
	for(i in 1:length(distMethods)) {
		distMethod <- distMethods[i]

		dispersions <- vector("double", length(nclustersRange))
		costsMeans <- vector("double", length(nclustersRange))

		for(nclusters in nclustersRange) {
			result <- kMedoidsClustering(gaPopulation, distMethod, nclusters)
			clusters <- result$clusters
			medoids <- result$medoids
			
			dispersion <- calculateClusteringDispersion(gaPopulation, clusters, medoids)
			costsMean <- calculateCostsMean(clusters, medoids, distMethod)

			dispersions[nclusters] <- dispersion
			costsMeans[nclusters] <- costsMean
		}

		methodsDispersions[[i]] <- dispersions
		methodsCostsMeans[[i]] <- costsMeans
	}

	return(list(methodsDispersions=methodsDispersions,
		methodsCostsMeans=methodsCostsMeans,
		distMethods=distMethods))
}

showTest1 <- function(D, left, right, customSeed) {
	popSizes <- c(50, 100)
	test1Results <- vector("list", length(popSizes))

	for(i in 1:length(popSizes)) {
		popSize <- popSizes[i]

		fileName <- paste("test1", functionNumber, D, popSize, sep="_")
		if(file.exists(fileName)) {
			print(paste("File:", fileName, "found!"))
			load(fileName)
		} else {
			print(paste("File:", fileName, "not found! Generating a new one..."))
			popSize <- 50
			set.seed(customSeed)
			GA <- runGA(D, popSize, left, right, fitnessFunction)
			gaPopulation <- GA@population
			test1Result <- test1(gaPopulation, nclusters)
			save(test1Result, file=fileName)
		}

		test1Results[[i]] <- test1Result
	}
}

# startTesting <- function(D) {
# 	# D <- 2
# 	left <- -100
# 	right <- 100
# 	customSeed <- 1234
# 	# showTest1(D, left, right, customSeed)

# 	popSize <- 100
# 	# dimensions <- c(10, 30, 50)
# 	# for(D in dimensions) {
# 		GA <- runGA(D, popSize, left, right, fitnessFunction)
# 	# }
# 	return(GA)
# }
