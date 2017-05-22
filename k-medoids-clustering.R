# sum of all distances to medoid 
# of elements inside cluster
calculateCost <- function(medoid, cluster, distMethod) {
	totalDistance <- 0
	for(elem in cluster) {
		distance <- dist(rbind(elem, medoid), method=distMethod)
		totalDistance <- totalDistance + distance
	}
	return(totalDistance)
}

kMedoidsClustering <- function(gaPopulation, distMethod, clustersSize) {
	# initialize clusters
	clusters <- vector("list", clustersSize)
	for(i in 1:clustersSize) {
		cluster <- list()
		clusters[[i]] <- cluster
	}

	# get initial random medoids
	# and remove them from gaPopulation
	medoids <- vector("list", clustersSize)
	set.seed(customSeed)
	idxs <- sample(1:popSize, clustersSize)
	for(i in 1:clustersSize) {
		idx <- idxs[i]
		medoid <- gaPopulation[idx,]
		medoids[[i]] <- medoid
	}
	gaPopulation <- gaPopulation[-idxs,]

	# convert gaPopulation (represented by a matrix)
	# into population variable (represented by list)
	restPopSize <- popSize - clustersSize
	population <- vector("list", restPopSize)
	for(idx in 1:restPopSize) {
		elem <- gaPopulation[idx,]
		population[[idx]] <- elem
	}

	# main algorithm loop
	while(1) {
		# assign elements of population into clusters
		for(elem in population) {
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

		# for each cluster find better medoid
		wasChanged <- FALSE
		for(i in 1:clustersSize) {
			cluster <- clusters[[i]]

			# if cluster is empty, nothing to find
			if(length(cluster) == 0) {
				next 
			}

			medoid <- medoids[[i]]
			currentCost <- calculateCost(medoid, cluster, distMethod)

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
				cost <- calculateCost(medoid, cluster, distMethod)
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

		# check, if elements can be assigned to other medoids
		# but if there was no changes inside clusters, end algorithm
		if(!wasChanged)
			break

		# unpack all elements from clusters into population variable
		population <- vector("list", restPopSize)
		idx <- 1
		for(k in 1:clustersSize) {
			cluster <- clusters[[k]]
			if(length(cluster) != 0) {
				population[idx:(length(cluster)+idx-1)] <- cluster
				idx <- idx + length(cluster)
				clusters[[k]] <- list()
			}
		}
	}

	return(list("clusters" = clusters, "medoids" = medoids))
}