# 	k-medoids-clustering.R
# Contains k-medoids-clustering algorithm and helper functions
#
# Purpose of Partitioning Around Medoids is to clusterize input
# dataset into clusters. Each cluster is represented by its medoid.
# Algorithm is searching best combination of clustering in order
# to minimize cost function, which is simply a sum of distances
# between elements inside clusters and medoids.
#
# Usage:
# 	data <- replicate(2, rnorm(20))
#	nclusters <- 3
# 	metric <- "euclidean"
# 	kMedoidsClustering(data, nclusters, metric)
# 
# Author: akowalew

# Calculates cost for particular cluster, as a sum of distances
# between its elements and medoid.
calcCost <- function(data, cluster, medoid, metric) {
	if(length(cluster) == 0)
		return(NA)

	v <- data[medoid, ]
	distances <- sapply(cluster, function(elem) {
		x <- data[elem, ]
		distance <- dist(rbind(x, v), method=metric)
	})

	sum(distances)
}

# Calculates cost for every cluster.
calcCosts <- function(data, clusters, medoids, metric) {
	costs <- sapply(1:length(clusters), function(k) {
		cluster <- clusters[[k]]
		medoid <- medoids[k]
		cost <- calcCost(data, cluster, medoid, metric)	
	})
}

# Calculates distances between given point 'x' and medoids.
calcDistances <- function(data, x, medoids, metric) {
	sapply(medoids, function(medoid) {
		v <- data[medoid, ]	
		dist(rbind(x, v), method=metric)
	})
}

# Returns number of medoid, to which given element has shortest distance
getNearestMedoidNumber <- function(data, elem, medoids, metric) {
	x <- data[elem, ]
	distances <- calcDistances(data, x, medoids, metric)
	k <- which.min(distances)
}

# Creates list, where each node contains vector of data elements
# which are assigned to particular clusters
makeClusters <- function(data, elements, medoids, metric) {
	clusters <- vector("list", length(medoids))
	for(elem in elements) {
		k <- getNearestMedoidNumber(data, elem, medoids, metric)
		clusters[[k]] <- c(clusters[[k]], elem)
	}
	return(clusters)
}

# Checks isolation of given cluster.
# Cluster is isolated, when after replacing medoid with its center
# 	we are getting lower cost value.
isClusterIsolated <- function(data, cluster, medoid, metric) {
	if(length(cluster) == 0)
		return(TRUE)

	currentCost <- calcCost(data, cluster, medoid, metric)

	# get mean point of whole group
	group <- data[c(cluster, medoid), ]
	xmean <- colMeans(group)

	# replace medoid with mean and compare its cost to current one
	newCost <- sum(sapply(cluster, function(elem) {
		x <- data[elem, ]
		distance <- dist(rbind(x, xmean), method=metric)
	}))

	isolated <- (newCost < currentCost)
}

# Returns vector representing an isolation of all clusters
calcIsolations <- function(data, clusters, medoids, metric) {
	isolations <- sapply(1:length(clusters), function(k) {
		cluster <- clusters[[k]]
		medoid <- medoids[k]

		isolated <- isClusterIsolated(data, cluster, medoid, metric)
	})
}

# Generates vector: membership of data with particular clusters
makeClustering <- function(data, clusters, medoids) {
	clustering <- vector(mode="integer", nrow(data))
	for(i in 1:length(clusters)) {
		cluster <- clusters[[i]]

		# if there is any element inside a cluster, assign it to data
		if(length(cluster) > 0)
			clustering[cluster] <- rep(i, length(cluster))

		# medoid is always assigned to cluster
		clustering[medoids[i]] <- i
	}
	return(clustering)
} 

# Implements PAM algorithm. It takes `data` matrix and is clustering them
# into `nclusters` partitions. All costs and distances are evaluated using
# specified `metric`.
# To improve speed, clustering results are represented by indices to exact 
# rows inside `data` matrix. 
kMedoidsClustering <- function(data, nclusters=5, metric="euclidean") {
	stopifnot(is.matrix(data))
	stopifnot(nclusters %% nclusters == 0) # nclusters must be int
	stopifnot(is.character(metric))
	stopifnot(nrow(data) > nclusters)

	# get initial random sample medoids
	medoids <- sample(1:nrow(data), nclusters)
	elements <- (1:nrow(data))[-medoids]

	repeat {
		clusters <- makeClusters(data, elements, medoids, metric)

		# find best medoid inside each cluster
		hasChanged <- FALSE
		for(k in 1:nclusters) {
			medoid <- medoids[k]
			cluster <- clusters[[k]]

			# if cluster is empty, skip iteration
			if(length(cluster) == 0)
				next

			bestCost <- calcCost(data, cluster, medoid, metric)
			bestIdx <- -1
			
			# calculate costs for each replacement of cluster element and medoid
			for(i in 1:length(cluster)) {
				# swap
				newMedoid <- cluster[i]
				cluster[i] <- medoid
				medoid <- newMedoid

				# test, whether new combination gives lower cost
				cost <- calcCost(data, cluster, medoid, metric)
				if(cost < bestCost) {
					bestCost <- cost
					bestIdx <- i
				}

				# un-swap
				oldMedoid <- cluster[i]
				cluster[i] <- medoid
				medoid <- oldMedoid
			}

			# if better option was found, replace permamently
			if(bestIdx != -1) {
				hasChanged <- TRUE
				newMedoid <- cluster[bestIdx]
				clusters[[k]][bestIdx] <- medoid
				medoids[[k]] <- newMedoid
			}
		}

		if(!hasChanged)
			break

		elements <- unlist(clusters)
	}

	clustering <- makeClustering(data, clusters, medoids)
	costs <- calcCosts(data, clusters, medoids, metric)
	isolations <- calcIsolations(data, clusters, medoids, metric)

	return(list(
		"medoids"=medoids,
		"clusters"=clusters,
		"clustering"=clustering,
		"costs"=costs,
		"isolations"=isolations))
}