calcCost <- function(data, cluster, medoid, metric) {
	v <- data[medoid, ]
	distances <- sapply(cluster, function(elem) {
		x <- data[elem, ]
		distance <- dist(rbind(x, v), method=metric)
	})
	sum(distances)
}

calcCosts <- function(data, clusters, medoids, metric) {
	costs <- sapply(1:length(clusters), function(k) {
		cluster <- clusters[[k]]
		medoid <- medoids[k]
		cost <- calcCost(data, cluster, medoid, metric)	
	})
}

calcDistances <- function(data, x, medoids, metric) {
	sapply(medoids, function(medoid) {
		v <- data[medoid, ]	
		dist(rbind(x, v), method=metric)
	})
}

makeClusters <- function(data, elements, medoids, metric) {
	clusters <- vector("list", length(medoids))
	for(elem in elements) {
		x <- data[elem, ]

		distances <- calcDistances(data, x, medoids, metric)
		k <- which.min(distances)
		clusters[[k]] <- c(clusters[[k]], elem)
	}
	return(clusters)
}

calcVariance <- function(data) {
	xmean <- colMeans(data)
	squares <- apply(data, 1, function(x) (x-xmean)*(x-xmean))
	variance <- sum(squares) / nrow(data)
}

calcClusterVariance <- function(data, cluster, medoid) {
	clusterSize <- length(cluster)
	if(clusterSize == 0)
		return(0)

	v <- data[medoid, ]
	squares <- sapply(cluster, function(elem) {
		x <- data[elem, ]
		square <- (x-v)*(x-v)
	})
	variance <- sum(squares) / clusterSize
}

# calculates variances for each cluster
calcClustersVariances <- function(data, clusters, medoids) {
	clustersVariances <- sapply(1:length(clusters), function(k) {
		cluster <- clusters[[k]]
		medoid <- medoids[k]
		clusterVariance <- calcClusterVariance(data, cluster, medoid)	
	})
}

# Checks for each cluster, whether its mean point
# after replacement with medoid provides lower cost
calcIsolations <- function(data, clusters, medoids, metric) {
	isolations <- sapply(1:length(clusters), function(k) {
		cluster <- clusters[[k]]
		medoid <- medoids[k]
		currentCost <- calcCost(data, cluster, medoid, metric)

		# get mean point of whole group
		group <- data[c(cluster, medoid), ]
		xmean <- colMeans(group)

		# replace medoid with mean and compare its cost to current one
		newCost <- sum(sapply(cluster, function(elem) {
			x <- data[elem, ]
			distance <- dist(rbind(x, xmean), method=metric)
		}))

		isIsolated <- (newCost < currentCost)
	})
}

calcClusteringDispersion <- function(data, clusters, medoids) {
	clustersVariance <- calcClustersVariances(data, clusters, medoids)
	variance <- calcVariance(data)
	dispersion <- (sum(clustersVariance) / variance)
}

# generates vector: membership of data to clusters
makeClustering <- function(data, clusters, medoids) {
	clustering <- vector(mode="integer", nrow(data))
	for(i in 1:nclusters) {
		cluster <- clusters[[i]]
		clusterSize <- length(cluster)
		if(clusterSize)
			clustering[cluster] <- rep(i, clusterSize)
		clustering[medoids[i]] <- i
	}
	return(clustering)
} 

kMedoidsClustering <- function(data, nclusters, metric) {
	# get initial random sample medoids
	medoids <- sample(1:nrow(data), nclusters)
	elements <- (1:nrow(data))[-medoids]

	repeat {
		# assign elements into medoids
		clusters <- makeClusters(data, elements, medoids, metric)

		# find best medoid inside each cluster
		hasChanged <- FALSE
		for(k in 1:nclusters) {
			medoid <- medoids[k]
			cluster <- clusters[[k]]
			clusterSize <- length(cluster)

			# if cluster is empty, skip iteration
			if(clusterSize == 0)
				next

			bestCost <- calcCost(data, cluster, medoid, metric)
			bestIdx <- -1
			
			# calculate costs for each replacement of cluster element and medoid
			for(i in 1:clusterSize) {
				# swap
				newMedoid <- cluster[i]
				cluster[i] <- medoid
				medoid <- newMedoid

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

	isolations <- calcIsolations(data, clusters, medoids, metric)
	clustersVariances <- calcClustersVariances(data, clusters, medoids)
	clustersStdDevs <- sqrt(clustersVariances)
	variance <- calcVariance(data)
	dispersion <- calcClusteringDispersion(data, clusters, medoids)
	stdDev <- sqrt(variance)
	costs <- calcCosts(data, clusters, medoids, metric)
	xmean <- colMeans(data)

	return(list("data"=data,
		"medoids"=medoids,
		"clusters"=clusters,
		"clustering"=clustering,
		"costs"=costs,
		"isolations"=isolations,
		"clustersVariances"=clustersVariances,
		"clustersStdDevs"=clustersStdDevs,
		"variance"=variance,
		"dispersion"=dispersion,
		"stdDev"=stdDev,
		"xmean"=xmean))
}