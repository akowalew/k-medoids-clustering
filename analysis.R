# 	analysis.R
# Implements statistical functions which may be run after clustering. 
#
# Author: akowalew

source("k-medoids-clustering.R")
source("data-loading.R")

# Calculates variance between elements in cluster and its medoid
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

# Returns a vector of variances for every cluster
calcClustersVariances <- function(data, clusters, medoids) {
	clustersVariances <- sapply(1:length(clusters), function(k) {
		cluster <- clusters[[k]]
		medoid <- medoids[k]
		clusterVariance <- calcClusterVariance(data, cluster, medoid)	
	})
}

# Calculates dispersion factor for final clustering result
calcClusteringDispersion <- function(data, clusters, medoids) {
	clustersVariance <- calcClustersVariances(data, clusters, medoids)
	variance <- calcVariance(data)
	dispersion <- (sum(clustersVariance) / variance)
}

# Calculates variance for whole dataset
calcVariance <- function(data) {
	xmean <- colMeans(data)
	squares <- apply(data, 1, function(x) (x-xmean)*(x-xmean))
	variance <- sum(squares) / nrow(data)
}

# Runs kMedoidsClustering algorithm and adds some statistics to result
doClustering <- function(data, nclusters, metric) {
	result <- kMedoidsClustering(data, nclusters, metric)

	clusters <- result$clusters
	medoids <- result$medoids

	# perform statistical analysis
	clustersVariances <- calcClustersVariances(data, clusters, medoids)
	clustersStdDevs <- sqrt(clustersVariances)
	variance <- calcVariance(data)
	dispersion <- calcClusteringDispersion(data, clusters, medoids)
	stdDev <- sqrt(variance)
	xmean <- colMeans(data)

	result$clustersVariances <- clustersVariances
	result$clustersStdDevs <- clustersStdDevs
	result$variance <- variance
	result$dispersion <- dispersion
	result$stdDev <- stdDev
	result$xmean <- xmean

	return(result)
}