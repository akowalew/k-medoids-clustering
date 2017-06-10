assignDataToMedoids <- function(data, medoids, metric) {
	apply(data, 1, function(x) {
		distances <- apply(medoids, 1, 
			function(medoid) dist(rbind(x, medoid), method=metric))
		which.min(distances)	
	})
}

kMedoidsClustering <- function(data, nclusters, metric) {
	clustering <- rep(0, nrow(data))
	medoidsIdxs <- sample(1:nrow(data), nclusters)

	repeat {
		medoids <- data[medoidsIdxs, ]

		restIdxs <- (1:nrow(data))[-medoidsIdxs]
		restData <- data[restIdxs, ]
		restClustering <- assignDataToMedoids(restData, medoids, metric)

		clustering[restIdxs] <- restClustering
		clustering[medoidsIdxs] <- 1:nclusters # always the same

		# for each cluster find best medoid
		newMedoidsIdxs <- sapply(1:nclusters, function(cluster) {
			# get elements assigned to cluster
			clusteredIdxs <- which(clustering == cluster)
			nclustered <- length(clusteredIdxs)

			# if cluster is empty - only medoid - return index to it
			if(nclustered == 1)
				return(clusteredIdxs[1])

			# calculate costs for each cluster
			costs <- sapply(1:nclustered, function(idx) {
				medoid <- data[clusteredIdxs[idx], ]	
				clusterData <- data[clusteredIdxs[-idx], ]

				if(!length(dim(clusterData))) { # if there is only one element
					distance <- dist(rbind(clusterData, medoid), method=metric)
					distances <- c(distance)
				} else {
					distances <- apply(clusterData, 1, 
						function(x) dist(rbind(x, medoid), method=metric))
				}

				cost <- sum(distances)
			})

			clusteredIdxs[which.min(costs)]
		})

		if(all(medoidsIdxs == newMedoidsIdxs))
			break;
			
		medoidsIdxs <- newMedoidsIdxs
	}

	medoids <- data[medoidsIdxs, ]

	return(list("medoids"=medoids, "clustering"=clustering))
}