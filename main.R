# 	main.R
# Entry point of solution
# Runs scripts from analysis.R
#
# Author: akowalew

rm(list=ls())
graphics.off()

set.seed(1234)

source("analysis.R")
source("data-loading.R")
load("cec2013_5_2_50_-100_100")

# # normal check
# data <- GA@population
# nclusters <- 5
# metric <- "euclidean"
# r <- doClustering(GA@population, nclusters, metric)

# check all cases
fnumbers <- c(5, 6, 14)
dimensions <- c(10, 30, 50)
nnclusters <- c(3, 5)
metrics <- c("euclidean", "manhattan")

results <- lapply(fnumbers, function(fnumber) {
	lapply(dimensions, function(dimension) {
		population <- getPopulation(dimension, fnumber)
		if(is.null(population))
			return(NULL)

		lapply(nnclusters, function(nclusters) {
			lapply(metrics, function(metric) {
				result <- doClustering(population, nclusters, metric)
				result$data <- population
				return(result)
			})	
		})
	})
})
save(results, file='k-medoids-results')
