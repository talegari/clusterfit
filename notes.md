# clusterfit

> an interface for clustering in R

## Motivation

Implementations are spread over multiple packages with different conventions which makes it hard for experimentation for clustering.

## Idea and Design

- Create a sklearn-like interface still keeping while keeping R's style
- Each clustering algo should have one function
- Passing data, number of clusters and obtaining cluster indexes should be standardized
- Unpassed arguments should become the part of the grid automatically
- Do not eager compute anything
- accept data matrix or dist object and make it seamless

## Artist's picture

```
# create a 'kmeans' object with params. Sets up a grid for k.
# arg called implementation is hidden, advisable to not meddle with it
cluster_kmeans <- clusterfit::kmeans(k = 2:4, init = "random", fuzzy = FALSE)
cluster_kmeans$fit(x = as.matrix(mtcars)) # fit a cluster for matrix/dist data
cluster_kmeans$clustering_ # obtain the clustering vector
cluster_kmeans$output_ # Standardized output across implementations

cluster_kmeans$evaluate()  # some internal metrics
cluster_kmeans$plot()      # ggplot object
cluster_kmeans$stability() # stability index

cluster_kmeans$validate(actualIndices) # external metrics

# create a 'pam' object with params. Sets up a grid for k and swap.
cluster_pam <- clusterfit::pam(k = 2:4, swap = c(TRUE, FALSE), fuzzy = TRUE)
cluster_pam$fit(x = dist(as.matrix(mtcars))) 
cluster_pam$clustering_
cluster_pam$output_

cluster_pam$evaluate()
cluster_kmeans$plot()

cluster_kmeans$validate(actualIndices)

# evaluate some internal metrics and rank aggregate (motivation 'optCluster')
chooseClustering(cluster_kmeans, cluster_pam)
# Get an ensemble of clusterings: think more about the design
ensembleClustering(cluster_kmeans, cluster_pam)
```
