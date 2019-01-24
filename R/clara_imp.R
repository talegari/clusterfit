#' @name CLARA
#' @title CLARA clustering
#' @description Implements CLARA clustering algorithm using
#'   \code{\link[cluster]{pam}}
#' @author Srikanth Komala Sheshachala (sri.teach@@gmail.com)
#' @details CLARA implementation:
#'
#'   \itemize{
#'
#'   \item PAM clustering is computed on multiple random samples of
#'   observations.
#'
#'   \item For a given clustering/medoids, cost is defined as the average
#'   dissimilarity/distance between observations(entire dataset) from the
#'   nearest medoid.
#'
#'   \item A clustering/medoids corresponding to the clustering with minimum
#'   cost is chosen.
#'
#'   }
#'
#'   The PAM fitting on multiple subsets is parallelized with \pkg{future}.
#' @param x (numeric matrix or dist) data
#' @param k (positive integer) Number of clusters
#' @param nSamples (positive integer, default: 5) Number of random samples
#' @param sampleFrac (positive fraction, default: 0.1) Fraction of observations
#'   in a sample
#' @param swap (flag, default: FALSE) Whether PAM should involve swap phase
#' @param pamonce (One among 0, 1, 2, default: 0) See pamonce argument in
#'   \code{\link[cluster]{pam}}
#' @return A list with three compoments:
#'
#'   \itemize{
#'
#'   \item clustering: An integer vector indicating the cluster number with
#'   length equal to number of observations
#'
#'   \item medoidsIndex: An integer vector of indices of medoids
#'
#'   \item cost: average dissimilarity/distance between observations(entire
#'   dataset) from the nearest medoid
#'
#'   }
#' @examples
#' set.seed(1)
#' clara(dist(mtcars), k = 4, sampleFrac = 0.4, nSamples = 10)
#' set.seed(2)
#' clara(stats::dist(mtcars, method = "maximum"), k = 4, sampleFrac = 0.4, nSamples = 10)
#' @export
CLARA <- function(x
                 , k
                 , nSamples    = 5
                 , sampleFrac  = 0.1
                 , swap        = FALSE
                 , pamonce     = 0){

  UseMethod("clara", x)
}


#' @export
CLARA.dist <- function(x
                       , k
                       , nSamples    = 5
                       , sampleFrac  = 0.1
                       , swap        = FALSE
                       , pamonce     = 0
                       ){

  # assertions ----
  assertthat::assert_that(inherits(x, "dist"))
  assertthat::assert_that(assertthat::is.count(k))
  assertthat::assert_that(assertthat::is.count(nSamples))
  assertthat::assert_that(sampleFrac > 0 && sampleFrac <= 1)
  assertthat::assert_that(assertthat::is.flag(swap))
  assertthat::assert_that(pamonce %in% c(0, 1, 2))
  largeKMsg <- paste0("Number of observations in sample should be "
                        , "greater than k (number of clusters)")
  assertthat::assert_that(sampleFrac * nrow(x) > k, msg = largeKMsg)

  # core ----
  # create seeds for each sample
  seeds <- sample(1e7, nSamples)

  # function: pam on the subset and compute cost
  pam_subset <- function(seed){

    # sample the distance object
    set.seed(seed)
    sampleIndex <- sample(ceiling(sampleFrac * attr(x, "Size")))
    dSubset     <- disto::dist_subset(x, sampleIndex)

    # fit PAM
    pamObject <- cluster::pam(x              = dSubset
                              , k            = k
                              , stand        = FALSE
                              , cluster.only = FALSE
                              , medoids      = NULL
                              , do.swap      = swap
                              , keep.diss    = FALSE
                              , keep.data    = FALSE
                              , pamonce      = pamonce
                              , trace.lev    = 0
                              )

    # compute average dissimilarity from closest medoids
    distances <- disto::dist_extract(x
                                     , i = pamObject$id.med
                                     , product = "outer"
                                     )

    avgDissim <- mean(matrixStats::colMins(distances))

    # return average dissimilarity and medoids
    list(avgDissim = avgDissim
         , medoids = sampleIndex[pamObject[["id.med"]]]
         )
  }

  # fit PAM over samples
  pamList <- purrr::transpose(future.apply::future_lapply(seeds, pam_subset))

  # pick medoids/sample with minimum cost ----
  pamList[["avgDissim"]] <- unlist(pamList[["avgDissim"]])
  minPos                 <- which.min(pamList[["avgDissim"]])
  minVal                 <- min(pamList[["avgDissim"]])

  # assign the cluster numbers to the entire dataset
  medoidsIndex   <- pamList[["medoids"]][[minPos]]
  distances      <- disto::dist_extract(x
                                        , i = medoidsIndex
                                        , product = "outer"
                                        )
  clustering     <- apply(distances, 2, which.min)

  # return ----
  return(list(clustering      = clustering
              , medoidsIndex  = medoidsIndex
              , cost          = minVal
              )
         )
}

#' @export
CLARA.matrix <- function(x
                         , k
                         , nSamples    = 5
                         , sampleFrac  = 0.1
                         , swap        = FALSE
                         , pamonce     = 0
                         ){

    # assertions ----
    assertthat::assert_that(inherits(x, "matrix") && is.numeric(x))
    assertthat::assert_that(assertthat::is.count(k))
    assertthat::assert_that(assertthat::is.count(nSamples))
    assertthat::assert_that(sampleFrac > 0 && sampleFrac <= 1)
    assertthat::assert_that(assertthat::is.flag(swap))
    assertthat::assert_that(pamonce %in% c(0, 1, 2))
    largeKMsg <- paste0("Number of observations in sample should be "
                        , "greater than k (number of clusters)")
    assertthat::assert_that(sampleFrac * nrow(x) > k, msg = largeKMsg)

    # core ----
    # create seeds for each sample
    seeds <- sample(1e7, nSamples)

    # function: pam on the subset and compute cost
    pam_subset <- function(seed){

      # sample the distance object
      set.seed(seed)
      sampleIndex <- sample(ceiling(sampleFrac * nrow(x)))
      xSubset     <- x[sampleIndex, ]

      # fit PAM
      pamObject <- cluster::pam(x              = xSubset
                                , k            = k
                                , stand        = FALSE
                                , cluster.only = FALSE
                                , medoids      = NULL
                                , do.swap      = swap
                                , keep.diss    = FALSE
                                , keep.data    = FALSE
                                , pamonce      = pamonce
                                , trace.lev    = 0
                                )

      # compute average dissimilarity from closest medoids
      distances <- proxy::dist(x[pamObject$id.med, ], x)

      avgDissim <- mean(matrixStats::colMins(distances))

      # return average dissimilarity and medoids
      list(avgDissim = avgDissim
           , medoids = sampleIndex[pamObject[["id.med"]]]
           )
    }

    # fit PAM over samples
    pamList <- purrr::transpose(future.apply::future_lapply(seeds, pam_subset))

    # pick medoids/sample with minimum cost ----
    pamList[["avgDissim"]] <- unlist(pamList[["avgDissim"]])
    minPos                 <- which.min(pamList[["avgDissim"]])
    minVal                 <- min(pamList[["avgDissim"]])

    # assign the cluster numbers to the entire dataset
    medoidsIndex   <- pamList[["medoids"]][[minPos]]
    distances      <- proxy::dist(x[medoidsIndex, ], x)
    clustering     <- apply(distances, 2, which.min)

    # return ----
    return(list(clustering      = clustering
                , medoidsIndex  = medoidsIndex
                , cost          = minVal
                )
           )
}

