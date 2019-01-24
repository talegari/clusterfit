# PAM object and methods
library("magrittr")

pam <- R6::R6Class(classname = "clusterfit"
  , public = list(

    method              = "pam"
    , object_             = NULL
    , clustering_       = NULL
    , clustering_fuzzy_ = NULL

    , # end of public list

    # initialize or new
    initialize = function(k = 2:5, swap = FALSE, ...){

      k    <- unique(k)
      swap <- unique(swap)

      private$nonGridArgs <- list(...)
      if(any(c("k", "do.swap") %in% names(private$nonGridArgsnonGridArgs))){
        stop("Arguments 'k', 'do.swap' should be passed with appropriate names and not as the part of '...'")
      }

      assertthat::assert_that(assertthat:::is.integerish(k))
      assertthat::assert_that(all(swap %in% c(TRUE, FALSE)))

      # init default values for not set non grid args
      if(!("metric" %in% names(private$nonGridArgs))){
        private$nonGridArgs$metric <- "euclidean"
      }

      if(!("keep.diss" %in% names(private$nonGridArgs))){
        private$nonGridArgs$keep.diss <- FALSE
      }

      if(!("keep.diss" %in% names(private$nonGridArgs))){
        private$nonGridArgs$keep.data <- FALSE
      }

      listForGrid <- list(k = k, swap = swap)

      private$argGrid   <- purrr::cross(listForGrid)
      private$argGridDF <- purrr::cross_df(listForGrid)
    }
    ,
    # fit the model for a data
    fit = function(x){

      if(inherits(x, "matrix")){
        private$inputType <- "matrix"
        private$x         <- x
      } else {
        assertthat::assert_that(inherits(x, "dist"))
        private$inputType <- "dist"
        private$xDist     <- x
      }

      # function to fit the model
      fitter <- function(input){
        if(private$inputType == "matrix"){
          do.call(cluster::pam
                 , c(list(x              = private$x
                          , k            = input$k
                          , medoids      = input$medoids
                          , do.swap      = input$swap
                          )
                     , private$nonGridArgs
                     )
                 )
        } else { # dist case
          do.call(cluster::pam
                   , c(list(x              = private$xDist
                            , k            = input$k
                            , medoids      = input$medoids
                            , do.swap      = input$swap
                            )
                       , private$nonGridArgs
                       )
                   )
        }

    }

      # fit across the arg grid
      self$object_        <- private$argGridDF
      self$object_$object <- future.apply::future_lapply(private$argGrid
                                                         , fitter
                                                         )

      # get clustering
      self$clustering_            <- private$argGridDF
      self$clustering_$clustering <- purrr::map(self$object_[["object"]]
                                                , function(x) x[["clustering"]]
                                                )

      # function to assign clustering_fuzzy using a pam object
      clustering_fuzzy <- function(model){
        if(private$inputType == "dist"){
          distances <- disto::dist_extract(object    = private$xDist
                                           , i       = model$id.med
                                           , product = "outer"
                                           )
        } else {
          distances <- proxy::dist(private$x[model$id.med, ], private$x)
        }

        res <- t(apply(distances, 2, function(acol) acol/sum(acol)))
        return(res)
      }

      # assign fuzzy clustering
      self$clustering_fuzzy_                  <- private$argGridDF
      self$clustering_fuzzy_$clustering_fuzzy <-
        purrr::map(self$object_[["object"]], function(x) clustering_fuzzy(x))

  }
    ,
    predict = function(newData, fuzzy = FALSE){

      if(private$inputType == "dist"){
        assertthat::assert_that(inherits(newData, "crossdist"))
        assertthat::assert_that(nrow(newData) == attr(private$xDist, "Size"))
      }
      if(private$inputType == "matrix"){
        assertthat::assert_that(inherits(newData, "matrix"))
        assertthat::assert_that(ncol(newData) == ncol(private$x))
      }
      assertthat::assert_that(assertthat::is.flag(fuzzy))

      get_prediction <- function(newData, model, fuzzy){

        if(inherits(newData, "crossdist")){
          medMat <- newData[model$id.med, ]
        } else {
          medMat <- proxy::dist(private$x[model$id.med, ], newData)
        }

        if(!fuzzy){
            apply(medMat, 2, which.min)
          } else {
            t(apply(medMat, 2, function(x) x/sum(x)))
          }
      }

      pred            <- private$argGridDF
      pred$prediction <- purrr::map(
        self$object_$object
        , function(x) get_prediction(newData, x, fuzzy)
        )
      return(pred)
    }

    ) # end of public
  ,
  private = list(x           = NULL
                 , argGrid   = NULL
                 , argGridDF = NULL
                 , inputType = NULL
                 , xDist     = NULL
                 , nonGridArgs = list()
                 )
)
