# clara object and methods
library("magrittr")

clara <- R6::R6Class(classname = "clusterfit"
  , public = list(

    method              = "clara"
    , object_           = NULL
    , clustering_       = NULL
    , clustering_fuzzy_ = NULL

    , # end of public list

    # initialize or new
    initialize = function(k                = 2:5
                          , nSamples       = 5:10
                          , sampleFrac     = c(0.01, 0.1, 0.2)
                          , ...
                          ){

      k          <- unique(k)
      nSamples   <- unique(nSamples)
      sampleFrac <- unique(sampleFrac)

      private$nonGridArgs <- list(...)
      if(any(c("k", "nSamples", "sampleFrac") %in% names(private$nonGridArgs))){
        stop("Arguments 'k', 'nSamples', 'sampleFrac' should be passed with appropriate names and not as the part of '...'. ")
      }

      assertthat::assert_that(assertthat:::is.integerish(k))
      assertthat::assert_that(all(sapply(nSamples, assertthat::is.count)))
      assertthat::assert_that(all(sapply(sampleFrac, function(x) x > 0 && x <= 1)))

      listForGrid <- list(k = k
                          , nSamples   = nSamples
                          , sampleFrac = sampleFrac
                          )

      private$argGrid   <- purrr::cross(listForGrid)
      private$argGridDF <- purrr::cross_df(listForGrid)
    }
    ,
    # fit the model for a data
    fit = function(x){

      if(inherits(x, "matrix")){
        private$inputType <- "matrix" # flag
        private$x            <- x
      } else {
        assertthat::assert_that(inherits(x, "dist"))
        private$inputType <- "dist"   # flag
        private$xDist        <- x
      }

      # function to fit the model
      fitter <- function(input){

        if(private$inputType == "dist"){
          # inner if for nonGridArgs
          if(length(private$nonGridArgs) > 0){
            do.call(CLARA
                   , c(list(x              = private$xDist
                          , k            = input$k
                          , nSamples     = input$nSamples
                          , sampleFrac   = input$sampleFrac
                          )
                        , private$nonGridArgs
                        )
                   )
          } else {
            do.call(CLARA
               , list(x              = private$xDist
                      , k            = input$k
                      , nSamples     = input$nSamples
                      , sampleFrac   = input$sampleFrac
                      )
               )
          }
        } else { # matrix case
          # inner if for nonGridArgs
          if(length(private$nonGridArgs) > 0){
            do.call(CLARA
                   , c(list(x              = private$x
                            , k            = input$k
                            , nSamples     = input$nSamples
                            , sampleFrac   = input$sampleFrac
                            )
                       , private$nonGridArgs
                       )
                   )
          } else {
            do.call(CLARA
               , list(x              = private$x
                      , k            = input$k
                      , nSamples     = input$nSamples
                      , sampleFrac   = input$sampleFrac
                      )
               )
          }
        }
    }

      # fit across the arg grid
      self$object_        <- private$argGridDF
      self$object_$object <- future.apply::future_lapply(private$argGrid, fitter)

      # get clustering
      self$clustering_            <- private$argGridDF
      self$clustering_$clustering <- purrr::map(self$object_[["object"]]
                                                , function(x) x[["clustering"]]
                                                )

      # function to assign clustering_fuzzy using a pam object
      clustering_fuzzy <- function(model){
        if(private$inputType == "dist"){
          distances <- disto::dist_extract(object    = private$xDist
                                           , i       = model$medoids
                                           , product = "outer"
                                           )
        } else {
          distances <- proxy::dist(private$x[model$medoids, ],private$x)
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
          medMat <- newData[model$medoids, ]
        } else {
          medMat <- proxy::dist(private$x[model$medoids, ], newData)
        }

        if(!fuzzy){
            apply(medMat, 2, which.min)
          } else {
            t(apply(medMat, 2, function(x) x/sum(x)))
          }
      }

      pred            <- private$argGridDF
      pred$prediction <- purrr::map(self$object_$object
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
                 , nonGridArgs = NULL
                 )
)
