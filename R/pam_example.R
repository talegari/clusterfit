# pam matrix with some optional arg
mtcars           <- scale(mtcars)
rownames(mtcars) <- NULL

test <- pam$new(k = 2:10, swap = c(TRUE, FALSE), stand = TRUE)
test$fit(as.matrix(mtcars)[-(1:10), ])

test$clustering_
test$clustering_$clustering[[1]]
test$clustering_$clustering[[2]]

test$object_
test$object_$object[[1]]

test$clustering_fuzzy_
test$clustering_fuzzy_$clustering_fuzzy[[3]]

test$predict(as.matrix(mtcars[1:10, ]))
test$predict(as.matrix(mtcars[1:10, ]))$prediction[[4]]
test$predict(as.matrix(mtcars[1:10, ]), fuzzy = TRUE)
test$predict(as.matrix(mtcars), fuzzy = TRUE)$prediction[[4]]

# pam with dist
mtcars           <- scale(mtcars)
rownames(mtcars) <- NULL

test <- pam$new(k = 2:10, swap = c(TRUE, FALSE), stand = TRUE)
test$fit(dist(as.matrix(mtcars)[-(1:10), ]))

test$clustering_
test$clustering_$clustering[[1]]
test$clustering_$clustering[[2]]

test$object_
test$object_$object[[1]]

test$clustering_fuzzy_
test$clustering_fuzzy_$clustering_fuzzy[[3]]

res <- test$predict(proxy::dist(as.matrix(mtcars[-(1:10), ])
                                , as.matrix(mtcars[1:10, ])
                                )
                    )
res
res$prediction[[4]]
res2 <- test$predict(proxy::dist(as.matrix(mtcars[-(1:10), ])
                                , as.matrix(mtcars[1:10, ])
                                )
                       , fuzzy = TRUE
                    )
res2
res2$prediction[[4]]
