# Calculate first n principal components
# Convenience function for computing first n principal components
# 
# @arguments data frame or matrix
# @arguments number of pcs to return
# @seealso \code{\link{prcomp}} for real machinery
# @keyword multivariate
#X pc(spiral)
#X
#X f <- frey[,apply(frey, 2,var) > 0.001]
#X localmds(f, initial=pc(f))
pc <- function(data, n=2) {
  predict(prcomp(data, scale.=TRUE))[, 1:n]
}