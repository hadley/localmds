# Local MDS
# 
# @keyword internal
# @seealso \code{\link{localmds.animate}}, \code{\link{localmds.multistart}}
localmds <- function(data, ...) UseMethod("localmds", data)

# Local MDS on a distance matrix.
# 
# Information about how localmds works here.
# 
# @arguments maximum number of interations
# @arguments starting configuration, defaults to random locations.  See also \code{\link{pc}} and \code{\link{localmds.multistart}}
# @arguments dimensionality of reduced space
# @arguments
# @arguments
# @arguments
# @arguments amount of static repulsive force, see also \code{\link{tausearch}}
# @arguments maximum number of iterations to perform
# @arguments metric threshold to use, if not specified neighbour will be based on \code{k} nearest neighbours, see \code{\link{threshold}} for more details
# @arguments should the computed (local) distance matrix be saved as an attribute on the result?
# @arguments number of neighbours to use to determine locality
# @keyword multivariate
#X localmds(eurodist)
#X ggobi(localmds(eurodist, savedist=TRUE))
localmds.dist <- function(data, initial=NULL, d=2, lambda=1, mu=1, nu=1, tau=0.1, maxiter=1000, threshold=NULL, savedist=FALSE,k=4, ...) {
  data <- as.matrix(data)
  data <- threshold(data, k=k, thresholds=threshold) 

  if (is.null(initial)) initial <- matrix(rnorm(d * nrow(data)),ncol=d)

  res <- rcboxcox(initial, data, lam=lambda, mu=mu, nu=nu, tau=tau, niter=maxiter)
  rownames(res) <- rownames(data)
  if (savedist) attr(res, "dist") <- data
  if (!missing(threshold)) attr(res,"threshold") <- threshold
  attr(res, "k") <- k
  class(res) <- c("localmds", class(res))
  
  res
}

# Local MDS for a graph
# Use localmds as a graph layout technique
# 
# @keyword internal
#X# localmds.graph(business)
localmds.graph <- function(data, ...) {
  if (!require("graph")) stop("graph package required.")
	data <- ftM2adjM(t(edgeMatrix(data)))
	localmds(data, ...)
}

# Local MDS on a data frame
# Provides a convenient method for converting a data.frame to a distance matrix and then applying localmds.
# 
# @arguments data frame
# @arguments dimensionality of projection
# @arguments distance metric to use, see \code{\link{dist}} for options
# @arguments other arguments passed to \code{\link{localmds.dist}}
# @keyword multivariate
#X plot(cmdscale(dist(swiss)))
#X plot(localmds(swiss, k=5))
localmds.data.frame <- function(data, d=2, metric="euclidean", ...) {
	localmds.matrix(data.matrix(data), d=d, metric=metric, ...)
}

# Local MDS on a matrix
# Provides a convenient method for converting a matrix to a distance matrix and then applying localmds.
# 
# @arguments matrix of data points
# @arguments dimensionality of projection
# @arguments distance metric to use, see \code{\link{dist}} for options
# @arguments other arguments passed to \code{\link{localmds.dist}}
# @keyword multivariate
#X s2 <- localmds(spiral)
#X plot(s2)
#X ggobi(s2)
#X s2a <- localmds(spiral, initial=s2, tau=1e-2)
localmds.matrix <- function(data, d=2, metric="euclidean", ...) {
  D0 <- as.matrix(stats::dist(data, method=metric))
  localmds.dist(D0, d=d, ...)
}