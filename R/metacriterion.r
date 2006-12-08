# Compute lc meta criterion
# Compute the lc meta criterion which is a 
# 
# The lc metacriterion computes for each point the 
# proportion of the closest k neighbours in the original
# space that are still among the closest k neighbours in the
# new reduced space.
#
# This can be optional adjusted by subtracting k/(n-1)
# to account for the fact that by chance alone as the neighbourhood
# size increases, the number of points that remain neighbours
# increases.
# 
# @arguments Original (full-dimensional) data set
# @arguments New (lower-dimensioanl) data set
# @arguments Numbers of neighbours to use (numeric vector)
# @arguments Distance metric to use
# @arguments Should adjustment for large k be used?
# @arguments Return pointwise or average lcmeta?
# @arguments Original distance matrix, computed from original data set if not specified
# @keyword multivariate
#X s2 <- localmds(spiral)
#X lcmeta(spiral, s2)
#X lcmeta(spiral, s2, 1:40)
#X lcmeta(spiral, s2, adjust=FALSE)
#X lcmeta(spiral, s2, pointwise=TRUE, adjust=FALSE)
#X lcmeta(spiral, s2, pointwise=TRUE)
lcmeta <- function(original, new, k=2:4, metric="euclidean", adjust=TRUE, pointwise=FALSE, do = as.matrix(dist(original, method=metric))) {
	n <- nrow(new)
	dn <- as.matrix(dist(new, method=metric))
	
	do.nn <- t(apply(do, 1, order)[-1, , drop=FALSE])
	dn.nn <- t(apply(dn, 1, order)[-1, , drop=FALSE])
	
	m <- do.call(cbind, lapply(k, function(i) {
		neighbours <- cbind(do.nn[, 1:i], dn.nn[, 1:i])
		apply(neighbours, 1, function(x) sum(duplicated(x)))
	})) 
	
	if (pointwise) {
		colnames(m) <- k
		m <- m / k 
	} else {
		names(m) <- k
		m <- colMeans(m) / k
	}
	if (adjust) m <- m - k / (n-1)
	
	m
}

# Find closest neighbour
# Used to check that lcmeta is working correctly
# 
# @keyword internal
#X s2 <- localmds(spiral)
#X man <- mean(closest(spiral) == closest(s2))
#X lc1 <- lcmeta(spiral, s2, 1, adjust=FALSE)
#X stopifnot(all.equal(man, lc1))
closest <- function(x) {
	apply(as.matrix(dist(x)), 1, function(x) order(x)[2])
}
