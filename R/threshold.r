# Compute thresholds for local distances
# This determines the cut off between local and non-local distances.
# 
# @arguments distance matrix
# @arguments number of neighbours to keep
# @arguments metric threshold to apply.  If \code{NULL} will compute threshold based on \code{k}
# @keyword multivariate
threshold <- function(d, k=4, thresholds=NULL) {
  if (is.null(thresholds)) thresholds <- apply(d,2,sort)[k+1,]
  d <- ifelse(d > thresholds, .myInf, d)
  pmin(d, t(d)) # make symmetric
}

