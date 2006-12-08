# Try different values of tau for best representation
# Better configurations may be obtained starting with large values of tau and getting smaller.  This function explores which values of tau are optimal.
#
# This involves running localmds many times, so may be slow.
# @arguments arguments passed on to localmds
# @arguments values to use for tau.  Automatically sorted in decreasing order
# @arguments neighbour to use for localmds and calculating meta criterion
# @arguments plot results?
# @value matrix of tau and lcmeta values
# @keyword multivariate
#X tausearch(spiral)
#X ks <- 4:10
#X diffk <- do.call(rbind, lapply(ks, function(k) tausearch(spiral, k=k, maxiter=100)))
#X qplot(-log(tau), lcmeta, data=diffk, id=k, colour=factor(k), type="line")
tausearch <- function(..., tau=10^(-(0:7)), k=4, plot=TRUE) {
  tau <- sort(tau, decreasing=TRUE)
  
  n <- length(tau)
  results <- numeric(n)
  config <- NULL
  for(i in 1:n) {
   config <- localmds(..., initial=config, k=k, tau=tau[i], savedist=TRUE)
 
   results[i] <- lcmeta(new=config, k=k, do=attr(config,"dist"))
   #cat(i, ": ", results[i], "\n", sep="")
  }
  res <- data.frame(tau=tau, lcmeta=results, k=k)

  if (plot) qplot(-log(tau), lcmeta, res, type=c("line", "point"))
  res
}