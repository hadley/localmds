# Multistart local MDS
# Try many random starts in order to find the best configuration.
# 
# @arguments All arguments passed on to \code{\link{localmds}}
# @arguments Number of times to try
# @arguments Display running progress?
# @arguments Default criterion function
# @keyword multivariate
#X localmds.multistart(spiral)
#X lc1 <- function (x) -lcmeta(spiral, x, k=2)
#X localmds.multistart(spiral, criterion=lc1)
localmds.multistart <- function(..., times=10, output=TRUE, criterion=energy) {
  best <- localmds(...)
  bestc <- criterion(best)
  
  if (output) .running(1, bestc)
  for (i in 1:(times - 1)) {
    new <- localmds(...)
    newc <- criterion(new)
    better <- newc < bestc
    if (better) {
      best <- new
      bestc <- newc
    }
    if (output) .running(i+1, newc, better)
  }
  
  best 
}

# Status messages
# Output status messages for localmds.multistart
# @keyword internal
.running <- function(i, value, better=FALSE) {
  nb <- if(better) " New best!" else ""
  cat("Iteration: ", i, "  Criterion: ", value, nb, "\n",sep="") 
}