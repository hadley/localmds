# Animate convergence of local MDS
# This function uses GGobi to produce an animation of local mds progress
#
# Points are coloured according to their \code{\link{lcmeta}} criterion
# value.  White points are poorly represented, dark blue points are well
# represented.  Lines joins local points.
#
# Each update is rescaled to range [0, 1] to ensure that the scale
# remains constant over time.
#
# @arguments All arguments passed on to \code{\link{localmds}}
# @arguments Number of iterations between update
# @argument Number of neighbours to use for \code{\link{lcmeta}}
# @keyword hplot
#X \dontrun{
#X localmds.animate(spiral, threshold=180)
#X localmds.animate(spiral)
#X localmds.animate(spiral, threshold=180, initial=pc(spiral))
#X s2 <- localmds(spiral, threshold=180, tau=1)
#X localmds.animate(spiral, initial=s2, threshold=180)
#X localmds.animate(frey)
#X }
localmds.animate <- function(..., initial=NULL, stepsize=10, k=3) {
  if (!require("rggobi")) stop("You must have the rggobi package installed to use the GGobi tools.")
  
  last <- localmds(..., initial=initial, maxiter=stepsize, savedist=TRUE, k=k)
  dlast <- apply(last, 2, function(x) (x - min(x)) / diff(range(last)))

  do <- attr(last, "dist")

  g <- ggobi(last)
  d <- g[1]
  
  while(TRUE) {
    last <- localmds(..., initial=last, maxiter=stepsize)
    dlast <- scale(last)
    d[[1]] <- dlast[,1]
    d[[2]] <- dlast[,2]

    lc <- as.vector(lcmeta(new=last, k=3, do=do, pointwise=TRUE, adjust=FALSE))
    glyph_colour(d) <- floor(lc * 9) + 1
    while(gtkEventsPending()) gtkMainIteration()
  }
}

# Open a localmds object in GGobi.
#
# If \code{\link{localmds}} has been run with \code{savedist=TRUE}
# this function will automatically add an edge set for all local
# distances.
# 
# The plot is also scaled between 0 and 1 to avoid having to
# change the scale over time.
#
# @arguments localmds object
# @arguments other unused arguments
# @keyword hplot
#X sp <- localmds(spiral, threshold=180, savedist=TRUE)
#X ggobi(sp)
#X ggobi(localmds(spiral, savedist=TRUE))
ggobi.localmds <- function(data, ...) {
  g <- ggobi(scale(data))
  colorscheme(g) <- "Blues 9"
  d <- g[1]
  
  e <- attr(data,"dist")
  if (!is.null(e)) {
    lc <- as.vector(lcmeta(new=data, k=attr(data, "k"), do=e, pointwise=TRUE, adjust=FALSE))
    #d$lc <- lc
    glyph_colour(d) <- floor(lc * 9) + 1

    a <- which(lower.tri(e), arr.ind=TRUE)
    src <- row.names(d)[a[,2]]
    dest <- row.names(d)[a[,1]] 
    lengths <- as.vector(as.dist(e))
    local <- lengths != .myInf

    g$edges <- data.frame(lengths[local])
    edges(g$edges) <- cbind(src[local], dest[local])
    
    d <- displays(g)[[1]]
    edges(d) <- g$edges

  }
 
 invisible(g)
}

# Scale localmds object
# Scales a localmds object to common scale [0, 1] across all dimensions.
#  
# @keyword internal
scale.localmds <- function(x, ...) {
 apply(x, 2, function(col) (col - min(x)) / diff(range(x))) 
}