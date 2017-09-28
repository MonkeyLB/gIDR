#' compute IDR for discrete categories
#'
#' @param idr local idr for each category.
#' @param cat.counts the number of observations in each category.
#' @return a numerical vector of the expected irreproducible discovery rate for categories that are as irreproducible or more irreproducible than the given categories.

get.IDR.discrete <- function(idr, cat.counts){

  idr.sum <- idr*cat.counts # idr*counts for each unique value
  o <- order(idr)
  IDR.cat <- c()
  IDR.cat[o] <- cumsum(idr.sum[o])/cumsum(cat.counts[o])

  invisible(IDR.cat)
}
