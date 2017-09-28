#' Convert cdf to pseudo values for one vector of cdfs of categories

#' @param x.cdf the vector with upper boundary of each category.
#' @param x.cdf.lo the vector with lower boundary of each category.
#' @param mu the mean of the reproducible component.
#' @param sigma the standard deviation of the reproducible component.
#' @param p the proportion of reproducible component.
#' @return  the pseduo values for the vector of cdfs of categroies.
#' @details Convert cdf to pseudo values for one vector of cdfs of categories.
#' The computation was done by first using unique values, then mapping back to
#' original vector to get a fast computation.


get.pseudo.mix.map2 <- function(x.cdf, x.cdf.lo, mu, sigma, p){

  x.cdf.unique <- sort(unique(x.cdf))

  x.pseudo.unique <- get.pseudo.mix(x.cdf.unique, mu, sigma, p)

  n.bin <- length(x.cdf.unique)
  lowest <- get.pseudo.mix.x(min(x.cdf.lo), mu, sigma, p)

  x.pseudo.unique.lo <- c(lowest, x.pseudo.unique[-n.bin])

  level.factor <- as.numeric(factor(x.cdf))

  x.pseudo <- x.pseudo.unique[level.factor]
  x.pseudo.lo <- x.pseudo.unique.lo[level.factor]

  invisible(list(x.pseudo=x.pseudo, x.pseudo.lo=x.pseudo.lo))
}
