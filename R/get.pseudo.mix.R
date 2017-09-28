#' Convert cdf to pseudo values for one cdf vector
#'
#' @param x.cdf the vector cdf for the observations.
#' @param mu the mean of the reproducible component.
#' @param sigma the standard deviation of the reproducible component.
#' @param p the proportion of reproducible component.
#' @return  the pseduo values for the vector of cdfs.

get.pseudo.mix <- function(x.cdf, mu, sigma, p){

  f.lo <- min(-4, mu-4*sigma)
  f.up <- max(mu+4*sigma, 4)

  a.lo <- p*pnorm(f.lo, mean=mu, sd=sigma) + (1-p)*pnorm(f.lo, mean=0, sd=1)-max(x.cdf)
  a.up <- p*pnorm(f.up, mean=mu, sd=sigma) + (1-p)*pnorm(f.up, mean=0, sd=1)-max(x.cdf)

  x.pseudo <- mapply(get.pseudo.mix.x, x.cdf, mu, sigma, p)
  invisible(x.pseudo)
}
