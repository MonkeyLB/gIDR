#' Convert cdf to pseudo values for one cdf value
#'
#' @param x.cdf the cdf for a single observation.
#' @param mu the mean of the reproducible component.
#' @param sigma the standard deviation of the reproducible component.
#' @param p the proportion of reproducible component.
#' @return  the pseduo value for the cdf value.

get.pseudo.mix.x <- function(x.cdf, mu, sigma, p){

  f <- function(x, x.cdf){
    p*pnorm(x, mean=mu, sd=sigma) + (1-p)*pnorm(x, mean=0, sd=1)- x.cdf
  }

  f.lo <- min(-4, mu-4*sigma)
  f.up <- max(mu+4*sigma, 4)

  temp <- uniroot(f, c(f.lo, f.up), tol=0.0001, x.cdf)
  invisible(temp$root)
}
