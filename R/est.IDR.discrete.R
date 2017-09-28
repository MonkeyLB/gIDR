#' EM algorithm to estimate the copula mixture model for discrete data

#' @param x replicate 1 : a vector of values from original observation
#' @param y replicate 2 : paired with \code{x}, a repeated measure having same length as \code{x}.
#' @param mu a starting value for the mean of the reproducible component.
#' @param sigma a starting value for the standard deviation of the reproducible component.
#' @param rho a starting value for the correlation coefficient of the reproducible component.
#' @param p a starting value for the proportion of reproducible component.
#' @param eps a small number to control convergence. For example, \code{esp=0.0001}.
#' @param n.missing \code{n.missing=0} indicates that we run the program as if there is no missing
#'       values; \code{n.missing=}positive number: the program estimates missing values.
#' @param miss.sym symbol for missing value. For example, \code{miss.sym=0}.
#' @param as.single.loglik an integer n to control the precision of
#'          numberical integration. If n=1, the computation is exact.
#'         If n>1, only integrate the bins with counts more than n,
#'          and treat bins with counts <= n as singletons in the
#'          likelihood computation
#' @param as.single.em similar to \code{as.single.loglik}, except that
#'                     it does the approximation in EM
#' @param common.only If \code{TRUE}, the computation will only use the
#'                    nonmissing entries. In this case, \code{n.missing}
#'                    should be set to 0. If \code{FALSE}, it will use
#'                    all the entries.
#' @param labels the cdf on the left side of the lower boundary of the bins
#' @return
#' \itemize{
#'   \item{para   }{estimated parameters: p, rho, mu, sigma, estimated n.missing.}
#'   \item{loglik  }{can be used for mapping category format back to original
#'                   e.g. \code{x.cat[level.factor[1]] == x[1], y.cat[level.factor[1]]==y[1]}.}
#'   \item{loglik.trace	 }{trajectory of log-likelihood.}
#'   \item{idr.cat  }{a numeric vector of the local idr for each category (i.e. estimated conditional probablility for each observation to belong to the irreproducible component).}
#'   \item{IDR.cat  }{a numerical vector of the expected irreproducible discovery rate for categories that are as irreproducible or more irreproducible than the given categories.}
#'   \item{idr.obs  }{a numeric vector of the local idr for each observation.}
#'   \item{IDR.obs  }{a numerical vector of the expected irreproducible discovery rate for observations that are as irreproducible or more irreproducible than the given observations.}
#' }
#'
#' @details
#' EM to compute the latent structure
#' steps:
#' 1. raw values are first transformed into pseudovalues
#' 2. EM is used to compute the underlining structure, which is a mixture of two normals
#'
#' @examples
#'
#' #load chip_seq data
#' data(chip_seq)
#' x = chip_seq[,1]
#' y = chip_seq[,2]
#'
#' # Initiation
#' mu <- 2.6
#' sigma <- 1.3
#' rho <- 0.8
#' p <- 0.7
#' eps <- 0.001
#' n.missing <- 0
#'
#' # Estimate parameters of mixture model
#' gidr.out <- est.IDR.discrete(x, y, mu, sigma, rho, p, eps, n.missing,
#'                             miss.sym = 0, as.single.loglik = 1,
#'                             as.single.em = 1, common.only=TRUE, labels=NULL)
#'
#' names(gidr.out)

est.IDR.discrete <- function(x, y, mu, sigma, rho, p, eps, n.missing, miss.sym, as.single.loglik, as.single.em, common.only=TRUE, labels=NULL){

  if(n.missing > 0)
    top.missing <- TRUE
  else
    top.missing <- FALSE

  x.ori <- x
  y.ori <- y

  # keep the non-missing ones
  if(common.only){
    which.nomissing <- x != miss.sym & y != miss.sym
    x <- x[which.nomissing]
    y <- y[which.nomissing]
    n.missing <- 0
    top.missing <- FALSE
  }

  if(!is.null(labels))
    xy.cat <- get.cat.xy(x, y, n.missing, miss.sym, top.missing, labels)
  else{
    xy.cat <- get.cat.xy(x, y, n.missing, miss.sym, top.missing)
  }
  xy <- xy.cat$cat.xy

  # initialization
  para <- list()
  para$mu <- mu
  para$sigma <- sigma
  para$rho <- rho
  para$p <- p
  para$m <- xy$m

  j <- 1
  to.run <- T
  loglik.trace <- c()
  loglik.inner.trace <- c()

  z1.pseudo <- get.pseudo.mix.map2(xy$x.cdf, xy$x.cdf.lo, para$mu, para$sigma,  para$p)
  z1.up <- z1.pseudo$x.pseudo
  z1.lo <- z1.pseudo$x.pseudo.lo

  z2.pseudo <- get.pseudo.mix.map2(xy$y.cdf, xy$y.cdf.lo, para$mu, para$sigma,  para$p)
  z2.up <- z2.pseudo$x.pseudo
  z2.lo <- z2.pseudo$x.pseudo.lo

  l.new.outer <- loglik.2binormal.discrete(cbind(z1.lo, z1.up), cbind(z2.lo, z2.up), para$m, para$mu, para$sigma, para$rho, para$p, as.single.loglik, top.missing)
  loglik.trace[1] <- l.new.outer

  while(to.run){

    i <- 1
    j <- j+1


    para.old <- para
    loglik.inner.trace[1] <- l.new.outer

    cat("before EM=", loglik.inner.trace[1], "\n")

    # run EM for a given pseudo-value
    to.run.EM <- TRUE
    while(to.run.EM){
      i <- i+1
      # EM for latent structure
      para <- em.step.2normal.discrete(cbind(z1.lo, z1.up), cbind(z2.lo, z2.up), para$m, para$mu, para$sigma, para$rho, para$p, as.single.em, top.missing)

      loglik.inner.trace[i] <- loglik.2binormal.discrete(cbind(z1.lo, z1.up), cbind(z2.lo, z2.up), para$m, para$mu, para$sigma, para$rho, para$p, as.single.loglik, top.missing)

      if(abs(loglik.inner.trace[i]-loglik.inner.trace[i-1])<=eps | i > 100 | para$p>0.999 | para$p<0.001 | para$rho>0.999 | para$rho < -0.999)
        to.run.EM <- FALSE
    }


    to.run.outer <- T  # always run the first iteration
    k <- 0
    while(to.run.outer & k<10){

      if(!is.null(labels))
        xy.cat <- get.cat.xy(x, y, n.missing, miss.sym, top.missing, labels)
      else{
        xy.cat <- get.cat.xy(x, y, n.missing, miss.sym, top.missing)
      }
      xy <- xy.cat$cat.xy
      z1.pseudo <- get.pseudo.mix.map2(xy$x.cdf, xy$x.cdf.lo, para$mu, para$sigma, para$p)
      z1.up <- z1.pseudo$x.pseudo
      z1.lo <- z1.pseudo$x.pseudo.lo

      z2.pseudo <- get.pseudo.mix.map2(xy$y.cdf, xy$y.cdf.lo, para$mu, para$sigma, para$p)
      z2.up <- z2.pseudo$x.pseudo
      z2.lo <- z2.pseudo$x.pseudo.lo

      l.new.outer <- loglik.2binormal.discrete(cbind(z1.lo, z1.up), cbind(z2.lo, z2.up), para$m, para$mu, para$sigma, para$rho, para$p, as.single.loglik, top.missing)


      if(l.new.outer < loglik.inner.trace[1]-eps){ # if likelihood is small, only go half step
        para$m <- (para$m+para.old$m)/2
        para$mu <- (para$mu+para.old$mu)/2
        para$rho <- (para$rho+para.old$rho)/2
        para$p <- (para$p+para.old$p)/2
        para$sigma <- (para$sigma+para.old$sigma)/2
        k <- k+1
      } else{
        to.run.outer <- FALSE # no need to adjust step
      }

      if(k>= 10){
        para$m <- para.old$m
        para$mu <- para.old$mu
        para$rho <- para.old$rho
        para$p <- para.old$p
        para$sigma <- para.old$sigma
        l.new.outer <- loglik.inner.trace[1]

        if(!is.null(labels))
          xy.cat <- get.cat.xy(x, y, n.missing, miss.sym, top.missing, labels)
        else{
          xy.cat <- get.cat.xy(x, y, n.missing, miss.sym, top.missing)
        }
        #      xy.cat <- get.cat.xy(x, y, para$m[1], miss.sym, top.missing)
        xy <- xy.cat$cat.xy
        z1.pseudo <- get.pseudo.mix.map2(xy$x.cdf, xy$x.cdf.lo, para$mu, para$sigma, para$p)
        z1.up <- z1.pseudo$x.pseudo
        z1.lo <- z1.pseudo$x.pseudo.lo

        z2.pseudo <- get.pseudo.mix.map2(xy$y.cdf, xy$y.cdf.lo, para$mu, para$sigma, para$p)
        z2.up <- z2.pseudo$x.pseudo
        z2.lo <- z2.pseudo$x.pseudo.lo

      }

    }

    loglik.trace[j] <- l.new.outer

    # stop when iteration>50
    if(j > 50 | abs(loglik.trace[j] - loglik.trace[j-1])< eps)
      to.run <- FALSE
  }

  if(top.missing)
    n.est.missing <- para$m[1]
  else
    n.est.missing <- 0

  # calculate idr and IDR for each category
  idr.cat <- 1-para$e.z
  IDR.cat <- get.IDR.discrete(idr.cat, xy$m)

  # calculate idr and IDR for each observation and map to the original obs
  if(common.only){
    idr.obs <- rep(NA, length(x.ori))
    IDR.obs <- rep(NA, length(x.ori))
    idr.obs[which.nomissing] <- idr.cat[xy.cat$level.factor]
    IDR.obs[which.nomissing] <- IDR.cat[xy.cat$level.factor]
  } else {
    num.obs <- (length(xy.cat$level.factor)-length(x.ori)+1):(length(xy.cat$level.factor))
    idr.obs <- idr.cat[xy.cat$level.factor][num.obs]
    IDR.obs <- IDR.cat[xy.cat$level.factor][num.obs]
  }

  return(list(para=list(p=para$p, rho=para$rho, mu=para$mu, sigma=para$sigma, n.missing=n.est.missing),
              loglik=l.new.outer,
              #e.z=para$e.z,
              loglik.trace=loglik.trace,
              #data=xy.cat,
              idr.cat=idr.cat,
              IDR.cat=IDR.cat,
              idr.obs=idr.obs,
              IDR.obs=IDR.obs))
}
