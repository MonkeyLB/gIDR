#' Convert original observations to emprical cdf
#'
#' @param x.ori replicate 1 : a vector of values from original observation
#' @param y.ori replicate 2 : paired with \code{x.ori}, a repeated measure having same length with \code{x.ori}.
#' @param miss.sym the numerical symbol to represent a missing observation.
#'                   It needs to be smaller than the smallest observed value.
#' @param n.missing the number of observations that miss both \code{x} and \code{y}.
#' @param multi.fac the cdf on the left side of the lower boundary of the bins
#' @param top.missing if missing value is present, \code{top.missing=T}; if not consider missing values, \code{top.missing=F}.
#' @return
#' \itemize{
#'   \item{cat.xy:   }{a data.frame with \code{nrow} = number of unique combinations of \code{(x, y)}.
#'                    It has the following columns:
#'                    \item{x, y: }{original values, only unique combinations to define a category.}
#'                    \item{factor.x, factor.y: }{level of factors, only unique values.}
#'                    \item{x.cdf, y.cdf: }{empirical cdf for each category, this is the upper boundary.}
#'                    \item{x.cdf.lo, y.cdf.lo: }{empirical cdf for each category, lower boundary.
#'                     (i.e. previous upper boundary). To avoid instability at boundary, the lowest \code{x.cdf.lo} (\code{y.cdf.lo}) is
#'                    \code{1/(n.xy+1)} and the highest \code{x.cdf} (\code{y.cdf}) is \code{n.xy/(n.xy+1)}.}
#'                    \item{m: }{number of observations in each category.}
#'                    }
#'    \item{level.factor:  }{can be used for mapping category format back to original
#'                   e.g. \code{x.cat[level.factor[1]] == x[1], y.cat[level.factor[1]]==y[1]}.}
#' }
#' @details Convert \code{(x, y)} from original observations to emprical cdf and unique pairs.
#' @export
#' @examples
#' data(chip_seq)
#' x = chip_seq[,1]
#' y = chip_seq[,2]
#' get.cat.xy(x, y, 0, miss.sym=0, top.missing=TRUE, multi.fac=0.0001)

get.cat.xy <- function(x.ori, y.ori, n.missing, miss.sym=0, top.missing=TRUE, multi.fac=0.0001){

  # add some safeguard to prevent going off limits
  xcdf <- ecdf(x.ori)
  ycdf <- ecdf(y.ori)

  if(sum(xcdf(x.ori)>0.999)>0)
    x.ori[xcdf(x.ori)>0.999] <- min(x.ori[xcdf(x.ori)>0.999])

  if(sum(xcdf(x.ori)<0.0001)>0)
    x.ori[xcdf(x.ori)<0.0001] <- max(x.ori[xcdf(x.ori)<0.0001])

  if(sum(ycdf(y.ori)>0.999)>0)
    y.ori[ycdf(y.ori)>0.999] <- min(y.ori[ycdf(y.ori)>0.999])

  if(sum(ycdf(y.ori)<0.0001)>0)
    y.ori[ycdf(y.ori)<0.0001] <- max(y.ori[ycdf(y.ori)<0.0001])



  # If top.missing is true, consider missing value. Otherwise not
  # In this case, we always keep a cell for missing values.
  # To prevent change of the number of cells when the estimate of
  # n.missing is close to 0, I set it to 1, so that the number of
  # cells always stay the same
  if(top.missing){
    if(n.missing < 1)
      n.missing <- 1

    x.missing <- rep(miss.sym, n.missing)
    y.missing <- rep(miss.sym, n.missing)

    x <- c(x.missing, x.ori)
    y <- c(y.missing, y.ori)
  } else {

    x <- x.ori
    y <- y.ori
  }


  n <- length(x)


  # find the level of x and y
  factor.x <- as.numeric(factor(x))
  factor.y <- as.numeric(factor(y))

  # combined factor
  factor.xy <- factor.y + factor.x*10^(floor(log10(max(factor.y)))+1)
  level.factor.xy <- as.numeric(factor(factor.xy))
  # find number of obs in each combined factor
  counts <- table(level.factor.xy)
  counts.x <- table(factor.x)
  counts.y <- table(factor.y)

  # keep unique combined levels on x and y
  is.duplicated <- duplicated(level.factor.xy)
  factor.x.unique <- factor.x[!is.duplicated]
  factor.y.unique <- factor.y[!is.duplicated]
  x.unique <- x[!is.duplicated]
  y.unique <- y[!is.duplicated]

  n.xy <- length(level.factor.xy)

  # sort by xy to be consistent with counts, i.e. the level of level.factor.xy
  o <- order(level.factor.xy[!is.duplicated])
  x.factor.sorted <- factor.x.unique[o]  # order x by the order of counts
  y.factor.sorted <- factor.y.unique[o]  # order y by the order of counts
  x.sorted <- x.unique[o]
  y.sorted <- y.unique[o]

  # factor for preventing out of bound
  n.cat <- length(x.unique)
  #  multi.fac <- min(multi.fac, 1/(n.cat+1))

  # cdf on marginals, with length=unique values of x and y
  # round to 10^-5 to avoid unnecessary small bins
  x.cdf <- round(cumsum(counts.x)/sum(counts.x),5)*(1-multi.fac*2)+multi.fac
  # First entry needs to be a small number to pretend to be the tail
  # and also needs to be not too small so that its inversion is bounded
  x.cdf.lo <- c(multi.fac, x.cdf[1:(n-1)])

  y.cdf <- round(cumsum(counts.y)/sum(counts.y),5)*(1-multi.fac*2)+multi.fac
  y.cdf.lo <- c(multi.fac, y.cdf[1:(n-1)])

  # cdf on marginals, fill in the vector of xy
  xy.cdf.x <- x.cdf[x.factor.sorted]
  xy.cdf.x.lo <- x.cdf.lo[x.factor.sorted]

  xy.cdf.y <- y.cdf[y.factor.sorted]
  xy.cdf.y.lo <- y.cdf.lo[y.factor.sorted]

  # categories of x and y with unique combined (x,y) entries and marginal cdf
  cat.xy <- data.frame(x=x.sorted, y=y.sorted, factor.x=x.factor.sorted, factor.y=y.factor.sorted, x.cdf=xy.cdf.x, y.cdf=xy.cdf.y, x.cdf.lo=xy.cdf.x.lo, y.cdf.lo=xy.cdf.y.lo, m=as.vector(counts))

  return(list(cat.xy=cat.xy, level.factor=level.factor.xy))

}
