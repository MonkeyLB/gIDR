#' compute log-likelihood for mixture of two discretized bivariate normals
#'

#' @param z.1.all pseudo values for first normal.
#' @param z.2.all pseduo values for second normal.
#' @param m.all number of observations in each category.
#' @param mu the mean of the reproducible component.
#' @param sigma the standard deviation of the reproducible component.
#' @param rho the correlation coefficient of the reproducible component.
#' @param p the proportion of reproducible component.
#' @param count.as.singleton the count being seen as singleton. Default is 1.
#' @param top.missing if there are missing observations, \code{top.missing=T}.
#' @return log-likehood for mixture of two discretized bivariate normals.
#'
loglik.2binormal.discrete <- function(z.1.all, z.2.all, m.all, mu, sigma, rho, p, count.as.singleton=1, top.missing=F){

  # need first remove m[1], z.1, z.2, and only leave the observed data
  if(top.missing){
    m <- m.all[-1]
    z.1 <- z.1.all[-1,]
    z.2 <- z.2.all[-1,]

    # the boundary of the missing one
    z1.lo.missing <- z.1.all[1,1]
    z1.up.missing <- z.1.all[1,2]
    z2.lo.missing <- z.2.all[1,1]
    z2.up.missing <- z.2.all[1,2]

  } else {
    m <- m.all
    z.1 <- z.1.all
    z.2 <- z.2.all
  }

  z1.lo <- z.1[,1]
  z1.up <- z.1[,2]
  z2.lo <- z.2[,1]
  z2.up <- z.2[,2]

  f <- function(z){
    return((1-p)*exp(d.binormal(z[1], z[2], 0, 1, 0)) + p*exp(d.binormal(z[1], z[2], mu, sigma, rho)))
  }

  f.vec <- function(z){
    return((1-p)*exp(d.binormal(z[,1], z[,2], 0, 1, 0)) + p*exp(d.binormal(z[,1], z[,2], mu, sigma, rho)))
  }

  nbin <- length(m)
  int.f <- rep(NA, nbin)
  l.m <- rep(NA, nbin) # log-liklihood contributed from each term
  m.total <- sum(m) # all observations

  log.pj <- rep(NA, length(m)) # log mass in each category

  # approximate using mid-points
  z1.mid <- (z1.up + z1.lo)/2
  z2.mid <- (z2.up + z2.lo)/2

  # for singletons, only likelihood is computed
  # if m=1

  i.as.single <- which(m<= count.as.singleton | z1.up==z1.lo | z2.up==z2.lo)
  # keep track the number of singletons
  n.as.single <- length(i.as.single)
  if(n.as.single>0){

    # treat as the observation is known
    l.m[i.as.single] <- m[i.as.single]*log(f.vec(cbind(z1.mid[i.as.single], z2.mid[i.as.single])))
    #   (log((1-p)*exp(d.binormal(z1.mid[i.as.single], z2.mid[i.as.single], 0, 1, 0)) + p*exp(d.binormal(z1.mid[i.as.single], z2.mid[i.as.single], mu, sigma, rho))))

    # # approximate categories that are treated as single
    # # (i.e. likelihood is computed by multiplying f by the area)
    #    log.pj[i.as.single] <- (l.m[i.as.single]/m[i.as.single]+log(z1.up[i.as.single]-z1.lo[i.as.single])+log(z2.up[i.as.single]-z2.lo[i.as.single]))
  }

  i.nonsingle <- which(m> count.as.singleton & (z1.up!=z1.lo) & (z2.up!=z2.lo))
  n.nonsingle <- length(i.nonsingle)

  f.adapt <- function(z.m){
    #    int.f.1 <- adapt(2, lo=c(z.m[1], z.m[2]), up=c(z.m[3], z.m[4]), functn=f, eps=1e-6)$value
    int.f.1 <- adaptIntegrate(f, lowerLimit=c(z.m[1], z.m[2]), upperLimit=c(z.m[3], z.m[4]))$integral
  }

  if(n.nonsingle>0){

    log.pj[i.nonsingle] <- log(apply(cbind(z1.lo[i.nonsingle], z2.lo[i.nonsingle], z1.up[i.nonsingle], z2.up[i.nonsingle]), 1, f.adapt))

    l.m[i.nonsingle] <- m[i.nonsingle]*(log.pj[i.nonsingle])
  }

  # compute P(Phi)=1-P(P_missing)
  if(top.missing)
    pj.nonmissing <- 1-apply(cbind(z1.lo.missing, z2.lo.missing, z1.up.missing, z2.up.missing), 1, f.adapt)
  else
    pj.nonmissing <- 1

  # integrate each bin for pj even in midpoint approximation
  #if(n.as.single>0){

  ##### change from version 5
  #  log.pj[i.as.single] <- log(apply(cbind(z1.lo[i.as.single], z2.lo[i.as.single], z1.up[i.as.single], z2.up[i.as.single]), 1, f.adapt))
  ##### end of change
  ### this step is unnecessary, because we don't use log.pj for singletons
  #pj.area <- (z1.up[i.as.single]-z1.lo[i.as.single])*(z2.up[i.as.single]-z2.lo[i.as.single])

  #log.pj[i.as.single[pj.area>0]] <- log(f.vec(cbind(z1.mid[i.as.single[pj.area>0]], z2.mid[i.as.single[ps.area>0]])))+log(pj.area[pj.area>0
  #}

  # change from version 6
  # observed total mass
  # log.pj.obs <- log(sum(exp(log.pj)))
  # end of change from version 6

  ##### change from version 5
  # normalize all pj
  #l.m.all <- sum(m*(log.pj-log.pj.obs)) + sum(l.m[i.as.single]-m[i.as.single]*log.pj[i.as.single])

  ### change from version 6
  #l.m.all <- sum(m*(log.pj-log.pj.obs)) + log_factorial(sum(m)) - sum(log_factorial_vec(m))

  # P0 is actually cancelled (will be useful elsewhere though) and P=1
  if(n.nonsingle>0)
    l.m.all <- sum(l.m) - sum(m)*log(pj.nonmissing) + log_factorial(sum(m)) - sum(log_factorial_vec(m[i.nonsingle])) -log_factorial(sum(m[i.as.single]))
  else
    l.m.all <- sum(l.m) - sum(m)*log(pj.nonmissing) + log_factorial(sum(m))-log_factorial(sum(m[i.as.single]))

  #### end of change from version 6
  ###### end of change
  return(l.m.all)
}
