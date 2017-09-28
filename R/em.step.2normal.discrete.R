# an E-step and an M-step
# z.1, z.2 should be the boundaries of the unique and nonoverlapping categories

#' The E-step and M-step of EM algorithm
#'
#' @param z.1 boundary of the unique and nonoverlapping categories based on pseduo values for the first replicate.
#' @param z.2 boundary of the unique and nonoverlapping categories based on pseduo values for the second replicate.
#' @param m number of observations in each category.
#' @param mu a starting value for the mean of the reproducible component.
#' @param sigma a starting value for the standard deviation of the reproducible component.
#' @param rho a starting value for the correlation coefficient of the reproducible component.
#' @param p a starting value for the proportion of reproducible component.
#' @param count.as.singleton the count being seen as singleton. Default is 1.
#' @param top.missing if there are missing observations, top.missing=T.
#' @return estimated parameters: p, rho, mu, sigma
#'
em.step.2normal.discrete  <- function(z.1, z.2, m, mu, sigma, rho, p, count.as.singleton=1, top.missing=F){

  mu.new <- mu
  sigma.new <- sigma
  rho.new <- rho
  p.new <- p

  # z is a vector
  p1f1 <- function(z){
    return(p*exp(d.binormal(z[1], z[2], mu, sigma, rho)))
  }

  # mixture of gaussian
  f <- function(z){
    return((1-p)*exp(d.binormal(z[1], z[2], 0, 1, 0)) + p*exp(d.binormal(z[1], z[2], mu, sigma, rho)))
  }


  # p1f1
  # z is a matrix with two columns
  p1f1.vec <- function(z){
    return(p*exp(d.binormal(z[,1], z[,2], mu, sigma, rho)))
  }

  # f=p1f1+p0f0
  f.vec <- function(z){
    return((1-p)*exp(d.binormal(z[,1], z[,2], 0, 1, 0)) + p*exp(d.binormal(z[,1], z[,2], mu, sigma, rho)))
  }

  z1.lo <- z.1[,1]
  z1.up <- z.1[,2]
  z2.lo <- z.2[,1]
  z2.up <- z.2[,2]

  z1.mid <- (z1.lo+z1.up)/2
  z2.mid <- (z2.lo+z2.up)/2

  nbin <- length(z1.lo)
  int.p1f1 <- c()
  int.f <- c()
  int.mu.p1f1 <- c()
  int.var.p1f1 <- c()
  int.rho.p1f1 <- c()
  pj <- c() # mass of the jth category

  i.as.single <- which(m<=count.as.singleton | z1.up==z1.lo | z2.up==z2.lo)
  i.nonsingle <- which(m>count.as.singleton & (z1.up !=z1.lo) & (z2.up !=z2.lo))

  n.as.single <- length(i.as.single)
  n.nonsingle <- length(i.nonsingle)

  # If treated as singleton, then this category is approximated using
  # the mid value of the category, f(mid), then integration is treated as the
  # f(mid)*area
  if(n.as.single>0){

    int.p1f1[i.as.single] <- p1f1.vec(cbind(z1.mid[i.as.single], z2.mid[i.as.single]))
    int.f[i.as.single] <- f.vec(cbind(z1.mid[i.as.single], z2.mid[i.as.single]))
    #    int.p1f1[i.as.single] <- p1f1.vec(cbind(z1.mid[i.as.single], z2.mid[i.as.single]))*(z1.up[i.as.single]-z1.lo[i.as.single])*(z2.up[i.as.single]-z2.lo[i.as.single])
    #    int.f[i.as.single] <- f.vec(cbind(z1.mid[i.as.single], z2.mid[i.as.single]))*(z1.up[i.as.single]-z1.lo[i.as.single])*(z2.up[i.as.single]-z2.lo[i.as.single])
    # we actually don't need this
    pj[i.as.single] <- int.f[i.as.single]*(z1.up[i.as.single]-z1.lo[i.as.single])*(z2.up[i.as.single]-z2.lo[i.as.single])

    ######## new change
    #    pj[i.as.single] <- int.f[i.as.single] # approximate mass
  }

  p1f1.adapt <- function(z){
    #      int.p1f1 <- adapt(2, lo=c(z[1], z[2]), up=c(z[3], z[4]), functn=p1f1, eps=1e-6)$value
    int.p1f1 <- adaptIntegrate(p1f1, c(z[1], z[2]), c(z[3], z[4]))$integral
  }

  f.adapt <- function(z){
    #      int.f <- adapt(2, lo=c(z[1], z[2]), up=c(z[3], z[4]), functn=f, eps=1e-6)$value
    int.f <- adaptIntegrate(f, c(z[1], z[2]), c(z[3], z[4]))$integral
  }

  if(n.nonsingle>0){
    int.p1f1[i.nonsingle] <- apply(cbind(z1.lo[i.nonsingle], z2.lo[i.nonsingle], z1.up[i.nonsingle], z2.up[i.nonsingle]), 1, p1f1.adapt)
    int.f[i.nonsingle] <- apply(cbind(z1.lo[i.nonsingle], z2.lo[i.nonsingle], z1.up[i.nonsingle], z2.up[i.nonsingle]), 1, f.adapt)

    pj[i.nonsingle] <- int.f[i.nonsingle]

  }

  ##### change from version 5
  #   if(n.as.single>0){
  #     pj[i.as.single] <- apply(cbind(z1.lo[i.as.single], z2.lo[i.as.single], z1.up[i.as.single], z2.up[i.as.single]), 1, f.adapt)
  #   }
  ### end of change from version 5
  #  if(top.missing){
  #      int.p1f1[1]<- p1f1.adapt(c(z1.lo[1], z2.lo[1], z1.up[1], z2.up[1]))
  #      int.f[1] <- f.adapt(c(z1.lo[1], z2.lo[1], z1.up[1], z2.up[1]))
  #      pj[1] <- int.f[1]

  #  }

  ##
  ## e-step
  ##
  e.z <- int.p1f1/int.f

  # the setup below is based on the first entries of m and int.f are
  # for missing data
  if(top.missing){
    m.obs <- sum(m[-1])

    # change in this version
    pj.obs <- pj[-1] # we don't need this, just for checking how different sum(pj) is from 1
    # m[1] <- m.obs*pj[1]/sum(pj.obs)
    m[1] <- m.obs*pj[1]/(1-pj[1])

    cat("pj[1]=", pj[1], "sum(pj.obs)", sum(pj.obs), "sum(pj)=", sum(pj), "\n\n")
  }


  ##
  ## m-step
  ##

  # mu from m-step
  mu.p1f1 <- function(z){
    return((z[1]+z[2])*p*exp(d.binormal(z[1], z[2], mu, sigma, rho)))
  }

  mu.p1f1.adapt <- function(z){
    #    int.mu.p1f1 <- adapt(2, lo=c(z[1], z[2]), up=c(z[3], z[4]), functn=mu.p1f1, eps=1e-6)$value
    int.mu.p1f1 <- adaptIntegrate(mu.p1f1, c(z[1], z[2]), c(z[3], z[4]))$integral
  }

  if(n.as.single>0){
    int.mu.p1f1[i.as.single] <- int.p1f1[i.as.single]*(z1.mid[i.as.single]+z2.mid[i.as.single])
  }

  if(n.nonsingle>0){
    int.mu.p1f1[i.nonsingle] <- apply(cbind(z1.lo[i.nonsingle], z2.lo[i.nonsingle], z1.up[i.nonsingle], z2.up[i.nonsingle]), 1, mu.p1f1.adapt)
  }


  # compute p and mu
  c1 <- sum(m*int.p1f1/int.f)

  p.new <- c1/sum(m)


  mu.new <- sum(m*int.mu.p1f1/int.f)/c1/2

  # compute rho and sigma^2

  # var from m-step
  var.p1f1 <- function(z){
    return(((z[1]-mu.new)^2+(z[2]-mu.new)^2)*p*exp(d.binormal(z[1], z[2], mu, sigma, rho)))
  }

  # rho from m-step
  rho.p1f1 <- function(z){
    return((z[1]-mu.new)*(z[2]-mu.new)*p*exp(d.binormal(z[1], z[2], mu, sigma, rho)))
  }

  var.p1f1.adapt <- function(z){

    #     int.var.p1f1 <- adapt(2, lo=c(z[1], z[2]), up=c(z[3], z[4]), functn=var.p1f1, eps=1e-6)$value
    int.var.p1f1 <- adaptIntegrate(var.p1f1, c(z[1], z[2]), c(z[3], z[4]))$integral
  }

  rho.p1f1.adapt <- function(z){
    #    int.rho.p1f1 <- adapt(2, lo=c(z[1], z[2]), up=c(z[3], z[4]), functn=rho.p1f1, eps=1e-6)$value
    int.rho.p1f1 <- adaptIntegrate(rho.p1f1, c(z[1], z[2]), c(z[3], z[4]))$integral
  }


  if(n.as.single>0){
    int.var.p1f1[i.as.single] <- int.p1f1[i.as.single]*((z1.mid[i.as.single]-mu.new)^2+(z2.mid[i.as.single]-mu.new)^2)
    int.rho.p1f1[i.as.single] <- int.p1f1[i.as.single]*(z1.mid[i.as.single]-mu.new)*(z2.mid[i.as.single]-mu.new)
  }


  if(n.nonsingle>0){
    int.var.p1f1[i.nonsingle]<- apply(cbind(z1.lo[i.nonsingle], z2.lo[i.nonsingle], z1.up[i.nonsingle], z2.up[i.nonsingle]), 1, var.p1f1.adapt)
    int.rho.p1f1[i.nonsingle]<- apply(cbind(z1.lo[i.nonsingle], z2.lo[i.nonsingle], z1.up[i.nonsingle], z2.up[i.nonsingle]), 1, rho.p1f1.adapt)

  }

  if(top.missing){
    int.var.p1f1[1] <- var.p1f1.adapt(c(z1.lo[1], z2.lo[1], z1.up[1], z2.up[1]))

    int.rho.p1f1[1] <- rho.p1f1.adapt(c(z1.lo[1], z2.lo[1], z1.up[1], z2.up[1]))
  }

  c2 <- sum(m*int.var.p1f1/int.f)

  sigma.new <- sqrt(c2/c1/2)

  rho.new <- 2*sum(m*int.rho.p1f1/int.f)/c2

  return(list(p=p.new, mu=mu.new, sigma=sigma.new, rho=rho.new, e.z=e.z, m=m))
}
