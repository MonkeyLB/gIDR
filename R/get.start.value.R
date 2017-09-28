#' use a continuous version to find starting values of mu, sigma, rho and n.missing
#'
#' @param x replicate one
#' @param y replicate two
#' @return starting values generated
#'
get.start.value <- function(x, y){

  n.01 <- sum(x==0 & y !=0)
  n.10 <- sum(x!=0 & y ==0)
  n.11 <- sum(x!=0 & y !=0)
  n.missing.s <- (n.01+n.10)/2
  p.s <- (n.11-n.missing.s)/(n.11+3*n.missing.s)
  return(list(p.s=p.s, n.missing.s=n.missing.s, n.01=n.01, n.10=n.10, n.11=n.11))
}
