#' compute a vector of log factorial
#'
#' @param n.vec a vector of integers
#' @return a vector of log factorial for given integers

log_factorial_vec <- function(n.vec){
  return(sapply(n.vec, log_factorial))
}
