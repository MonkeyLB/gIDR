#' compute the log factorial
#' @param n integer
#' @return log factorial of the given integer

log_factorial <- function(n){
  return(sum(log(1:n)))
}
