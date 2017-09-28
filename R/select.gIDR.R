#' Select observations according to gIDR
#'
#' Select observations that exceeding a given IDR level
#'
#' @param x a n by m numeric matrix, where m= num of replicates, n=num of observations. Numerical values representing the significance of the observations, where larger values represent higher significance, for example, -log(p-value). Currently, m=2.
#' @param IDR.x Irreproducibile discovery rate for each entry of x. It is computed from \code{est.IDR.discrete()}.
#' @param IDR.level IDR cutoff, a numerical value between [0, 1]. Typical value is 0.01 or 0.05.
#' @return
#' \itemize{
#'    \item{x:   }{selected reproducible observations given the cutoff.}
#'    \item{n:   }{number of selected reproducible observations.}#'
#'     }
#'
select.gIDR <- function (x, IDR.x, IDR.level) {
    is.selected <- (IDR.x < IDR.level & !is.na(IDR.x))
    x.selected <- x[is.selected, ]
    n.selected <- nrow(x.selected)
    return(list(x = x.selected, n = n.selected))
}
