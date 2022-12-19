##' Vectorized computation of product of normal densities for use in pmvnprd
##'
##' @param y
##' @param z
##' @param rho
fn <- Vectorize(function(y, z, rho){
  prod(pnorm((z - sqrt(rho) * y * sqrt(2))/sqrt(1 - rho)))/sqrt(pi)
}, "y")
##' Computes marginalised likelihood contribution
##'
##' @param x value
##' @param mu mean
##' @param c covariance
##' @param v variance
##' @param q_rule Quadrature rule
pmvnprd <- function(x, mu, c, v, q_rule){
  rho <- c / v
  fastGHQuad::ghQuad(fn, q_rule, z = (x - mu) / sqrt(v), rho = rho)
}
