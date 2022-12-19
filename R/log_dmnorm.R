##' Non-censored likelihood parts
##'
##' @param y obs
##' @param mu mean
##' @param sigma2e residual variance
##' @param sigma2b random effect variance
##'
log_dmnorm <- function(y, mu, sigma2e, sigma2b){
  ybar <- mean(y)
  mu <- mean(mu)
  d <- length(y)
  -(sum((y - mu)^2) - d^2*sigma2b*(ybar - mu)^2/(sigma2e + d*sigma2b))/(2*sigma2e) - log((1 + d*sigma2b/sigma2e)*sigma2e^d)/2
}
