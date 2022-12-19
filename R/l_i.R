##' Computes the log-likelihood contribution for a subset
##' of data with fixed x
##'
##' @param data Dataframe with columns y, x and cens.
##' @param alpha Intercept
##' @param beta Slope
##' @param sigma2b Random effects variance
##' @param sigma2e Residual variance
##' @param q_rule Quadrature rule
l_i <- function(data, alpha, beta, sigma2b, sigma2e, q_rule){
  n <- nrow(data)
  n_c <- sum(data$cens)
  mu <- alpha + beta * data$x_c
  mu1 <- mu[1:n_c]
  y1 <- data$log_y[1:n_c]
  mu2 <- mu[(n_c + 1):n]
  y2 <- data$log_y[(n_c + 1):n]
  ll <- 0
  if ((n_c > 0) & (n_c < n)){
    ll <- ll + log(pmvnprd(y1,
                           mu = mu1 + sum(y2 -mu2) * sigma2b / (sigma2e + (n - n_c) * sigma2b),
                           c = sigma2b * sigma2e / (sigma2e + (n - n_c) * sigma2b),
                           v = sigma2e + sigma2b * sigma2e / (sigma2e + (n - n_c) * sigma2b),
                           q_rule)
    )
  }
  if (n_c < n){
    ll <- ll + log_dmnorm(y2, mu2, sigma2e, sigma2b)
  }
  if (n_c == n){
    ll <- log(pmvnprd(y1, mu = mu1, c = sigma2b, v = sigma2b + sigma2e, q_rule))
  }
  ll
}
