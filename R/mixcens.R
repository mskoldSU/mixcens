##' Fit censored (from the left) Linear Mixed-Effects Model
##'
##' @param x vector of covariates
##' @param y vector of responses
##' @param cens vector of censoring indicators, if TRUE y is censored from the left
##' @param family lognormal or normal
##' @param use_lmec Should lmec::lmec be used for fitting?
##' @examples
##' fit <- with(lindane_landsort, mixcens(year, value, cens))
##' summary(fit)
##' plot(fit)
##'
##' @import NADA survival
##' @export
##'
mixcens <- function(x, y, cens, family = "lognormal", use_lmec = FALSE){
  call <- match.call()
  family <- match.arg(family, c("lognormal", "normal"))
  C <- diag(4)
  C[1, 2] <- -mean(x)

  q_rule <- fastGHQuad::gaussHermiteData(20)

  if (family == "lognormal"){
    log_y <- log(y)
  } else {
    log_y <- y
    y <- exp(log_y)
  }

  data <- data.frame(y = y,
                     log_y = log_y,
                     x = x,
                     x_c = x - mean(x),
                     cens = cens) |>
    na.omit()

  if (use_lmec == TRUE){
    Z <- matrix(rep_len(1, nrow(data)), ncol = 1)
    X <- cbind(Z, data$x)
    cluster <- as.numeric(as.factor(data$x))
    fit <- lmec::lmec(data$log_y, data$cens, X, Z, cluster, method = "ML", maxstep = 40)
    varFix <- diag(4)
    varFix[1:2, 1:2] <- fit$varFix
    output <- list(call = call,
                   family = family,
                   data = data,
                   beta = c(fit$beta, log(fit$Psi), 2*log(fit$sigma)),
                   beta0 = NULL,
                   varFix = varFix,
                   c = fit$step,
                   loglik = NULL,
                   convergence = NULL
    )
    class(output) <- "mixcens"
    return(invisible(output))
  }

  # # Inits for numerical fitting
  init_fit <- NADA::cenreg(data$y, data$cens, data$x_c) |> summary()
  init <- c(init_fit$table[1, 1], init_fit$table[2,1])
  init_fit0 <- NADA::cenreg(data$y, data$cens, 1) |> summary()
  init0 <- c(init_fit0$table[,1], init_fit0$table[2,1])

  data_lst <- split(data, data$x) |>
    lapply(function(data) data[order(-data$cens), ])

  l <- function(data_lst, pars){
    -sum(sapply(data_lst, function(data) l_i(data, pars[1], pars[2], pars[3], pars[4], q_rule)))
  }

  fitted0 <- optim(par = init0,
                   fn = function(x) l(data = data_lst, pars = c(x[1], 0, exp(x[2]), exp(x[3]))),
                   hessian = TRUE, method = "BFGS")
  fitted <- optim(par = c(fitted0$par[1],0, fitted0$par[2:3]),
                  fn = function(x) l(data = data_lst, pars = c(x[1], x[2], exp(x[3]), exp(x[4]))),
                  hessian = TRUE, method = "BFGS")

  output <- list(call = call,
                 family = family,
                 data = data,
                 beta = as.numeric(C %*% fitted$par),
                 beta0 = as.numeric(C %*% c(fitted0$par[1], 0, fitted0$par[2:3])),
                 varFix = C %*% solve(fitted$hessian) %*% t(C),
                 c = fitted$counts,
                 loglik = 2*(fitted0$value - fitted$value),
                 convergence = fitted$convergence)
  class(output) <- "mixcens"
  invisible(output)
}
