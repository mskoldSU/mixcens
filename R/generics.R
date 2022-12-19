##' Generic
##'
##' @param object A mixcens object
##' @method vcov mixcens
##' @export
vcov.mixcens <- function(object, pars = c("effects", "all"), ...){
  pars <- match.arg(pars)
  if (pars == "effects"){
    V <- object$varFix[1:length(object$beta), 1:length(object$beta)]
  }
  if (pars == "all"){
    V <- object$varFix
  }
  V
}

##' Generic
##'
##' @param object A mixcens object
##' @method summary mixcens
##' @export
summary.mixcens <- function(object, ...){
  V <- vcov(object, pars = "all")
  pars <- object$beta
  se <- sqrt(diag(V))
  z_values <- pars / se
  p_values <- 2 * pnorm(abs(z_values), lower.tail = FALSE)
  coefficients <- matrix(c(pars[1:2],
                           se[1:2],
                           z_values[1:2],
                           p_values[1:2]),
                         nrow = 2)
  colnames(coefficients) <- c("Estimate", "Std. Error", "z value", "Pr(>|z|)")
  rownames(coefficients) <- c("(Intercept)", "(Slope)")
  object$coefficients <- coefficients

  class(object) <- "mixcensSummary"
  object
}

##' Generic
##'
##' @param object A mixcens object
##'
##' @method print mixcens
##'
##' @export
print.mixcens <- function(object, ...){
  cat("Call: ")
  print(object$call)
  cat("Coefficients:\n (Intercept) \t Slope\n",
      paste(signif(object$beta[1], 3), "\t\t", signif(object$beta[2], 3)))
}

##' Generic
##'
##' @param object A mixcensSummary object
##'
##' @method print mixcensSummary
##'
##' @export
print.mixcensSummary <- function(object, ...){
  cat("Censored linear mixed model fit by 'mixcens'\n\n")
  cat(paste("Total", sum(object$data$cens), "out of", nrow(object$data), "observations censored.\n\n"))
  printCoefmat(object$coefficients)
  cat(paste("\n\nRandom effect standard deviation: ", signif(sqrt(exp(object$beta[3])), 3)), "\n")
  cat(paste("Residual standard deviation: ", signif(sqrt(exp(object$beta[4])), 3)), "\n")
}

##' Generic
##'
##' @param object A mixcens object
##'
##' @method plot mixcens
##' @import ggplot2
##' @export
##'
plot.mixcens <- function(object, type = "response", ...){
  data <- object$data
  data_nocens <- data[data$cens == FALSE, ]
  data_cens <- data[data$cens == TRUE, ]
  data_max_cens <- data_cens |>
    dplyr::group_by(x) |>
    dplyr::summarise(ymax = suppressWarnings(max(y)), ymin = 0)
  data_pred <- data.frame(x = unique(data$x),
                          y = predict(object, newdata = unique(data$x)))
  ggplot() +
    geom_point(data = data_nocens, aes(x = x, y = y), alpha = .5) +
    geom_linerange(data = data_max_cens, aes(x = x, ymax = ymax, ymin = ymin)) +
    geom_linerange(data = data_cens, aes(y = y, xmin = x - 0.2, xmax = x  +  0.2), alpha = .5) +
    geom_line(data = data_pred, aes(x = x, y = y), size = 1, color = "steelblue") +
    labs(y = "" , x = "") +
    theme_classic() + scale_y_continuous(limits = c(0, NA), expand = expansion(mult = c(0, .1)))
}

##' Generic
##'
##' @param object A mixcens object
##'
##' @method predict mixcens
##'
##' @export
##'
predict.mixcens <- function(object, type = "response", newdata = object$data$x, ...){
  mu <- object$beta[1] + object$beta[2] * newdata
  if ((type == "response") & (object$family == "lognormal"))
    mu <- exp(mu)
  mu
}

##' Generic
##'
##' @param object A mixcensSummary object
##'
##' @method coef mixcensSummary
##'
##' @export
##'
coef.mixcensSummary <- function(object, ...){
  object$coefficients
}

#' @importFrom generics tidy
#' @export
generics::tidy

##' Generic
##'
##' @param object A mixcens object
##'
##' @method tidy mixcens
##'
##' @export
##'
tidy.mixcens <- function(object, ...){
  summary_obj <- summary(object)
  coefs <- coef(summary_obj) |> as.data.frame()
  term <- rownames(coefs)
  rownames(coefs) <- NULL
  colnames(coefs) <- c("estimate", "std.error", "statistic", "p.value")
  cbind(term, coefs)
}

